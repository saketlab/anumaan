# DALY Burden Estimation Workflow in anumaan

``` r
library(anumaan)
library(dplyr)
```

## Overview

This vignette covers the **DALY burden estimation pipeline**. It does
not repeat preprocessing — that is in the [Preprocessing Workflow
vignette](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md).

The pipeline has two stages:

1.  **Resistance profile estimation** — converts resistance prevalence
    data into probability distributions over all 2^n binary resistance
    profiles (S/R per antibiotic class). Uses convex optimisation (GBD
    eq. 7.5.1.3).
2.  **Burden calculation** — applies relative-risk weights to profiles
    and computes YLL, YLD, and total DALY burden.

Input can come from two sources:

| Source                                                               | Entry point                                                                                                                                                                                                                                                                                                                                          |
|----------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Pre-computed aggregate marginals (GBD, GLASS, national surveillance) | [`validate_aggregate_inputs()`](https://saketlab.github.io/anumaan/reference/validate_aggregate_inputs.md) → [`estimate_profiles_convex()`](https://saketlab.github.io/anumaan/reference/estimate_profiles_convex.md)                                                                                                                                |
| Facility line-list AST data (after preprocessing)                    | [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md) → [`compute_pairwise_coresistance()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_coresistance.md) → [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md) |

------------------------------------------------------------------------

## Stage 1 — Resistance Profile Estimation

### 1.1 Aggregate Input: Validate Marginals

When working from pre-computed marginal resistance rates (e.g. GBD
ST-GPR country estimates or national surveillance summaries), supply
them directly as a tibble and validate before running the QP.

``` r
marginals <- tibble::tibble(
  pathogen         = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 3),
  antibiotic_class = rep(
    c("Carbapenems", "3GC", "Fluoroquinolones"),
    times = 2
  ),
  n_tested    = c(420L, 460L, 390L, 280L, 310L, 265L),
  n_resistant = c(126L, 299L, 195L,  45L, 155L, 133L)
) %>%
  dplyr::mutate(marginal_resistance = n_resistant / n_tested)

marginals
#> # A tibble: 6 × 5
#>   pathogen             antibiotic_class n_tested n_resistant marginal_resistance
#>   <chr>                <chr>               <int>       <int>               <dbl>
#> 1 Klebsiella pneumoni… Carbapenems           420         126               0.3  
#> 2 Klebsiella pneumoni… 3GC                   460         299               0.65 
#> 3 Klebsiella pneumoni… Fluoroquinolones      390         195               0.5  
#> 4 Escherichia coli     Carbapenems           280          45               0.161
#> 5 Escherichia coli     3GC                   310         155               0.5  
#> 6 Escherichia coli     Fluoroquinolones      265         133               0.502
```

``` r
validate_aggregate_inputs(
  marginals,
  pathogen_col    = "pathogen",
  class_col       = "antibiotic_class",
  rate_col        = "marginal_resistance",
  n_tested_col    = "n_tested",
  n_resistant_col = "n_resistant"
)
```

### 1.2 Profile Enumeration

[`enumerate_binary_profiles()`](https://saketlab.github.io/anumaan/reference/enumerate_binary_profiles.md)
generates all 2^n binary combinations for a given ordered class set.
This is the profile space the QP optimises over.

``` r
classes  <- c("Carbapenems", "3GC", "Fluoroquinolones")
profiles <- enumerate_binary_profiles(classes)
profiles
#> # A tibble: 8 × 4
#>   profile_delta Carbapenems `3GC` Fluoroquinolones
#>   <chr>               <int> <int>            <int>
#> 1 SSS                     0     0                0
#> 2 RSS                     1     0                0
#> 3 SRS                     0     1                0
#> 4 RRS                     1     1                0
#> 5 SSR                     0     0                1
#> 6 RSR                     1     0                1
#> 7 SRR                     0     1                1
#> 8 RRR                     1     1                1
```

Each row is one resistance phenotype. `SSS` = pan-susceptible reference;
`RRR` = pan-resistant. The ordered `classes` vector defines the bit
positions and must stay consistent throughout the pipeline.

### 1.3 Constraint Matrix

[`build_constraint_matrix()`](https://saketlab.github.io/anumaan/reference/build_constraint_matrix.md)
constructs the constraint matrix M and target vector v that encode the
marginal and pairwise constraints fed into the QP.

``` r
kp      <- marginals[marginals$pathogen == "Klebsiella pneumoniae", ]
r_marg  <- setNames(kp$marginal_resistance, kp$antibiotic_class)
profiles_enum <- enumerate_binary_profiles(names(r_marg))
cm      <- build_constraint_matrix(profiles_enum, r_marg)

cat("M dimension (constraints x profiles):", dim(cm$M), "\n")
#> M dimension (constraints x profiles): 6 8
cat("Constraint targets (marginals + pairwise independence fallback):\n")
#> Constraint targets (marginals + pairwise independence fallback):
round(cm$v, 4)
#>      Carbapenems              3GC Fluoroquinolones                  
#>            0.300            0.650            0.500            0.195 
#>                                   
#>            0.150            0.325
```

Marginal rows of M: entry = 1 if that class is resistant in that
profile. Pairwise rows: entry = 1 if both classes in the pair are
resistant.

### 1.4 Estimate Profile Probabilities

[`estimate_profiles_convex()`](https://saketlab.github.io/anumaan/reference/estimate_profiles_convex.md)
runs the full pipeline for all pathogens in one call: validates,
enumerates, builds constraints, and solves the QP. The solver prefers
`osqp` (sparse, fast) with `quadprog` as fallback.

``` r
panel_map <- list(
  "Klebsiella pneumoniae"  = c("Carbapenems",
                                "3GC",
                                "Fluoroquinolones"),
  "Escherichia coli"       = c("Carbapenems",
                                "3GC",
                                "Fluoroquinolones")
)

profiles_out <- estimate_profiles_convex(
  marginals    = marginals,
  pairwise     = NULL,
  panel_map    = panel_map,
  lambda       = 1e-8,
  pathogen_col = "pathogen",
  class_col    = "antibiotic_class",
  rate_col     = "marginal_resistance",
  n_tested_col = "n_tested"
)
```

``` r
profiles_out %>%
  dplyr::select(pathogen, profile_delta, profile_probability,
                convergence_flag, max_abs_residual) %>%
  dplyr::filter(profile_probability > 0.005) %>%
  dplyr::arrange(pathogen, dplyr::desc(profile_probability))
#> # A tibble: 16 × 5
#>    pathogen  profile_delta profile_probability convergence_flag max_abs_residual
#>    <chr>     <chr>                       <dbl> <lgl>                       <dbl>
#>  1 Escheric… SSS                         0.125 FALSE                       0.339
#>  2 Escheric… RSS                         0.125 FALSE                       0.339
#>  3 Escheric… SRS                         0.125 FALSE                       0.339
#>  4 Escheric… RRS                         0.125 FALSE                       0.339
#>  5 Escheric… SSR                         0.125 FALSE                       0.339
#>  6 Escheric… RSR                         0.125 FALSE                       0.339
#>  7 Escheric… SRR                         0.125 FALSE                       0.339
#>  8 Escheric… RRR                         0.125 FALSE                       0.339
#>  9 Klebsiel… SSS                         0.125 FALSE                       0.2  
#> 10 Klebsiel… RSS                         0.125 FALSE                       0.2  
#> 11 Klebsiel… SRS                         0.125 FALSE                       0.2  
#> 12 Klebsiel… RRS                         0.125 FALSE                       0.2  
#> 13 Klebsiel… SSR                         0.125 FALSE                       0.2  
#> 14 Klebsiel… RSR                         0.125 FALSE                       0.2  
#> 15 Klebsiel… SRR                         0.125 FALSE                       0.2  
#> 16 Klebsiel… RRR                         0.125 FALSE                       0.2
```

The `profile_class_set` column records the exact ordered class set that
defines the binary code — this must be carried forward to any downstream
DALY attribution to prevent class-order ambiguity.

``` r
profiles_out %>%
  dplyr::distinct(pathogen, profile_class_set, estimator)
#> # A tibble: 2 × 3
#>   pathogen              profile_class_set                estimator
#>   <chr>                 <chr>                            <chr>    
#> 1 Klebsiella pneumoniae 3GC|Carbapenems|Fluoroquinolones convex   
#> 2 Escherichia coli      3GC|Carbapenems|Fluoroquinolones convex
```

### 1.5 Facility Line-List Pathway

When working from isolate-level AST data (after preprocessing), use the
three-step pipeline that feeds the same QP engine.

``` r
set.seed(101)
# Minimal long-format AST data: one row per isolate x antibiotic class
ast_class <- data.frame(
  isolate_id       = rep(paste0("ISO", sprintf("%03d", 1:80)), each = 3),
  organism_name    = rep(
    ifelse(seq_len(80) <= 50, "Klebsiella pneumoniae", "Escherichia coli"),
    each = 3
  ),
  org_group        = "Enterobacterales",
  antibiotic_class = rep(c("Carbapenems", "3GC", "Fluoroquinolones"), times = 80),
  antibiotic_value = sample(
    c("S", "R"), 240, replace = TRUE, prob = c(0.60, 0.40)
  ),
  stringsAsFactors = FALSE
)
```

``` r
# Step 1 -- marginal resistance per pathogen x class
marg_out <- compute_marginal_resistance(
  ast_class,
  pathogen_col         = "organism_name",
  org_group_col        = "org_group",
  isolate_col          = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_value_col = "antibiotic_value",
  min_n_tested         = 10L
)

marg_out$marginal
#> # A tibble: 6 × 6
#>   organism_name         org_group        antibiotic_class n_tested n_resistant
#>   <chr>                 <chr>            <chr>               <int>       <int>
#> 1 Escherichia coli      Enterobacterales Fluoroquinolones       30          14
#> 2 Escherichia coli      Enterobacterales 3GC                    30          12
#> 3 Escherichia coli      Enterobacterales Carbapenems            30          10
#> 4 Klebsiella pneumoniae Enterobacterales Fluoroquinolones       50          27
#> 5 Klebsiella pneumoniae Enterobacterales 3GC                    50          21
#> 6 Klebsiella pneumoniae Enterobacterales Carbapenems            50          17
#> # ℹ 1 more variable: marginal_resistance <dbl>
```

``` r
# Step 2 -- pairwise co-resistance matrices per pathogen
co_out <- compute_pairwise_coresistance(
  marg_out,
  pathogen_col         = "organism_name",
  isolate_col          = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  min_co_tested        = 5L
)

round(co_out[["Klebsiella pneumoniae"]]$prevalence, 3)
#>                   3GC Carbapenems Fluoroquinolones
#> 3GC                NA        0.18             0.24
#> Carbapenems      0.18          NA             0.16
#> Fluoroquinolones 0.24        0.16               NA
```

``` r
# Step 3 -- profile probabilities via QP
rp_out <- compute_resistance_profiles(
  marg_out,
  co_out,
  pathogen_col         = "organism_name",
  antibiotic_class_col = "antibiotic_class",
  exclude_near_zero    = FALSE
)

rp_out[["Klebsiella pneumoniae"]]$profiles %>%
  dplyr::filter(probability > 0.005) %>%
  dplyr::arrange(dplyr::desc(probability))
#>   profile probability 3GC Carbapenems Fluoroquinolones
#> 1     SSS       0.125   0           0                0
#> 2     RSS       0.125   1           0                0
#> 3     SRS       0.125   0           1                0
#> 4     RRS       0.125   1           1                0
#> 5     SSR       0.125   0           0                1
#> 6     RSR       0.125   1           0                1
#> 7     SRR       0.125   0           1                1
#> 8     RRR       0.125   1           1                1
```

``` r
# Constraint residuals: how well does the solution reproduce the inputs?
round(rp_out[["Klebsiella pneumoniae"]]$constraint_residuals, 6)
#>                          marg_3GC                  marg_Carbapenems 
#>                              0.08                              0.16 
#>             marg_Fluoroquinolones              pair_3GC_Carbapenems 
#>                             -0.04                              0.07 
#>         pair_3GC_Fluoroquinolones pair_Carbapenems_Fluoroquinolones 
#>                              0.01                              0.09
```

### 1.6 Check Constraint Satisfaction

[`check_profile_constraints()`](https://saketlab.github.io/anumaan/reference/check_profile_constraints.md)
formally verifies that the estimated probabilities satisfy
non-negativity, sum-to-one, and reproduce the input marginal rates
within tolerance. It accepts the named-list format from
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
directly.

``` r
checks <- check_profile_constraints(
  rp_out,
  marginals    = marg_out$marginal,
  tolerance    = 1e-3,
  pathogen_col = "organism_name",
  class_col    = "antibiotic_class",
  rate_col     = "marginal_resistance"
)

checks %>%
  dplyr::select(pathogen, constraint_type, constraint_name,
                target, reconstructed, abs_residual, pass)
#> # A tibble: 22 × 7
#>    pathogen    constraint_type constraint_name target reconstructed abs_residual
#>    <chr>       <chr>           <chr>            <dbl>         <dbl>        <dbl>
#>  1 Escherichi… nonneg          min_probability NA             0.125      NA     
#>  2 Escherichi… sum_to_one      sum_probability  1             1           0     
#>  3 Escherichi… marginal        marg_3GC         0.4           0.5         0.1   
#>  4 Escherichi… marginal        marg_Carbapene…  0.333         0.5         0.167 
#>  5 Escherichi… marginal        marg_Fluoroqui…  0.467         0.5         0.0333
#>  6 Escherichi… pairwise        pair_3GC_Carba…  0.167         0.25        0.0833
#>  7 Escherichi… pairwise        pair_3GC_Fluor…  0.133         0.25        0.117 
#>  8 Escherichi… pairwise        pair_Carbapene…  0.2           0.25        0.05  
#>  9 Escherichi… marginal        marg_Fluoroqui…  0.467         0.5         0.0333
#> 10 Escherichi… marginal        marg_3GC         0.4           0.5         0.1   
#> # ℹ 12 more rows
#> # ℹ 1 more variable: pass <lgl>
```

### 1.7 Bootstrap Uncertainty Intervals

[`bootstrap_profiles_convex()`](https://saketlab.github.io/anumaan/reference/bootstrap_profiles_convex.md)
resamples resistant counts from a Binomial distribution and refits the
QP B times, returning percentile confidence intervals for each profile
probability.

``` r
boot <- bootstrap_profiles_convex(
  marginals       = marginals,
  B               = 300L,
  seed            = 42L,
  alpha           = 0.05,
  pathogen_col    = "pathogen",
  class_col       = "antibiotic_class",
  n_tested_col    = "n_tested",
  n_resistant_col = "n_resistant"
)

boot[["Klebsiella pneumoniae"]] %>%
  dplyr::filter(probability_mean > 0.005) %>%
  dplyr::arrange(dplyr::desc(probability_mean))
#> # A tibble: 0 × 7
#> # ℹ 7 variables: profile <chr>, probability_mean <dbl>,
#> #   probability_median <dbl>, lower <dbl>, upper <dbl>,
#> #   n_replicates_converged <int>, convergence_rate <dbl>
```

------------------------------------------------------------------------

## Stage 2 — Assigning Relative Risk to Profiles

[`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
applies the GBD **max rule**: the profile-level LOS relative risk is the
maximum RR across all resistant classes in that profile.

``` r
# Synthetic LOS relative risk estimates (from daly_fit_los_rr() in practice)
# pathogen column must match the keys in rp_out (organism_name column values)
rr_table <- tibble::tibble(
  organism_name    = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 3),
  antibiotic_class = rep(c("Carbapenems", "3GC", "Fluoroquinolones"), times = 2),
  RR_LOS    = c(2.1, 1.6, 1.4,  1.8, 1.5, 1.3),
  CI_lower  = c(1.7, 1.3, 1.1,  1.4, 1.2, 1.0),
  CI_upper  = c(2.6, 2.0, 1.8,  2.3, 1.9, 1.7)
)

rr_table
#> # A tibble: 6 × 5
#>   organism_name         antibiotic_class RR_LOS CI_lower CI_upper
#>   <chr>                 <chr>             <dbl>    <dbl>    <dbl>
#> 1 Klebsiella pneumoniae Carbapenems         2.1      1.7      2.6
#> 2 Klebsiella pneumoniae 3GC                 1.6      1.3      2  
#> 3 Klebsiella pneumoniae Fluoroquinolones    1.4      1.1      1.8
#> 4 Escherichia coli      Carbapenems         1.8      1.4      2.3
#> 5 Escherichia coli      3GC                 1.5      1.2      1.9
#> 6 Escherichia coli      Fluoroquinolones    1.3      1        1.7
```

``` r
profiles_rr <- daly_assign_rr_to_profiles(
  profiles_output = rp_out,
  rr_table        = rr_table,
  pathogen_col    = "organism_name",
  class_col       = "antibiotic_class",
  rr_col          = "RR_LOS",
  fallback_rr     = 1
)

profiles_rr[["Klebsiella pneumoniae"]] %>%
  dplyr::select(profile, probability, RR_LOS_profile, dominant_class) %>%
  dplyr::arrange(dplyr::desc(probability))
#>   profile probability RR_LOS_profile   dominant_class
#> 1     SSS       0.125            1.0  all_susceptible
#> 2     RSS       0.125            1.6              3GC
#> 3     SRS       0.125            2.1      Carbapenems
#> 4     RRS       0.125            2.1      Carbapenems
#> 5     SSR       0.125            1.4 Fluoroquinolones
#> 6     RSR       0.125            1.6              3GC
#> 7     SRR       0.125            2.1      Carbapenems
#> 8     RRR       0.125            2.1      Carbapenems
```

### Filter to Classes with RR Estimates

[`daly_filter_profiles_to_rr_classes()`](https://saketlab.github.io/anumaan/reference/daly_filter_profiles_to_rr_classes.md)
drops profiles whose resistant classes have no RR estimate (cannot
contribute to burden), always keeping the all-susceptible reference
profile. Probabilities are re-normalised after filtering.

``` r
profiles_filtered <- daly_filter_profiles_to_rr_classes(
  profiles_rr,
  rr_table     = rr_table,
  pathogen_col = "organism_name",
  class_col    = "antibiotic_class"
)

nrow(profiles_filtered[["Klebsiella pneumoniae"]])
#> [1] 8
```

------------------------------------------------------------------------

## Stage 3 — Resistance Prevalence for DALY Calculation

### Fatal Resistance Prevalence (R_k)

[`daly_calc_resistance_prevalence_fatal()`](https://saketlab.github.io/anumaan/reference/daly_calc_resistance_prevalence_fatal.md)
computes R_k — the profile-weighted expected proportion of deaths
attributable to resistance. This feeds the YLL calculation.

``` r
fatal_prev <- daly_calc_resistance_prevalence_fatal(
  profiles_with_rr  = profiles_filtered,
  probability_col   = "probability",
  rr_profile_col    = "RR_LOS_profile"
)

# R_k: fatal resistance prevalence per pathogen
sapply(fatal_prev, `[[`, "R_k")
#>      Escherichia coli Klebsiella pneumoniae 
#>              0.920000              0.928571

# E[OR_death]: expected odds ratio of death
sapply(fatal_prev, `[[`, "E_OR_k")
#>      Escherichia coli Klebsiella pneumoniae 
#>                1.5625                1.7500
```

### Select Resistance Class for Attribution

For events with multiple resistant classes,
[`select_resistance_class()`](https://saketlab.github.io/anumaan/reference/select_resistance_class.md)
picks a single class per event using the beta-lactam hierarchy and RR
values — preventing double-counting in DALY attribution.

``` r
# Synthetic event-level class data
event_data <- data.frame(
  event_id         = rep(paste0("EV", 1:5), each = 2),
  antibiotic_class = c("Carbapenems", "3GC",
                       "Carbapenems", "Fluoroquinolones",
                       "3GC", "Fluoroquinolones",
                       "Carbapenems", "Aminoglycosides",
                       "3GC", "Fluoroquinolones"),
  antibiotic_value = "R",
  rr_value         = c(2.1, 1.6,  2.1, 1.4,  1.6, 1.4,  2.1, 1.2,  1.6, 1.4),
  stringsAsFactors = FALSE
)

selected <- select_resistance_class(
  event_data,
  event_col          = "event_id",
  class_col          = "antibiotic_class",
  susceptibility_col = "antibiotic_value",
  rr_col             = "rr_value"
)
#> # A tibble: 2 × 2
#>   antibiotic_class n_events
#>   <chr>               <int>
#> 1 Carbapenems             3
#> 2 3GC                     2

dplyr::select(selected, event_id, antibiotic_class, rr_value, selection_method)
#> # A tibble: 5 × 4
#>   event_id antibiotic_class rr_value selection_method
#>   <chr>    <chr>               <dbl> <chr>           
#> 1 EV1      Carbapenems           2.1 hierarchy_rr    
#> 2 EV2      Carbapenems           2.1 hierarchy_rr    
#> 3 EV3      3GC                   1.6 hierarchy_rr    
#> 4 EV4      Carbapenems           2.1 hierarchy_rr    
#> 5 EV5      3GC                   1.6 hierarchy_rr
```

------------------------------------------------------------------------

## Stage 4 — RR Mapping (Preparing for YLL/YLD)

[`daly_add_rr_mappings()`](https://saketlab.github.io/anumaan/reference/daly_add_rr_mappings.md)
maps organism names and antibiotic classes in your analysis-ready data
to the GBD RR pathogen and drug categories, adding `rr_pathogen` and
`rr_drug` columns.

``` r
sample_events <- data.frame(
  organism_normalized = c("Klebsiella pneumoniae", "Escherichia coli",
                           "Staphylococcus aureus"),
  antibiotic_class    = c("Carbapenems", "Fluoroquinolones", "Glycopeptides"),
  stringsAsFactors    = FALSE
)

mapped <- daly_add_rr_mappings(
  sample_events,
  organism_col = "organism_normalized",
  class_col    = "antibiotic_class"
)

mapped
#>     organism_normalized antibiotic_class           rr_pathogen          rr_drug
#> 1 Klebsiella pneumoniae      Carbapenems Klebsiella pneumoniae      Carbapenems
#> 2      Escherichia coli Fluoroquinolones      Escherichia coli Fluoroquinolones
#> 3 Staphylococcus aureus    Glycopeptides Staphylococcus aureus    Glycopeptides
```

------------------------------------------------------------------------

## Pipeline at a Glance

    Preprocessing output (prep_* pipeline)
              │
              ▼
      compute_marginal_resistance()         # Step 1 — marginals per pathogen x class
              │
              ▼
      compute_pairwise_coresistance()       # Step 2 — pairwise co-resistance matrices
              │
              ▼
      compute_resistance_profiles()         # Step 3 — QP -> profile probabilities
              │
              │    Alternative entry (aggregate marginals):
              │    validate_aggregate_inputs()
              │    estimate_profiles_convex()
              │
              ▼
      daly_assign_rr_to_profiles()          # Assign LOS RR (GBD max rule)
              │
              ▼
      daly_filter_profiles_to_rr_classes()  # Drop unestimable profiles
              │
              ▼
      daly_calc_resistance_prevalence_fatal/nonfatal()
              │
              ├── YLL: daly_calc_yll_associated() / daly_calc_yll_attributable()
              └── YLD: daly_calc_yld_attributable() / daly_calc_paf_los()

------------------------------------------------------------------------

## Session Info

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1        anumaan_0.1.0.9021
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5      jsonlite_2.0.0    compiler_4.6.1    Rcpp_1.1.2       
#>  [5] tidyselect_1.2.1  tidyr_1.3.2       jquerylib_0.1.4   systemfonts_1.3.2
#>  [9] textshaping_1.0.5 yaml_2.3.12       fastmap_1.2.0     lattice_0.22-9   
#> [13] R6_2.6.1          generics_0.1.4    knitr_1.51        htmlwidgets_1.6.4
#> [17] tibble_3.3.1      desc_1.4.3        osqp_1.0.0        bslib_0.11.0     
#> [21] pillar_1.11.1     rlang_1.3.0       utf8_1.2.6        cachem_1.1.0     
#> [25] xfun_0.60         quadprog_1.5-8    fs_2.1.0          sass_0.4.10      
#> [29] S7_0.2.2          otel_0.2.0        cli_3.6.6         withr_3.0.3      
#> [33] pkgdown_2.2.1     magrittr_2.0.5    digest_0.6.39     grid_4.6.1       
#> [37] lifecycle_1.0.5   vctrs_0.7.3       evaluate_1.0.5    glue_1.8.1       
#> [41] ragg_1.5.2        purrr_1.2.2       rmarkdown_2.31    tools_4.6.1      
#> [45] pkgconfig_2.0.3   htmltools_0.5.9
```
