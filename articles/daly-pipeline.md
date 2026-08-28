# DALY Burden Estimation Workflow in anumaan

``` r
library(anumaan)
library(dplyr)
```

## Overview

This vignette covers the **DALY burden estimation pipeline**. It does
not repeat preprocessing — that is in the [Preprocessing Workflow
vignette](https://saketlab.github.io/anumaan/articles/preprocessing-workflow.md).

For Bayesian multivariate probit resistance-profile estimation,
`anumaan` supports both CPU and OpenCL execution backends through the
`compute` argument to
[`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).
This is a computational choice only: it does not alter the statistical
model, validation framework, or DALY estimands.

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

## Resistance Profile Estimation

### Option 1: Convex Optimisation (Pathway 1)

#### 1.1 Aggregate Input: Validate Marginals

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

#### 1.2 Profile Enumeration

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

#### 1.3 Constraint Matrix

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

#### 1.4 Estimate Profile Probabilities

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

#### 1.5 Facility Line-List Pathway

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

#### 1.6 Check Constraint Satisfaction

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

#### 1.7 Bootstrap Uncertainty Intervals

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

### Option 2: Bayesian Multivariate Probit (Pathway 2)

Pathway 2 fits a Bayesian hierarchical multivariate probit model
directly to facility-level, event-level AST data. Unlike Pathway 1, it
can incorporate fixed-effect covariates (age, gender, location, …) and
optional random effects (hospital, admission, …), and estimate
resistance *jointly* across antibiotic classes rather than
reconstructing joint profiles purely from marginal + pairwise
constraints.

This requires the optional `cmdstanr` + CmdStan dependency (see
`DESCRIPTION`) and fits an actual MCMC model, so the code chunks below
are shown with `eval = FALSE` for portability – the printed output is
real, captured from an actual run on the synthetic data below (2 chains,
200 warmup/200 sampling iterations – deliberately small and fast for
illustration; a real analysis uses far more, e.g. 4 chains x 3000
warmup/1000 sampling).

#### 2.1 Prepare Event-Level Class Data

One row per organism-event. Antibiotic class columns hold `0`
(susceptible), `1` (resistant), or `NA` (not tested) – a different
encoding from Pathway 1’s S/R text values, because this is the direct
input to the Stan model.

``` r
set.seed(123)
n_events <- 150
centers  <- c("Hospital A", "Hospital B", "Hospital C")

event_class_data <- tibble::tibble(
  event_id       = paste0("EV", sprintf("%04d", 1:n_events)),
  pathogen       = "Klebsiella pneumoniae",
  center_name    = sample(centers, n_events, replace = TRUE),
  Age_normalised = round(runif(n_events, 0, 90), 1),
  gender         = sample(c("Male", "Female"), n_events, replace = TRUE),
  final_outcome  = sample(c("Discharged", "Died"), n_events, replace = TRUE, prob = c(0.8, 0.2))
)

true_p <- c(Carbapenems = 0.35, Fluoroquinolones = 0.55, Aminoglycosides = 0.30)
for (cls in names(true_p)) {
  vals <- rbinom(n_events, 1, true_p[[cls]])
  vals[sample.int(n_events, size = floor(0.15 * n_events))] <- NA  # some untested
  event_class_data[[cls]] <- vals
}

class_cols <- c("Carbapenems", "Fluoroquinolones", "Aminoglycosides")
```

#### 2.2 Fit the Model

`fixed_effects` is required. `random_effects` may be a legacy character
vector, a named list-of-blocks specification, or
[`list()`](https://rdrr.io/r/base/list.html) for a fixed-effects-only
model. When no random effects are fitted, `profile_group_col` is
required to state the column used for profile aggregation and
validation; it does not make that column a random effect.
`residual_structure = "identity"` treats classes as conditionally
independent given the modelled mean (default, more stable);
`"correlated"` estimates a full residual correlation matrix via an
LKJCholesky prior, but needs adequate pairwise co-testing overlap to be
identifiable (`fit$eligibility_report$pairwise`).

``` r
fit <- fit_bayesian_multivariate_probit(
  event_class_data   = event_class_data,
  class_cols         = class_cols,
  fixed_effects      = c("Age_normalised", "gender"),
  random_effects     = c("center_name"),
  pathogen           = "Klebsiella pneumoniae",
  outcome_col        = "final_outcome",
  residual_structure = "identity",
  prior_config       = list(beta_sd = 1.5, tau_sd = 1.0),
  sampler_config     = list(chains = 2, iter_warmup = 200, iter_sampling = 200,
                             seed = 123, parallel_chains = 2, adapt_delta = 0.9)
)

fit$diagnostics
#> # A tibble: 1 x 36
#>   n_chains iter_warmup iter_sampling n_re_levels n_observed_pairs n_events
#>      <int>       <int>         <int>       <int>            <int>    <int>
#> 1        2         200           200           1              384      149
#> # i 30 more variables: n_classes <int>, max_rhat_structural <dbl>,
#> #   min_ess_bulk_structural <dbl>, min_ess_tail_structural <dbl>,
#> #   n_divergent <int>, converged_structural <lgl>, diagnostic_status <chr>, ...
```

For a fixed-effects-only comparison, retain the same fixed-effect design
and make the grouping purpose explicit without adding a random-effect
block:

``` r
fixed_only_fit <- fit_bayesian_multivariate_probit(
  event_class_data = event_class_data,
  class_cols = class_cols,
  fixed_effects = c("Age_normalised", "gender", "center_name"),
  random_effects = list(),
  profile_group_col = "center_name",
  pathogen = "Klebsiella pneumoniae",
  residual_structure = "identity"
)
```

With these deliberately tiny/fast settings, `converged_structural` comes
back `FALSE` (max R-hat ~1.03, a few divergences) – expected for a
200/200 smoke-test fit, not a sign of a broken model. This mirrors the
`exp_00_smoke_test` convention used for fast iteration before committing
to a full run with production-scale `sampler_config`.

#### 2.3 Convert Posterior Draws to Resistance-Profile Probabilities

``` r
profiles <- compute_event_profile_probabilities(
  fit,
  n_posterior_draws_for_profiles = 200L,
  outcome_col = "final_outcome",
  seed = 123L
)

names(profiles)
#> [1] "event_profiles"  "aggregate_draws"

profiles$event_profiles
#> [compute_event_profile_probabilities] 200 draws | 149 events | 3 hp-pairs | 1 RE level(s)
#> # A tibble: 1,192 x 14
#>   center_name pathogen               event_idx profile_class_set   profile_delta
#>   <chr>       <chr>                      <int> <chr>               <chr>
#> 1 Hospital C  Klebsiella pneumoniae         1   Carbapenems|Fluoro… SSS
#> 2 Hospital C  Klebsiella pneumoniae         1   Carbapenems|Fluoro… RSS
#> # i 1,190 more rows, and 9 more variables: profile_probability <dbl>, ...
```

`event_profiles` is one row per event x profile (up to 2^D rows per
event); `aggregate_draws` carries per-draw
`R_ALL`/`R_KNOWN_OUTCOME`/`R_NF` used for credible intervals in the next
step.

#### 2.4 Aggregate for DALY

``` r
agg <- aggregate_profiles_for_daly(profiles, hospital_col = "center_name", pathogen_col = "pathogen")

agg
#> [aggregate_profiles_for_daly] 24 hospital-pathogen-profile rows
#> # A tibble: 24 x 44
#>   center_name pathogen        profile_class_set    profile_delta R_ALL_mean R_ALL_lower R_ALL_upper ...
#>   <chr>       <chr>           <chr>                <chr>              <dbl>       <dbl>       <dbl>
#> 1 Hospital A  Klebsiella pne… Carbapenems|Fluoro…   RRR
#> # i 23 more rows
```

`agg` is analogous to Pathway 1’s `profiles_out`/`rp_out` – both are the
resistance-profile output this vignette produces; downstream burden
calculation
([`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
onward) is not yet covered here.

#### 2.5 Validate Calibration

Three calibration checks compare the model’s posterior predictions back
against the observed data it was fit on, at increasing levels of joint
complexity: single-class marginals, pairs of classes, and complete
observed profiles.

``` r
val_marginal <- validate_marginal_calibration(fit, n_posterior_draws_for_validation = 200L, seed = 123L)
val_pairwise <- validate_pairwise_calibration(fit, n_posterior_draws_for_validation = 200L, seed = 123L)
val_complete <- validate_complete_profile_calibration(fit, n_posterior_draws_for_validation = 200L,
                                                       seed = 123L, min_complete_events = 5L)

compute_profile_validation_status(
  marginal_tbl          = val_marginal,
  pairwise_tbl          = val_pairwise,
  complete_profile_tbl  = val_complete
)
#> $status
#> [1] "pass"
#> $reasons
#> character(0)
#> $thresholds_used
#> $thresholds_used$max_mean_abs_error_marginal
#> [1] 0.1
#> ...
```

[`compute_profile_validation_status()`](https://saketlab.github.io/anumaan/reference/compute_profile_validation_status.md)
is the single authoritative pass/fail call – it applies default
thresholds (mean absolute error, interval coverage) across all three
tables and returns one combined `status`, rather than requiring the
caller to eyeball three separate tibbles.

#### 2.6 Predictive Checks

Beyond calibration against the fitted data, prior and posterior
predictive checks simulate data *from* the model to check the priors are
reasonable before fitting, and that the fitted model reproduces the
observed data’s statistical structure after fitting.

``` r
prior_pred <- simulate_probit_prior_predictive(fit, n_states = 100L, seed = 123L)
compute_prior_predictive_status(prior_pred)
#> $summary$fraction_probability_lt_0.001
#> [1] 0.478
#> $summary$fraction_all_resistant
#> [1] 0.0857

ppc_stats <- compute_probit_ppc_statistics(
  fit, n_states = 100L, seed = 123L,
  statistics = c("marginal", "resistant_count", "pairwise")
)
compute_posterior_predictive_status(ppc_stats)
#> $status
#> [1] "pass"
#> $family_status$marginal$status
#> [1] "ok"
#> $family_status$pairwise$status
#> [1] "ok"
```

`plot_probit_diagnostics(fit, output_dir, experiment_id, pathogen)`
(needs `bayesplot`) writes trace/rank/pair plots for the monitored
parameters to a PDF, and
[`plot_probit_posterior_predictive_checks()`](https://saketlab.github.io/anumaan/reference/plot_probit_posterior_predictive_checks.md)
visualises the `ppc_stats` tables – both are for interactive review, not
part of the scripted pipeline.

#### 2.7 Unified Dispatcher

[`estimate_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/estimate_resistance_profiles.md)
wraps both pathways behind one interface – switch pathway with `method`:

``` r
# Pathway 1 (convex)
estimate_resistance_profiles(data = marginals, method = "convex", panel_map = panel_map)

# Pathway 2 (Bayesian) -- runs steps 2.1-2.4 above in one call
estimate_resistance_profiles(
  data = event_class_data, method = "bayesian",
  class_cols = class_cols, fixed_effects = c("Age_normalised", "gender"),
  random_effects = c("center_name"), pathogen = "Klebsiella pneumoniae",
  outcome_col = "final_outcome", residual_structure = "identity",
  sampler_config = list(chains = 2, iter_warmup = 200, iter_sampling = 200, seed = 123)
)
#> [estimate_resistance_profiles] Pathway 2 complete.
#> $profiles      # same shape as aggregate_profiles_for_daly() output above
#> $eligibility
#> $diagnostics
#> $fitted_models
#> $config_used
```

Both pathways’ `profiles`/`agg` output are the deliverable of this
vignette – downstream burden calculation (RR assignment, YLL/YLD) is not
yet covered here.

------------------------------------------------------------------------

## Pipeline at a Glance

    Preprocessing output (prep_* pipeline)
              │
              ├─────────────────────────────┬─────────────────────────────────────┐
              ▼ Option 1: Convex (Pathway 1) ▼                Option 2: Bayesian Probit (Pathway 2)
      compute_marginal_resistance()          fit_bayesian_multivariate_probit()
              │  (Step 1: marginals)                 │  (facility event-level AST + covariates)
              ▼                                       ▼
      compute_pairwise_coresistance()        compute_event_profile_probabilities()
              │  (Step 2: pairwise co-R)             │  (posterior draws -> profile probs)
              ▼                                       ▼
      compute_resistance_profiles()          aggregate_profiles_for_daly()
              │  (Step 3: QP -> profiles)            │  (+ validate_*_calibration(),
              │                                       │    predictive checks)
              │    Alternative entry (aggregate       │
              │    marginals): validate_aggregate_    │
              │    inputs() -> estimate_profiles_     │
              │    convex()                           │
              │                                       │
              │         (or, either pathway via one call:
              │          estimate_resistance_profiles(method = "convex" | "bayesian"))
              │                                       │
              └───────────────────┬───────────────────┘
                                   ▼
                      Resistance profile probabilities (this vignette's scope ends here;
                      downstream burden calculation -- RR assignment, YLL/YLD -- not yet covered)

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
#> [1] dplyr_1.2.1        anumaan_0.1.0.9028
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5      jsonlite_2.0.0    compiler_4.6.1    tidyselect_1.2.1 
#>  [5] Rcpp_1.1.2        tidyr_1.3.2       jquerylib_0.1.4   systemfonts_1.3.2
#>  [9] textshaping_1.0.5 yaml_2.3.12       fastmap_1.2.0     lattice_0.22-9   
#> [13] R6_2.6.1          generics_0.1.4    knitr_1.51        tibble_3.3.1     
#> [17] desc_1.4.3        osqp_1.0.0        bslib_0.12.0      pillar_1.11.1    
#> [21] rlang_1.3.0       utf8_1.2.6        cachem_1.1.0      xfun_0.60        
#> [25] quadprog_1.5-8    fs_2.1.0          sass_0.4.10       S7_0.2.2         
#> [29] otel_0.2.0        cli_3.6.6         withr_3.0.3       pkgdown_2.2.1    
#> [33] magrittr_2.0.5    digest_0.6.39     grid_4.6.1        lifecycle_1.0.5  
#> [37] vctrs_0.7.3       evaluate_1.0.5    glue_1.8.1        ragg_1.5.2       
#> [41] purrr_1.2.2       rmarkdown_2.31    tools_4.6.1       pkgconfig_2.0.3  
#> [45] htmltools_0.5.9
```
