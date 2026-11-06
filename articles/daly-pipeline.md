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

1.  **Resistance profile estimation** — converts event-level AST data
    into probability distributions over all 2^D binary resistance
    profiles (S/R per antibiotic class).
2.  **Burden calculation** — applies relative-risk weights to profiles
    and computes YLL, YLD, and total DALY burden.

Resistance profile estimation historically had two pathways: a convex-
optimisation route over pre-computed aggregate marginals (“Pathway 1”),
and the Bayesian multivariate probit route below (“Pathway 2”) fit
directly to facility event-level AST data. Pathway 1 had zero callers
anywhere in production use and was retired; this vignette covers Pathway
2 only.

------------------------------------------------------------------------

## Resistance Profile Estimation

### Bayesian Multivariate Probit

facility-level, event-level AST data. Unlike Pathway 1, it can
incorporate fixed-effect covariates (age, gender, location, …) and
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
centers <- c("Hospital A", "Hospital B", "Hospital C")

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
  vals[sample.int(n_events, size = floor(0.15 * n_events))] <- NA # some untested
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
  event_class_data = event_class_data,
  class_cols = class_cols,
  fixed_effects = c("Age_normalised", "gender"),
  random_effects = c("center_name"),
  pathogen = "Klebsiella pneumoniae",
  outcome_col = "final_outcome",
  residual_structure = "identity",
  prior_config = list(beta_sd = 1.5, tau_sd = 1.0),
  sampler_config = list(
    chains = 2, iter_warmup = 200, iter_sampling = 200,
    seed = 123, parallel_chains = 2, adapt_delta = 0.9
  )
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
val_complete <- validate_complete_profile_calibration(fit,
  n_posterior_draws_for_validation = 200L,
  seed = 123L, min_complete_events = 5L
)

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

#### 2.6 Post-Fit Checklist: What to Run, and Why

Sections 2.2-2.5 introduced each Pathway 2 function individually; this
section ties them into one ordered checklist, since it’s easy to run
only some of them and miss a gap. Every function below reads `fit`, the
object returned by
[`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
– none of them require re-fitting.

**1. Did the sampler converge?** (`fit$diagnostics`, populated at fit
time – no separate call needed). Check `diagnostic_status` (`"pass"` /
`"warning_rhat"` / `"fail_energy"` / `"fail_divergent"` / …),
`converged_structural`, `max_rhat`, `min_ess_bulk`, `ebfmi_min`. This is
a *computational* question – did the four chains actually agree with
each other about where the posterior mass is – and it is answered
independently of everything below. A model can fail here while still
answering later checks well (the reverse is not reassuring: don’t trust
calibration numbers from an unconverged fit).

**2. Estimate the resistance profiles.**
[`compute_event_profile_probabilities()`](https://saketlab.github.io/anumaan/reference/compute_event_profile_probabilities.md)
computes, per event, the posterior probability of every possible
resistance profile – analytically (`Phi(mu_d)` per class) for
`residual_structure = "identity"`, or via Gibbs sampling of the
truncated multivariate normal conditional on that event’s own observed
AST results for `"correlated"`. It never overwrites an observed cell –
only untested classes are imputed, so a fully-observed event’s profile
is a fact about the data, not a claim from the model.
[`aggregate_profiles_for_daly()`](https://saketlab.github.io/anumaan/reference/aggregate_profiles_for_daly.md)
rolls the per-event, per-draw output up to one row per hospital x
pathogen x profile (`R_ALL_mean/lower/upper`, DALY-eligibility flags
gated by `sampler_acceptable` and panel support).

**3. Are the model’s stated probabilities calibrated against reality?**
[`validate_marginal_calibration()`](https://saketlab.github.io/anumaan/reference/validate_marginal_calibration.md),
[`validate_pairwise_calibration()`](https://saketlab.github.io/anumaan/reference/validate_pairwise_calibration.md),
and
[`validate_complete_profile_calibration()`](https://saketlab.github.io/anumaan/reference/validate_complete_profile_calibration.md)
(Section 2.5) compare observed rates to model-implied probabilities at
increasing joint complexity – single class, class pairs, and full
profiles. All three correctly use `Omega` for correlated fits where
relevant (pairwise and complete-profile; marginal probabilities don’t
depend on Omega by construction).
[`mask_and_validate_ast()`](https://saketlab.github.io/anumaan/reference/mask_and_validate_ast.md)
adds a masked-AST holdout check – by default (`refit = FALSE`) this is
an *in-sample* diagnostic (the masked cells were still seen during
fitting, labelled `"in_sample_no_refit"` in its output, not a genuine
holdout), while `refit = TRUE` actually refits without the masked cells
for a real out-of-sample score, at proportionally higher cost.
[`compute_profile_validation_status()`](https://saketlab.github.io/anumaan/reference/compute_profile_validation_status.md)
rolls all four into one `profile_validation_status` verdict.

**Quick reference: which check answers which question**

| Question | Function | Needs `Omega`? |
|----|----|----|
| Did the sampler converge? | `fit$diagnostics` | – |
| What’s the estimated resistance profile? | [`compute_event_profile_probabilities()`](https://saketlab.github.io/anumaan/reference/compute_event_profile_probabilities.md) + [`aggregate_profiles_for_daly()`](https://saketlab.github.io/anumaan/reference/aggregate_profiles_for_daly.md) | correlated fits only |
| Is the model calibrated on single classes? | [`validate_marginal_calibration()`](https://saketlab.github.io/anumaan/reference/validate_marginal_calibration.md) | no |
| Is the model calibrated on class pairs? | [`validate_pairwise_calibration()`](https://saketlab.github.io/anumaan/reference/validate_pairwise_calibration.md) | correlated fits only |
| Is the model calibrated on full profiles? | [`validate_complete_profile_calibration()`](https://saketlab.github.io/anumaan/reference/validate_complete_profile_calibration.md) | correlated fits only |
| Can it predict held-out AST results? | [`mask_and_validate_ast()`](https://saketlab.github.io/anumaan/reference/mask_and_validate_ast.md) | correlated fits only |

A model can score well on some of these and poorly on others at the same
time – e.g. a correlated-residual fit with unresolved sampler
convergence issues can still show much better calibration on class pairs
and full profiles than a well-converged identity-residual fit, because
convergence and calibration are answering different questions (numerical
reliability of the specific parameter estimates vs. structural adequacy
of the model form). Treat “did it converge” and “is it calibrated” as
two separate verdicts, not one combined score.

------------------------------------------------------------------------

## Pipeline at a Glance

    Preprocessing output (prep_* pipeline)
              │
              ▼
      fit_bayesian_multivariate_probit()
              │  (facility event-level AST + covariates)
              ▼
      compute_event_profile_probabilities()
              │  (posterior draws -> profile probs)
              ▼
      aggregate_profiles_for_daly()
              │  (+ validate_*_calibration(), mask_and_validate_ast())
              ▼
      Resistance profile probabilities (this vignette's scope ends here;
      downstream burden calculation -- RR assignment, YLL/YLD -- not yet covered)

------------------------------------------------------------------------

## Session Info

``` r
sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
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
#> [1] dplyr_1.2.1        anumaan_0.1.0.9039
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.3       cli_3.6.6         knitr_1.52        rlang_1.3.0      
#>  [5] xfun_0.62         otel_0.2.0        generics_0.1.4    textshaping_1.0.5
#>  [9] jsonlite_2.0.0    glue_1.8.1        htmltools_0.5.9   ragg_1.5.2       
#> [13] sass_0.4.10       rmarkdown_2.32    tibble_3.3.1      evaluate_1.0.5   
#> [17] jquerylib_0.1.4   fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5  
#> [21] compiler_4.6.2    fs_2.1.0          pkgconfig_2.0.3   systemfonts_1.3.2
#> [25] digest_0.6.39     R6_2.6.1          tidyselect_1.2.1  pillar_1.11.1    
#> [29] magrittr_2.0.5    bslib_0.12.0      tools_4.6.2       pkgdown_2.2.1    
#> [33] cachem_1.1.0      desc_1.4.3
```
