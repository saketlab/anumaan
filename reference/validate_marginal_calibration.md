# Observed-versus-Model Marginal Resistance Validation

For every hospital x pathogen x antibiotic-class cell with adequate
observed testing support, compares the observed marginal resistance rate
against the model's fitted probability \\\Phi(\mu\_{ed})\\, averaged
over the same tested-event cohort and over posterior draws. This is a
model-validation calculation: it uses the fitted probability for TESTED
cells (unlike
[`compute_event_profile_probabilities()`](https://saketlab.github.io/anumaan/reference/compute_event_profile_probabilities.md),
which never overwrites an observed value). A calibration residual near 0
and interval coverage near the nominal CI level indicate the model
reproduces the observed marginal resistance rates; it does not by itself
validate joint (pairwise/profile) structure.

## Usage

``` r
validate_marginal_calibration(
  fitted_model,
  n_posterior_draws_for_validation = 2000L,
  seed = 123L,
  ci_level = 0.95,
  min_tested = NULL,
  min_resistant = NULL,
  min_susceptible = NULL
)
```

## Arguments

- fitted_model:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- n_posterior_draws_for_validation:

  Integer. Posterior draws used. Default `2000L`.

- seed:

  Integer. Random seed (draw subsampling only). Default `123L`.

- ci_level:

  Numeric. Credible interval coverage. Default `0.95`.

- min_tested, min_resistant, min_susceptible:

  Integer or `NULL`. Eligibility thresholds. When all `NULL` (default),
  eligibility is taken from `fitted_model$eligibility_report$marginal`
  (the same approved rules used for profile panels) rather than
  re-deriving new thresholds. Supplying any of these overrides that
  report and applies the supplied thresholds instead.

## Value

Tibble, one row per eligible hospital x pathogen x class cell:
`n_events`, `n_tested`, `n_resistant`, `n_susceptible`,
`observed_resistance`, `model_resistance_mean/lower/upper`,
`absolute_error`, `calibration_residual`, `interval_contains_observed`.
