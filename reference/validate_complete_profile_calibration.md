# Observed-versus-Model Complete-Profile Validation

For hospital x pathogen panels (from
[`.resolve_profile_class_panel()`](https://saketlab.github.io/anumaan/reference/dot-resolve_profile_class_panel.md),
the same approved-eligibility panel used for DALY profiles) with at
least `min_complete_events` events that have **every** panel class
actually observed (no imputation involved), compares the empirical
complete-profile frequency distribution against the model-implied
profile probability distribution for that same event cohort (product of
per-class \\\Phi(\mu\_{ed})\\, as if the classes had not been observed –
i.e. the model's unconditional prediction, deliberately NOT using
[`compute_event_profile_probabilities()`](https://saketlab.github.io/anumaan/reference/compute_event_profile_probabilities.md)'s
observed-cell-preserving logic, since the point here is to check whether
the model's predictions agree with what was actually measured). Panels
without enough complete events are skipped and the reason is recorded in
the `status` column rather than silently omitted.

## Usage

``` r
validate_complete_profile_calibration(
  fitted_model,
  n_posterior_draws_for_validation = 2000L,
  seed = 123L,
  ci_level = 0.95,
  min_complete_events = 30L
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

- min_complete_events:

  Integer. Minimum number of fully-observed-panel events required to
  evaluate a hospital-pathogen panel. Default `30L`.

## Value

Tibble with one row per hospital x pathogen x profile for evaluated
panels (`status = "evaluated"`: `empirical_count/frequency`,
`model_frequency_mean/lower/upper`, `absolute_error`,
`interval_contains_observed`) and one row per skipped panel (`status`
starts with `"skipped_"`, numeric columns `NA`).
