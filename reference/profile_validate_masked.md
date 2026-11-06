# Masked-Cell Validation (In-Sample, No Refit)

Masks a reproducible, eligibility-preserving subset of genuinely
observed cells in `Y` and predicts them from the **already-fitted**
model (`fit_result`) – honestly an in-sample conditional-prediction
diagnostic (`validation_mode = "in_sample_no_refit"`), since the masked
cells still contributed to fitting `beta`/`eta`, not a genuine held-out
test. Useful for a quick miscalibrated-subgroup scan; prefer
[`profile_validate_heldout`](https://saketlab.github.io/anumaan/reference/profile_validate_heldout.md)
for the real holdout.

## Usage

``` r
profile_validate_masked(
  fit_result,
  design,
  Y,
  ids,
  fraction_to_mask = 0.05,
  max_mask_per_cell = NULL,
  min_tested_after_mask = 10L,
  min_resistant_after_mask = 3L,
  min_susceptible_after_mask = 3L,
  n_draws = 500L,
  seed = 123L
)
```

## Arguments

- fit_result:

  List from
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  (reads `$fit`).

- design:

  List from
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (reads `$stan_data`: `X`, `B`, `V`).

- Y, ids:

  The SAME `Y` matrix and `ids` data frame passed to
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (via
  [`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md))
  – the observed side of every comparison comes from these, not
  re-derived from raw data.

- fraction_to_mask:

  Numeric in (0, 1\]. Target fraction of observed cells to mask per
  group x class cell, before the eligibility floor and
  `max_mask_per_cell` cap. Default `0.05`.

- max_mask_per_cell:

  Integer or `NULL`. Optional additional cap.

- min_tested_after_mask, min_resistant_after_mask,
  min_susceptible_after_mask:

  Integer. Panel-eligibility floors enforced during masking. Default
  `10L`, `3L`, `3L`.

- n_draws:

  Integer. Posterior draws used (subsampled if the fit has more).
  Default `500L`.

- seed:

  Integer. Mask selection AND posterior-draw-subsampling seed. Default
  `123L`.

## Value

Named list: `predictions` (one row per masked cell: `group`,
`event_row`, `class`, `true_value`, `predicted_resistance_probability`,
`predicted_class`, `log_score`, `brier_score`, `validation_mode`),
`summary` (one-row aggregate: counts, mean log/Brier score, `log_loss`
(the LogLoss convention, lower-is-better – the NEGATIVE of
`mean_log_score`, both reported so the sign is never ambiguous),
accuracy, `auroc` (`NA`, not 0, with one-class truth),
`average_precision`, `calibration_intercept`/ `calibration_slope` (`NA`
when unestimable)), and `calibration` (binned predicted-vs-observed
calibration curve).
