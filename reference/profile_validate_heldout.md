# Masked-Cell Validation (Genuine Holdout, Refits the Model)

Masks a reproducible, eligibility-preserving subset of genuinely
observed cells, sets them to `NA` in a copy of `Y`, rebuilds the Stan
design on the masked data via
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md),
refits via
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md),
and predicts the held-out cells from THAT refit
(`validation_mode = "refit_masked"`). This is the preferred, more
expensive procedure once computationally affordable; see
[`profile_validate_masked`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md)
for the cheaper in-sample diagnostic.

## Usage

``` r
profile_validate_heldout(
  Y,
  X,
  ids,
  output_dir,
  fraction_to_mask = 0.05,
  max_mask_per_cell = NULL,
  min_tested_after_mask = 10L,
  min_resistant_after_mask = 3L,
  min_susceptible_after_mask = 3L,
  sampler_config = list(),
  n_draws = 500L,
  seed = 123L
)
```

## Arguments

- Y, ids:

  The SAME `Y` matrix and `ids` data frame passed to
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (via
  [`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md))
  – the observed side of every comparison comes from these, not
  re-derived from raw data.

- X:

  The SAME covariate design matrix used for the original fit
  (`nrow(X) == nrow(Y)`) – rebuilt against the masked `Y` via
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md).

- output_dir:

  Character. Passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  for the refit's artifacts.

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

- sampler_config:

  Named list. Passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  for the refit. Default [`list()`](https://rdrr.io/r/base/list.html)
  (that function's own defaults).

- n_draws:

  Integer. Posterior draws used (subsampled if the fit has more).
  Default `500L`.

- seed:

  Integer. Mask selection AND posterior-draw-subsampling seed. Default
  `123L`.

## Value

Named list: `predictions`, `summary` (same shape as
[`profile_validate_masked`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md),
`validation_mode = "refit_masked"`), and `refit_result` (the full return
value of
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
on the masked data, so the refit's own sampler diagnostics can be
checked via
[`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md)
before trusting these scores).
