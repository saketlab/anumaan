# Temporal Held-Out Validation

Trains on events before `cutoff_date` and predicts events on or after
it, with NOTHING known about the held-out event's AST results (the
`"unconditional"` estimand) – the SAME prediction and scoring mechanics
as
[`profile_validate_newpatient`](https://saketlab.github.io/anumaan/reference/profile_validate_newpatient.md),
just a date-based split instead of a patient-based one. Reports the
training/ holdout patient-overlap count explicitly when `patient_id` is
supplied, rather than assuming temporal and patient independence.

## Usage

``` r
profile_validate_temporal(
  Y,
  X,
  ids,
  event_date,
  cutoff_date,
  patient_id = NULL,
  output_dir,
  sampler_config = list(),
  n_draws = 500L,
  seed = 123L,
  design_spec = list()
)
```

## Arguments

- Y, X, ids:

  Same as
  [`profile_validate_heldout`](https://saketlab.github.io/anumaan/reference/profile_validate_heldout.md).

- event_date:

  Date vector (or coercible via `as.Date`), length `nrow(Y)`, in the
  SAME row order as `Y`/`ids`.

- cutoff_date:

  Date (or coercible). Events with `event_date < cutoff_date` train;
  `>= cutoff_date` are held out.

- patient_id:

  Character/integer vector or `NULL`. When supplied, the returned
  `patient_overlap` documents how many holdout-period patients also
  appear in the training period – temporal validation is not
  automatically patient-independent. Default `NULL` (no overlap report).

- output_dir:

  Character. Passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  for the training-only refit.

- sampler_config:

  Named list. Passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  for the training-only refit.

- n_draws, seed:

  Passed to
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md).

- design_spec:

  Named list, forwarded to
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  for the training-only refit – see
  [`profile_validate_newpatient`](https://saketlab.github.io/anumaan/reference/profile_validate_newpatient.md)'s
  `design_spec` doc for the exact fields and why this matters.

## Value

Named list: `predictions`, `summary`, `calibration`,
`joint_profile_scores` (see
[`profile_validate_newpatient`](https://saketlab.github.io/anumaan/reference/profile_validate_newpatient.md)
for its exact shape – identical here), `refit_result`,
`n_holdout_events`, `n_train_events`, `cutoff_date` (echoed), and
`patient_overlap` (list: `n_overlapping_patients`,
`fraction_of_holdout_patients` – both `NA` when `patient_id` is not
supplied).
