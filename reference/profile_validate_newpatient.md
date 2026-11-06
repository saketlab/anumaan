# New-Patient Held-Out Validation (Within Represented Hospitals)

Holds out EVERY event for a selected set of patients entirely from
training (refits on the remaining events via
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)),
then predicts each held-out event using ONLY its covariates – the
`"unconditional"` estimand, i.e. nothing known about that patient's AST
results, matching what is actually available for a genuinely new patient
at deployment – and scores against whichever classes were actually
tested for those events. This is validation WITHIN the hospitals already
represented in training; it is NOT leave-one-hospital-out (this family
does not implement that as a default task – a hospital's fixed effect is
unobserved for an unseen hospital, needing an explicitly different
prediction rule/model).

## Usage

``` r
profile_validate_newpatient(
  Y,
  X,
  ids,
  patient_id,
  output_dir,
  holdout_patient_ids = NULL,
  holdout_fraction = 0.1,
  sampler_config = list(),
  n_draws = 500L,
  seed = 123L
)
```

## Arguments

- Y, X, ids:

  Same as
  [`profile_validate_heldout`](https://saketlab.github.io/anumaan/reference/profile_validate_heldout.md)
  – the FULL (not yet masked/split) data.

- patient_id:

  Character/integer vector, length `nrow(Y)`, in the SAME row order as
  `Y`/`ids` – the patient identifier for every event.

- output_dir:

  Character. Passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  for the training-only refit.

- holdout_patient_ids:

  Character/integer vector or `NULL`. Explicit patient ids to hold out.
  Default `NULL`: sampled via `holdout_fraction` instead.

- holdout_fraction:

  Numeric in (0, 1). Used only when `holdout_patient_ids` is `NULL`:
  this fraction of DISTINCT patients (not events) is sampled for
  holdout. Default `0.1`.

- sampler_config:

  Named list. Passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  for the training-only refit.

- n_draws, seed:

  Passed to
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  for the held-out prediction (also the patient-sampling seed, when
  `holdout_patient_ids` is `NULL`).

## Value

Named list: `predictions`, `summary`, `calibration` (same shape as
[`profile_validate_masked`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md)),
`refit_result` (the training-only fit – check its sampler diagnostics
via
[`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md)
before trusting these scores), `n_holdout_events`, `n_train_events`,
`holdout_patient_ids` (echoed, for provenance).
