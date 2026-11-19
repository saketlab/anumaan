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
  strata = NULL,
  holdout_patient_ids = NULL,
  holdout_fraction = 0.1,
  sampler_config = list(),
  n_draws = 500L,
  seed = 123L,
  design_spec = list(),
  min_train_per_stratum = 1L
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

- strata:

  Character/integer vector, length `nrow(Y)` (same row order as
  `patient_id`), or `NULL`. When supplied and `holdout_patient_ids` is
  `NULL`, `holdout_fraction` of DISTINCT patients is sampled WITHIN EACH
  STRATUM separately (e.g. one value per event giving that event's
  hospital) rather than pooled across the whole cohort – avoids a
  hospital-blind draw leaving a small hospital with zero (or
  disproportionately many) held-out patients. A patient must have the
  SAME stratum value across all of their events (checked explicitly;
  errors otherwise, naming the patient). Default `NULL`: pooled
  sampling, the original behaviour, unchanged.

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

- design_spec:

  Named list, any of `class_names`, `max_order`, `beta_sd`,
  `interaction_sd` – forwarded to
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  for the training-only refit. Default
  [`list()`](https://rdrr.io/r/base/list.html): the refit silently uses
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)'s
  own defaults, which can be a different interaction order or prior
  scale from the model actually being validated (e.g. a fit with
  non-default priors). Pass the original fit's own design specification
  here to validate the SAME model, not a default-prior stand-in. Unknown
  field names are an error; the effective specification is returned as
  `effective_design_spec`.

- min_train_per_stratum:

  Integer \>= 1. Minimum number of distinct patients that must remain in
  training per stratum (or overall, when `strata` is `NULL`). Strata
  with no more patients than this are kept entirely in training and
  reported as not evaluable, and `holdout_fraction` is capped so no
  stratum drops below it. Default `1L`.

## Value

Named list: `predictions`, `summary`, `calibration` (same shape as
[`profile_validate_masked`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md)),
`joint_profile_scores` (data frame, one row per held-out event whose
true complete profile is actually known – i.e. every class was genuinely
tested for that event: `joint_log_loss`, `multiclass_brier` scoring the
WHOLE profile as one unit, which cell-level scores in `summary` cannot
do – empty with 0 rows if no held-out event had a complete true
profile), `refit_result` (the training-only fit – check its sampler
diagnostics via
[`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md)
before trusting these scores), `n_holdout_events`, `n_train_events`,
`holdout_patient_ids` (echoed, for provenance).
