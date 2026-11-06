# Aggregate Posterior Profile Estimates for DALY Use

Bridges
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)'s
per-event, per-draw state-probability array into the hospital x pathogen
x profile summary DALY calculation needs: posterior mean + credible
interval, panel composition, event/draw support counts, and explicit
`eligible_for_YLL`/`eligible_for_YLD` flags with reasons.

## Usage

``` r
profile_aggregate_for_daly(
  est,
  ids,
  pathogen,
  cohort_definition,
  sampler_check = NULL,
  ci_level = 0.95,
  min_n_events = 10L,
  min_n_draws = 100L
)
```

## Arguments

- est:

  List from
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md),
  called with `estimand = "observed_preserving_conditional"` (the DALY-
  appropriate estimand – never overwrites a tested cell). A warning is
  issued if `est$estimand` is anything else.

- ids:

  Data frame, the SAME `ids` passed to
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  (`event_id`, `group`).

- pathogen:

  Character, length 1.
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  always operates on one pathogen's panel at a time; this just labels
  the output.

- cohort_definition:

  Character, length 1. Free-text description of which events were fit
  (e.g. `"Discharged only"`, `"all events, no outcome restriction"`) –
  echoed into every row.

- sampler_check:

  List or `NULL`. Return value of
  [`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md).
  `sampler_check$status == "pass"` is required for
  `eligible_for_YLL`/`eligible_for_YLD`. `NULL` (default): both are
  always `FALSE` with reason `"sampler_not_checked"`.

- ci_level:

  Numeric. Default `0.95`.

- min_n_events, min_n_draws:

  Integer. Eligibility floors, same meaning as the legacy function's.
  Default `10L`/`100L`.

## Value

A data frame, one row per group x profile (state): `group`, `pathogen`,
`cohort_definition`, `profile_delta` (state label, e.g. `"RRSS"`),
`R_mean`, `R_lower`, `R_upper`, `n_events`, `n_draws`,
`n_profile_classes`, `sampler_status`, `eligible_for_YLL`,
`exclusion_reason_YLL`, `eligible_for_YLD`, `exclusion_reason_YLD`
(identical to the YLL columns here – there is only one cohort to be
eligible for), `estimator` (`"joint_bernoulli"`).

## Details

This does NOT split into `R_ALL`/`R_KNOWN_OUTCOME`/`R_NF` cohorts after
fitting. In this pipeline, cohort selection (e.g. a Discharged/
non-fatal restriction) happens BEFORE fitting, via
[`profile_select_cohort`](https://saketlab.github.io/anumaan/reference/profile_select_cohort.md)
– the fit is already trained on exactly the cohort a DALY estimand
wants, never on a broader population split apart afterward. Reporting a
post-hoc "R_ALL" alongside an "R_NF" from a model that was only ever fit
on non-fatal events would imply a second population the model was never
trained on and has no actual information about. Instead, this reports
ONE `R_mean`/`R_lower`/`R_upper` per group x profile, for whatever
cohort was actually fit, carried in the `cohort_definition` column so
the number is never mistaken for a different population's. A caller who
genuinely needs both an all-events and a non-fatal-only DALY estimand
gets both by running
[`profile_select_cohort`](https://saketlab.github.io/anumaan/reference/profile_select_cohort.md)
twice and fitting twice – two real fits, not one fit relabeled twice.

A second departure: `sampler_acceptable` has no default here. The legacy
function defaulted it to `TRUE` (opt-out: eligible unless told
otherwise). This function requires an actual
[`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md)
result and treats a missing one as not eligible, with
`exclusion_reason_YLL/YLD = "sampler_not_checked"` – never silently
assumed acceptable.
