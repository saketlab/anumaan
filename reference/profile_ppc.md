# Posterior Predictive Checks for a Joint Bernoulli Fit

Observed-vs-model comparison at three levels – marginal (per class),
pairwise (per class pair, all four RR/RS/SR/SS cells), and
complete-profile (per full `2^D` state) – plus a `"not_computable"`
status for any group x check cell with too little direct support to
evaluate at all. Nothing below the relevant `min_*` floor is silently
omitted; it is reported with `status = "not_computable"` and `NA`
numeric columns.

## Usage

``` r
profile_ppc(
  fit_result,
  design,
  Y,
  ids,
  n_draws = 500L,
  min_tested = 10L,
  min_cotested = 10L,
  min_complete = 10L,
  min_count_group = 10L,
  min_partial_pattern = 10L,
  covariate_subgroups = NULL,
  min_subgroup = 10L,
  ci_level = 0.95,
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

- n_draws:

  Integer. Posterior draws used (subsampled if the fit has more).
  Default `500L`.

- min_tested, min_cotested, min_complete:

  Integer. Below these, a group x check cell is `"not_computable"`
  rather than evaluated. Default `10L` each.

- min_count_group:

  Integer. Below this many COMPLETE events in a group, the
  `count_groups` family is `"not_computable"` for that group. Default
  `10L`.

- min_partial_pattern:

  Integer. Below this many events sharing an EXACT partial test mask (in
  a group), that (group, mask) combination is skipped from
  `partial_patterns` entirely (never fabricated as a zero-support row –
  there are too many possible masks to report every absent one; only
  masks actually meeting this floor are listed, with every sub-profile
  of a LISTED mask reported, including its zero-count ones). Default
  `10L`.

- covariate_subgroups:

  Named list of vectors/factors, or `NULL` (default). Each vector must
  have length `nrow(Y)`, in the SAME row order as `Y`/`ids` – e.g.
  `list(specimen_type = ..., age_band = ...)`. When supplied, adds the
  `covariate_subgroups` table: the SAME marginal-level check as
  `marginal`, stratified by each named vector's levels instead of
  `group` (pooled across groups – these are prespecified clinical
  subgroup checks, not another hospital breakdown). Default `NULL` (the
  returned `covariate_subgroups` table has zero rows).

- min_subgroup:

  Integer. Below this many tested events in a subgroup level x class
  cell, that cell is `"not_computable"`. Default `10L`.

- ci_level:

  Numeric. Default `0.95`.

- seed:

  Integer. Draw-subsampling AND replicated-outcome simulation seed.
  Default `123L`.

## Value

A list of data frames:

- `marginal`:

  one row per group x class: `status`, `n_tested`, `n_resistant`,
  `n_susceptible`, `observed_resistance`, `expected_mean/lower/upper` (=
  `model_mean/lower/upper`), `predictive_mean/lower/upper`,
  `absolute_error`, `interval_contains_observed`,
  `predictive_interval_contains_observed`

- `pairwise`:

  one row per group x class-pair: `status`, `n_cotested`,
  `n_RR/n_RS/n_SR/n_SS`, `observed_RR_prop` etc.,
  `expected_RR_mean/lower/upper` (= `model_RR_mean/lower/upper`, and
  `_RS`/`_SR`/ `_SS`), `predictive_RR_mean/lower/upper`,
  `absolute_error` (RR cell), `interval_contains_observed`,
  `predictive_interval_contains_observed`

- `complete_profile`:

  one row per group x state (one `"not_computable"` row per group below
  `min_complete`, state `NA`): `status`, `n_complete`, `state`,
  `empirical_count/frequency`, `expected_mean/lower/upper` (=
  `model_mean/lower/upper`), `predictive_mean/lower/upper`,
  `absolute_error`, `interval_contains_observed`,
  `predictive_interval_contains_observed`

- `count_groups`:

  one row per group x resistance-count `0..D` (COMPLETE events only –
  total number of resistant classes): `status`, `n_complete`, `count`,
  `empirical_count/frequency`, `expected_mean/lower/upper`,
  `predictive_mean/lower/upper`, `absolute_error`,
  `interval_contains_observed`, `predictive_interval_contains_observed`

- `partial_patterns`:

  one row per (group, exact test mask meeting `min_partial_pattern`,
  sub-profile over the TESTED classes in that mask) – the denominator is
  EVERY event sharing that EXACT mask, not events co-tested as part of a
  larger panel (that is `pairwise`'s job) and not limited to one pair:
  an event tested on 3 of 5 panel classes contributes to exactly one
  mask's sub-profile breakdown here, at whatever order that mask is.
  Columns: `status`, `mask` (tested class names, "; "-joined), `n_mask`,
  `sub_profile` (R/S string over just the tested classes, in panel
  order), `empirical_count/ frequency`, `expected_mean/lower/upper`,
  `predictive_mean/lower/upper`, `absolute_error`,
  `interval_contains_observed`, `predictive_interval_contains_observed`.
  EVERY sub-profile of a listed mask is reported, including zero-count
  ones.

- `covariate_subgroups`:

  empty (zero rows) unless `covariate_subgroups` is supplied; one row
  per (subgroup name, level, class), same columns as `marginal` but
  `group` replaced by `subgroup_name`/`subgroup_level`.

## Details

TWO different quantities are reported for every evaluated cell, and they
answer different questions – do not treat them interchangeably:

- `expected_mean`/`expected_lower`/`expected_upper` (also exposed as
  `model_mean`/`model_lower`/ `model_upper` for backward compatibility):

  the posterior mean and credible interval of the model's EXPECTED
  probability – an exact closed-form sum over state probabilities (the
  same softmax-over-states construction
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)/
  the Stan model use for the likelihood). This interval reflects
  PARAMETER uncertainty only: a model that is certain `p = 0.5` reports
  a near-zero-width interval here, even though the actual count among
  replicated events would still vary.

- `predictive_mean`/`predictive_lower`/ `predictive_upper`:

  a genuine posterior PREDICTIVE interval, simulated per posterior draw:
  for every event, one complete profile is sampled from that draw's
  joint state distribution, the event's ORIGINAL observation mask is
  reapplied (only classes actually tested in the real data contribute a
  replicated value), and the SAME replicated profile is used to compute
  the replicated marginal, pairwise, and complete-profile statistics for
  that event within that draw – so the joint dependence the model was
  fit to capture is preserved across all three levels, not simulated
  independently per level. This interval reflects BOTH parameter
  uncertainty and genuine outcome-level sampling variation, per Stan's
  own posterior-predictive-vs-posterior-expectation distinction
  (<https://mc-stan.org/docs/stan-users-guide/posterior-prediction.html>).

`interval_contains_observed` (unchanged from earlier versions of this
function) checks against the EXPECTED interval, for backward
compatibility with
[`profile_release`](https://saketlab.github.io/anumaan/reference/profile_release.md)'s
rollup; `predictive_interval_contains_observed` checks against the
PREDICTIVE interval – the two can disagree, and a cell that fails the
tighter expected-probability check may still pass the (necessarily
wider) predictive check.

This also avoids the confirmed defect in the J2-L four-cell-sensitivity
experiment's `fc_partial_specs()` (grouping comparison rows by the
observed outcome pattern itself, which makes "observed proportion" 1.0
by construction, not a real check): every "observed" count here –
`n_tested`/`n_resistant`, `n_cotested`/`n_RR` etc., empirical state
frequency – is a plain tabulation from `Y` against a FIXED category (a
class, a pair, a state), never derived by matching each row's own label
back to itself.
