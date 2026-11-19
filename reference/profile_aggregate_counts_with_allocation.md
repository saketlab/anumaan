# Aggregate Hospital-Profile Counts Under Discrete Missing-Data Allocation

A complete event's compatible-state set is always a single state (see
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)'s
`"missing_allocation"` doc), so it deterministically allocates to its
own true state on every draw – there is no need to separately add an
"n_observed" term. Passing complete AND partial events together (the
same convention
[`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)
already uses) reproduces the observed count automatically, with the
partial events' allocation uncertainty layered on top.

## Usage

``` r
profile_aggregate_counts_with_allocation(
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
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  with `estimand = "missing_allocation"` (reads `$allocated_state`,
  `$state_labels`).

- ids:

  Data frame with a `group` column, same row order/length as the `Y`
  passed to
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  (one row per event).

- pathogen, cohort_definition:

  Character, echoed into the output for provenance (same convention as
  [`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)).

- sampler_check:

  List or `NULL`. Same convention as
  [`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md):
  reads `$status`; `"pass"` is the only acceptable value, anything else
  (including `NULL`) makes every row ineligible with a named reason.

- ci_level:

  Numeric in (0, 1). Quantile interval level for the realised-count
  distribution. Default `0.95`.

- min_n_events, min_n_draws:

  Integer floors below which a group x profile row is marked ineligible
  (same meaning as
  [`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)).
  Defaults `10L`, `100L`.

## Value

Data frame, one row per group x profile: `N_mean`, `N_lower`, `N_upper`
(the realised completed-count distribution's mean and `ci_level`
quantile interval – integer-valued draws summarised, not a continuous
proportion), `n_events`, `n_draws`, `n_profile_classes`,
`sampler_status`, `computational_checks_passed_YLL`/
`computational_checks_passed_YLD` and their `exclusion_reason_*` (same
computational-check logic and column names as
[`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md),
for a consistent downstream filter regardless of which aggregator
produced a row – and see that function's doc for what these flags do NOT
establish: PPC adequacy, held-out validation, missingness sensitivity,
death/ episode attribution, double counting, or RR/duration/disability
input availability. This is a computational pre-check, not a burden-use
release decision.), and
`estimator = "joint_bernoulli_discrete_allocation"` so the two
aggregators' outputs are never silently confused if concatenated.

## Details

For hospital/group \\h\\ and profile \\s\\, at each posterior draw
\\m\\: \$\$N\_{hs}^{(m)} = \sum\_{e \in h}
I(\mathrm{allocated\\state}\[m, e\] = s)\$\$ Collecting \\N\_{hs}^{(1)},
\ldots, N\_{hs}^{(M)}\\ gives the realised completed-count distribution
for that hospital/profile; this function reports its mean and `ci_level`
quantile interval. This single procedure automatically mixes parameter
uncertainty (the draw-to-draw allocation weights shift because
\\\theta^{(m)}\\ shifts) and allocation uncertainty (even at fixed
weights, which specific state gets drawn is itself random) – the
two-term variance decomposition discussed alongside this function's
design, obtained here by direct simulation rather than a closed-form
formula.
