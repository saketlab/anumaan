# Classify overall posterior predictive fit from a discrepancy-statistics table

Does NOT classify from a single tail probability. Groups supported
statistics into families (`marginal`, `pairwise`, `profile`
\[complete-profile and resistant-count statistics – both address joint
multidrug-resistance burden\], `hospital_heterogeneity`, `cluster`
\[admission and patient clustering\]) and flags a family when either (a)
any statistic in it is SEVERE (`ppc_tail_probability < tail_severe` or
`> 1 - tail_severe`), or (b) the FRACTION of that family's supported
statistics that are extreme exceeds `max_fraction_core_extreme`.
Thresholds are configurable, not universal hardcoded truths, and are
applied as discrepancy FLAGS, not frequentist rejection tests.

## Usage

``` r
compute_posterior_predictive_status(ppc_statistics, thresholds = list())
```

## Arguments

- ppc_statistics:

  Tibble returned by
  [`compute_probit_ppc_statistics`](https://saketlab.github.io/anumaan/reference/compute_probit_ppc_statistics.md).

- thresholds:

  Named list overriding any of the defaults: `tail_warning = 0.025`,
  `tail_severe = 0.005`, `max_fraction_core_extreme = 0.20`. These are
  INITIAL AMR-PROJECT DEFAULTS, not universal statistical truths –
  validated only against this package's own synthetic recovery scenarios
  (`test-probit-predictive-synthetic-recovery.R`). Expect to recalibrate
  them once real multi-pathogen experiment results are available;
  override via this argument rather than editing the defaults in place,
  so the thresholds actually used for any given run remain visible in
  `$thresholds_used`.

## Value

List with `status` (one of `"pass"`, `"warning_marginal_ppc"`,
`"warning_pairwise_ppc"`, `"warning_profile_ppc"`,
`"warning_hospital_heterogeneity_ppc"`, `"warning_cluster_ppc"`,
`"fail_major_ppc_misfit"`, `"insufficient_ppc_support"`), `reasons`
(character vector of every triggered family-level flag),
`thresholds_used`, and `family_status` (per-family
n/n_extreme/n_severe/fraction_extreme).

## Details

Statistics with essentially ZERO posterior-predictive variance
(`replicated_sd` at or near 0 – e.g. a discrete statistic that lands on
the same value in every replicated state, such as `max(resistant_count)`
pinned at its ceiling for a small class panel) are excluded from
extreme/severe classification (though still returned in `ppc_statistics`
with their true, unmodified `ppc_tail_probability`):
`mean(T_rep >= T_obs)` is exactly 1 whenever every replicate exactly
EQUALS the observed value, which is a well-known degeneracy of
posterior-predictive tail probabilities for near-constant discrete
statistics, not evidence of misfit – a discrepancy statistic that cannot
vary under the posterior predictive distribution carries no
discriminating information about model fit.
