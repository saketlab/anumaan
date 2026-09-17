# Reduce one item's (observed_value, T_rep across states) to the Part 7 output schema: statistic_name, stratum, observed_value, replicated_mean/ sd/q025/q50/q975, ppc_tail_probability, ppc_two_sided, n_replications, support_status. ppc_tail_probability/ppc_two_sided are posterior predictive TAIL-PROBABILITY-LIKE quantities. The reported one-sided value is the inclusive upper tail, `(1 + sum(T_rep >= T_obs)) / (B + 1)`; the two-sided value is `min(1, 2 * min(lower_tail, upper_tail))`, with the lower tail calculated analogously. These finite-replicate corrected summaries are NOT classical calibrated p-values.

Reduce one item's (observed_value, T_rep across states) to the Part 7
output schema: statistic_name, stratum, observed_value, replicated_mean/
sd/q025/q50/q975, ppc_tail_probability, ppc_two_sided, n_replications,
support_status. ppc_tail_probability/ppc_two_sided are posterior
predictive TAIL-PROBABILITY-LIKE quantities. The reported one-sided
value is the inclusive upper tail,
`(1 + sum(T_rep >= T_obs)) / (B + 1)`; the two-sided value is
`min(1, 2 * min(lower_tail, upper_tail))`, with the lower tail
calculated analogously. These finite-replicate corrected summaries are
NOT classical calibrated p-values.

## Usage

``` r
.ppc_summarize_statistic(
  statistic_name,
  stratum,
  observed_value,
  T_rep,
  support_status,
  ci_level
)
```
