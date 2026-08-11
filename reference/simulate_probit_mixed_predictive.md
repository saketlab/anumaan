# Simulate mixed predictive replicate datasets: new random-effect levels for a requested block, retaining fitted hyperparameters

For each requested block \\r\\ in `blocks_to_replicate`, retains the
posterior draws of \\\tau_r^{(s)}\\ and the block's own fitted
cross-class correlation \\R\_{block\[r\]}^{(s)}\\ (a generated quantity
already stored in `fitted_model$draws`), but draws
`n_new_levels_per_block` ENTIRELY NEW synthetic level effects
\$\$u\_{r,\text{new}}^{(s)} \sim MVN(0, \Sigma_r^{(s)})\$\$ per
posterior state, where \\\Sigma_r^{(s)} = \mathrm{diag}(\tau_r^{(s)})\\
R\_{block\[r\]}^{(s)}\\\mathrm{diag}(\tau_r^{(s)})\\. Every real fitted
event is deterministically (seeded) reassigned to one of the
`n_new_levels_per_block` new synthetic levels for each replicated block,
so a complete replicate outcome matrix can be generated exactly as in
[`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md).
Blocks NOT in `blocks_to_replicate` use their FITTED posterior random
effects (`unreplicated_block_behavior = "retain_fitted"`, default) or
cause an error (`"error"`) if any declared block is left unhandled.

## Usage

``` r
simulate_probit_mixed_predictive(
  fitted_model,
  blocks_to_replicate,
  n_states = 500L,
  seed = 123L,
  n_new_levels_per_block = 1L,
  unreplicated_block_behavior = c("retain_fitted", "error")
)
```

## Arguments

- fitted_model:

  A fitted model object as returned by
  [`fit_bayesian_multivariate_probit`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- blocks_to_replicate:

  Character vector; a subset of
  `fitted_model$random_effects_prep$block_names` to draw NEW levels for.

- n_states:

  Integer; number of posterior states to use. Default 500.

- seed:

  Integer seed (deterministic).

- n_new_levels_per_block:

  Integer; how many new synthetic levels to draw per replicated block.
  Default 1.

- unreplicated_block_behavior:

  `"retain_fitted"` (default) uses fitted posterior random effects for
  every declared block not in `blocks_to_replicate`; `"error"` requires
  every declared block to be listed in `blocks_to_replicate`.

## Value

An object of class `c("amr_mixed_predictive_draws", "amr_ppc_draws")`
with `setup`, `generate_state`, and the other fields shared with
[`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md)'s
return value, plus `blocks_to_replicate`, `unreplicated_blocks`,
`n_new_levels_per_block`, `unreplicated_block_behavior`.

## Details

This is NOT standard posterior predictive checking (which conditions on
the observed groups' own fitted random effects) – do not conflate the
two. `blocks_to_replicate` is currently exercised for a single block at
a time in this package's own tests (the hospital-level use case
explicitly scoped by Part 13 of the predictive-checking specification);
multi-block replication is implemented generically but not yet covered
by synthetic recovery tests.

## References

Stan Development Team. "Posterior and Prior Predictive Checks." Stan
User's Guide.
<https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html>
