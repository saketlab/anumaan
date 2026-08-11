# Simulate posterior predictive replicate datasets from a fitted probit model

For each of `n_states` saved posterior states \\s\\ and each fitted
event \\e\\, reconstructs \$\$\mu_e^{(s)} = X_e \beta^{(s)} + \sum_r
u\_{r,g_r(e)}^{(s)}\$\$ using the existing generic random-effect
machinery
([`prepare_random_effects`](https://saketlab.github.io/anumaan/reference/prepare_random_effects.md)/[`re_contribution`](https://saketlab.github.io/anumaan/reference/re_contribution.md)
– never hardcoded block names), then draws a COMPLETE replicate outcome
matrix \\Y\_{rep}^{(s)}\\ unconditionally from \\P(Y \mid \mu^{(s)},
\Omega^{(s)})\\: `Bernoulli(Phi(mu))` per class for identity residual
models, or `Z ~ MVN(mu, Omega); Y = I(Z > 0)` for correlated residual
models (see Stan User's Guide,
<https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html>).

## Usage

``` r
simulate_probit_posterior_predictive(
  fitted_model,
  n_states = 500L,
  seed = 123L,
  preserve_observation_mask = TRUE,
  return_replicates = FALSE,
  replicate_random_effects = FALSE,
  random_effect_blocks_to_replicate = NULL
)
```

## Arguments

- fitted_model:

  A fitted model object as returned by
  [`fit_bayesian_multivariate_probit`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
  (or the equivalent `fit_light` object saved by the analysis layer with
  `$fit` stripped – this function never touches `$fit`, only stored
  draws and design matrices).

- n_states:

  Integer; number of saved posterior states to use (a random subsample,
  seeded, if fewer than the total number of saved draws). Default 500;
  use 1000-2000 for shortlisted models.

- seed:

  Integer seed for state subsampling and stochastic generation.

- preserve_observation_mask:

  Logical, default `TRUE`. If `TRUE`, the returned generator's `Y_rep`
  slot masks the complete generated replicate to exactly the same
  tested/untested cell pattern as the real observed AST matrix, so that
  discrepancy statistics computed downstream compare like with like. The
  COMPLETE (unmasked) replicate is always separately available as
  `Y_rep_complete` for statistics that genuinely need it (e.g.
  complete-profile statistics restricted to adequately-supported
  fully-observed subsets) – the two are never compared without explicit
  qualification.

- return_replicates:

  Logical, default `FALSE`. If `TRUE`, materialises full
  `S x N_events x D` replicate arrays (`Y_rep_array`, masked per
  `preserve_observation_mask`, and `Y_rep_complete_array`, always
  unmasked) – small debug datasets only; refuses above a fixed
  element-count ceiling.

- replicate_random_effects, random_effect_blocks_to_replicate:

  Reserved for API compatibility with the (distinct) mixed predictive
  check. Not implemented here – passing a non-default value errors with
  a pointer to
  [`simulate_probit_mixed_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_mixed_predictive.md)
  rather than silently being ignored.

## Value

An object of class `"amr_ppc_draws"`: a list with `setup` (internal
draw-setup object), `generate_state` (a function of one integer argument
`s` returning `list(Y_rep_complete, Y_rep)` for that state),
`n_states_used`, `seed_used`, `preserve_observation_mask`,
`residual_structure`, `return_replicates`, and (only when
`return_replicates = TRUE`) `Y_rep_array`/`Y_rep_complete_array`.

## Details

This function NEVER calls the internal DALY conditional-completion
helper
([`.gibbs_conditional_profile_probs()`](https://saketlab.github.io/anumaan/reference/dot-gibbs_conditional_profile_probs.md))
– that function samples missing cells conditional on an event's own
observed AST signs, which is a fundamentally different (and, for this
purpose, wrong) generative mechanism. See the file header of
`probit_posterior_predictive.R` and `test-probit-posterior-predictive.R`
for the explicit distinction and regression test.

By default (`return_replicates = FALSE`), no `S x N_events x D` array is
ever materialised – this function returns a lightweight generator object
(class `"amr_ppc_draws"`) that
[`compute_probit_ppc_statistics`](https://saketlab.github.io/anumaan/reference/compute_probit_ppc_statistics.md)
drives one state at a time, streaming, to bound peak memory at
O(N_events x D) rather than O(S x N_events x D). Each state's replicate
is generated under a temporary, locally-seeded RNG state derived
deterministically from `(seed, s)` alone (restoring the caller's global
RNG state afterwards) – so a given state's draw depends only on `seed`,
`s`, and that state's posterior parameters, never on call order, on
interleaving with any other generator object, or on any other RNG use
elsewhere in the calling session. Draws are cached per state, so calling
the generator repeatedly for the same `s` is idempotent.

`return_replicates = TRUE` additionally materialises the full replicate
arrays and is intended for small debugging datasets only; a safety
ceiling refuses to silently run this on a large realistic dataset
(reduce `n_states` or leave `return_replicates = FALSE` instead).

## References

Stan Development Team. "Posterior and Prior Predictive Checks." Stan
User's Guide.
<https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html>
