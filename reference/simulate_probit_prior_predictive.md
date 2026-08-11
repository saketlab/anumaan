# Simulate prior predictive replicate datasets for a probit model

Draws `n_states` independent parameter states directly from the declared
priors (no NUTS/no conditioning on outcomes) and generates a complete
replicate outcome matrix from each, using the REAL fixed-effect design
matrix and random-effect grouping structure of `fitted_model` as
conditioning predictors (per Stan User's Guide-style prior predictive
simulation). Reuses the identical generative mechanism as
[`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md)
(Bernoulli(Phi(mu)) for identity residual,
`Z ~ MVN(mu, Omega); Y = I(Z>0)` for correlated residual) – only the
source of `theta` differs (prior vs. posterior draws).

## Usage

``` r
simulate_probit_prior_predictive(
  fitted_model,
  n_states = 1000L,
  seed = 123L,
  prior_config_override = NULL,
  preserve_observation_mask = TRUE,
  return_replicates = FALSE
)
```

## Arguments

- fitted_model:

  A fitted model object as returned by
  [`fit_bayesian_multivariate_probit`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
  – supplies the real `X_event` design matrix, random-effect grouping
  structure, class panel, and `prior_config_used` (the priors to
  replicate exactly). Only design/configuration fields are read; no
  posterior draws are used.

- n_states:

  Integer; number of independent prior states to draw. Default 1000.

- seed:

  Integer seed (deterministic).

- prior_config_override:

  Optional named list overriding any of `beta_sd`, `tau_sd`, `lkj_eta`
  from `fitted_model$prior_config_used` – e.g. to prior-predictive-check
  a CANDIDATE prior before fitting.

- preserve_observation_mask, return_replicates:

  See
  [`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md).

## Value

An `"amr_prior_predictive_draws"`/`"amr_ppc_draws"` object with the same
shape as
[`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md)'s
return value.

## Details

Returns an object of class
`c("amr_prior_predictive_draws", "amr_ppc_draws")` – structurally
compatible with
[`compute_probit_ppc_statistics`](https://saketlab.github.io/anumaan/reference/compute_probit_ppc_statistics.md)
(which accepts any `"amr_ppc_draws"`-classed object), so the same
AMR-specific discrepancy-statistic machinery can summarise the prior
predictive distribution's own properties (see
[`compute_prior_predictive_status`](https://saketlab.github.io/anumaan/reference/compute_prior_predictive_status.md),
which additionally derives the prior-specific plausibility summaries in
Part 11 of the predictive- checking specification: fraction of
event-class probabilities near 0/1, fraction of degenerate
all-resistant/all-susceptible profiles, and hospital-level spread).

## References

Stan Development Team. "Posterior and Prior Predictive Checks." Stan
User's Guide.
<https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html>
