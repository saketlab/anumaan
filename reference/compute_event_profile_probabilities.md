# Compute Posterior Resistance Profile Probabilities via MVN Simulation

For each posterior draw, constructs the event-level linear predictor
\\\mu_e\\ using fixed effects and correlated random effects from the
fitted model, simulates latent \\Z_e = \mu_e +
L\_\Omega\\\varepsilon_e\\, \\\varepsilon_e \sim N(0,I_D)\\, converts to
binary resistance profiles, and accumulates profile probabilities. All
\\2^D\\ profiles appear in the output; profiles not observed in
simulation receive probability 0.

## Usage

``` r
compute_event_profile_probabilities(
  fitted_model,
  n_posterior_draws_for_profiles = 2000L,
  outcome_col = NULL,
  nonfatal_values = c("Discharged", "Survived", "Alive", "discharged", "survived",
    "alive"),
  seed = 123L
)
```

## Arguments

- fitted_model:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- n_posterior_draws_for_profiles:

  Integer. Number of posterior draws to use for profile simulation.
  Subsampled without replacement when total draws exceed this value.
  Each draw generates one simulated latent profile per event; rare
  profiles may be underestimated when this is small. For finer Monte
  Carlo resolution of rare profiles, consider increasing this value or
  adding `n_predictive_replicates_per_draw` in a future extension.
  Default `2000L`.

- outcome_col:

  Character or `NULL`. Patient outcome column in
  `fitted_model$event_metadata`. When `NULL`, all events are treated as
  having a known outcome and R_NF is not computed separately.

- nonfatal_values:

  Character vector. Outcome values for the non-fatal cohort (R_NF).
  Default covers common discharge/survival labels.

- seed:

  Integer. Random seed for MVN simulation. Default `123L`.

## Value

Named list: `event_profiles` (event-level posterior means) and
`aggregate_draws` (per-draw R_ALL and R_NF per hospital x pathogen x
profile, used by
[`aggregate_profiles_for_daly()`](https://saketlab.github.io/anumaan/reference/aggregate_profiles_for_daly.md)
for credible intervals). Both tibbles contain all \\2^D\\ profiles per
hospital-pathogen pair.

## Details

**Estimand:** The posterior predictive distribution over the observed
event case-mix – the profile distribution you would see if you drew a
new event uniformly from the set of observed events (same covariate
distribution as the data). This is labelled
`"observed_stewardship_event_mix"`.
