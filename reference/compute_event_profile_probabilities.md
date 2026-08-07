# Compute Observed-Plus-Imputed Resistance Profile Probabilities

For each event, retains every observed (tested) AST cell exactly and
computes a probability distribution only over the genuinely untested
(`NA`) antibiotic-class cells in that event's profile panel. Any
enumerated profile inconsistent with the event's observed cells receives
probability 0 – this function never overwrites an observed R/S result.

## Usage

``` r
compute_event_profile_probabilities(
  fitted_model,
  n_posterior_draws_for_profiles = 2000L,
  outcome_col = NULL,
  nonfatal_values = c("Discharged", "Survived", "Alive", "discharged", "survived",
    "alive"),
  seed = 123L,
  n_gibbs_burnin = 10L,
  n_gibbs_kept = 20L
)
```

## Arguments

- fitted_model:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- n_posterior_draws_for_profiles:

  Integer. Number of posterior draws averaged over for both the
  event-level profile probabilities and the draw-level aggregates.
  Subsampled without replacement when total draws exceed this value.
  Default `2000L`.

- outcome_col:

  Character or `NULL`. Patient outcome column in
  `fitted_model$event_metadata`. When `NULL`, all events are treated as
  having a known outcome and R_NF is not computed separately.

- nonfatal_values:

  Character vector. Outcome values for the non-fatal cohort (R_NF).
  Default covers common discharge/survival labels.

- seed:

  Integer. Random seed. Default `123L`.

- n_gibbs_burnin, n_gibbs_kept:

  Integer. Only used when
  `fitted_model$residual_structure == "correlated"` – see
  [`.gibbs_conditional_profile_probs()`](https://saketlab.github.io/anumaan/reference/dot-gibbs_conditional_profile_probs.md).
  Defaults `10L`/`20L`.

## Value

Named list: `event_profiles` (event-level posterior mean
observed-plus-imputed profile probabilities, with
`n_classes_observed`/`n_classes_missing` per event) and
`aggregate_draws` (per-draw `R_ALL` \[truly all events in the hp pair\],
`R_KNOWN_OUTCOME` \[known-outcome cohort\], and `R_NF` \[nonfatal
cohort\], used by
[`aggregate_profiles_for_daly()`](https://saketlab.github.io/anumaan/reference/aggregate_profiles_for_daly.md)
for credible intervals). Both tibbles contain all \\2^D\\ profiles per
hospital-pathogen pair and carry `profile_generation_method` and
`panel_reason` columns.

## Details

**Identity residual structure**
(`fitted_model$residual_structure == "identity"`): classes are
conditionally independent given the linear predictor, so the
missing-cell probabilities are computed analytically as \\P(Y\_{ed}=1
\mid \theta) = \Phi(\mu\_{ed})\\ for each missing class \\d\\, and the
profile probability over the missing dimensions is the exact product of
per-class Bernoulli probabilities – no latent-variable simulation is
used, so there is no added latent-profile-simulation noise. Rows carry
`profile_generation_method = "conditional_analytic_identity"`.

**Correlated residual structure**: the missing latent dimensions are
sampled conditional on the observed resistance SIGNS of the tested
classes via a Gibbs sampler on the truncated multivariate normal – see
[`.gibbs_conditional_profile_probs()`](https://saketlab.github.io/anumaan/reference/dot-gibbs_conditional_profile_probs.md).
This replaces the earlier unconditional latent-\\Z\\ simulation, which
resampled tested cells along with untested ones and was therefore not
observed-plus-imputed. Rows carry
`profile_generation_method = "conditional_gibbs_correlated"` and ARE
eligible for downstream DALY use (subject to the same panel-support and
sampler-acceptability gates as identity-residual profiles – see
[`aggregate_profiles_for_daly()`](https://saketlab.github.io/anumaan/reference/aggregate_profiles_for_daly.md)).
Because Gibbs is only run for `n_gibbs_burnin + n_gibbs_kept` iterations
per draw rather than to full convergence, correlated-residual profiles
carry additional Monte Carlo error beyond identity-residual profiles at
the same `n_posterior_draws_for_profiles`; increase `n_gibbs_kept`
(and/or `n_gibbs_burnin`) for lower-noise estimates at higher compute
cost.

Both residual structures produce, per event per posterior draw, a
profile probability vector over the panel's \\2^{D_p}\\ enumerated
profiles that (a) sums to exactly 1 and (b) assigns exactly 0 to any
profile inconsistent with that event's observed cells – enforced by
construction in both cases (an exact product for identity; Gibbs never
visiting an inconsistent pattern for correlated), and checked with a
hard [`stop()`](https://rdrr.io/r/base/stop.html) for the identity path
since it has no other source of Monte Carlo noise to blur a genuine bug.

**Class panels**: the antibiotic-class panel enumerated for each
hospital x pathogen pair is drawn from the approved eligibility rules
computed at fit time (`fitted_model$eligibility_report`: marginal
n_tested/n_resistant/n_susceptible thresholds, plus pairwise co-testing
sufficiency for correlated-residual fits) – not from "tested at least
once". See
[`.resolve_profile_class_panel()`](https://saketlab.github.io/anumaan/reference/dot-resolve_profile_class_panel.md).

**Estimand:** The posterior distribution over the observed event
case-mix, conditional on each event's own observed AST results. This is
labelled `"observed_stewardship_event_mix"`.
