# Estimate Resistance Profiles From a Fitted Joint Bernoulli Model

Produces per-event joint-profile estimates under one of three EXPLICITLY
DIFFERENT estimands – `estimand` has no default, by design:

## Usage

``` r
profile_estimate(
  fit_result,
  design,
  Y,
  ids,
  estimand,
  n_draws = 500L,
  ci_level = 0.95,
  seed = 123L
)
```

## Arguments

- fit_result:

  List from
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  (reads `$fit`), or a raw `cmdstanr` `CmdStanMCMC`.

- design:

  List from
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (reads `$stan_data` and `$design$states`).

- Y, ids:

  The SAME `Y` matrix and `ids` data frame passed to
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md).

- estimand:

  Character, one of `"unconditional"`,
  `"observed_preserving_conditional"`, `"missing_allocation"`. No
  default – the caller must choose.

- n_draws:

  Integer. Posterior draws used (subsampled if the fit has more).
  Default `500L`. Note: for `"unconditional"`/
  `"observed_preserving_conditional"`, the returned
  `profile_probability_draws` array is `n_draws x nrow(Y) x 2^D` –
  memory scales with all three, so reduce `n_draws` for a very large `Y`
  with a large panel.

- ci_level:

  Numeric in (0, 1). Credible-interval level for
  `profile_probability_lower`/`upper` and
  `class_probability_lower`/`upper`. Default `0.95` (giving the
  `[Q_.025, Q_.975]` marginal CrI). Ignored for
  `estimand = "missing_allocation"` (its uncertainty is already the full
  set of discrete per-draw allocations, not a further interval). These
  are MARGINAL per-event/per-class/per-state CrIs, each computed by
  aggregating WITHIN a draw first (summing state probabilities into a
  class margin, say) and taking quantiles ACROSS draws second – never
  the reverse, and never a simultaneous region across multiple
  events/classes/states at once.

- seed:

  Integer. Default `123L`.

## Value

A list with `estimand` (echoes the argument, for provenance) and
`state_labels` (character vector, length `2^D`, e.g. `"RRSS"`), plus:

- `"unconditional"` / `"observed_preserving_conditional"`:

  `profile_probability`: matrix, events x `2^D` states, the
  posterior-mean probability of each state per event;
  `profile_probability_lower`/`profile_probability_upper`: same shape,
  the `ci_level` marginal CrI per event x state;
  `profile_probability_draws`: the full `n_draws x events x 2^D` array;
  `class_probability`/
  `class_probability_lower`/`class_probability_upper`: matrix, events x
  `D` classes – the per-class resistance margin (states summed WITHIN
  each draw, then meaned/quantiled across draws), with its own CrI, not
  derived by summing the per-state CrIs (that would NOT be a valid
  interval for the sum)); `class_probability_draws`: the full
  `n_draws x events x D` array behind `class_probability` (e.g. for
  feeding
  [`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md)'s
  `derived_estimates` MCSE check).

- `"missing_allocation"`:

  `allocated_state`: integer matrix, `n_draws x events`, each entry a
  state index into `state_labels`.

## Details

- `"unconditional"`:

  The model's own implied profile distribution from an event's
  covariates alone, ignoring whatever was actually tested for that
  event. A population-level baseline, never a per-patient "best guess"
  when AST results exist for that patient.

- `"observed_preserving_conditional"`:

  Never overwrites a tested cell: distributes probability mass only over
  the states compatible with the event's own observed R/S pattern (via
  [`profile_compatible_states`](https://saketlab.github.io/anumaan/reference/profile_compatible_states.md))
  – exactly the same exact marginalization the Stan model's own
  likelihood already uses for untested classes, not an approximation. A
  fully-tested event collapses to 100 reduces exactly to the
  unconditional distribution (every state is "compatible" when nothing
  was tested). This is the estimand DALY-style burden calculations need:
  never claim to know more about a patient's resistance pattern than was
  actually measured.

- `"missing_allocation"`:

  Like `"observed_preserving_conditional"`, but instead of a continuous
  probability over compatible states, DRAWS one specific state per event
  per posterior draw from that same conditional distribution – carrying
  allocation uncertainty forward as a set of discrete posterior draws,
  for downstream calculations that need an actual profile assignment per
  patient (e.g. counting how many patients fall in profile RRRR) with
  Monte Carlo uncertainty across draws, rather than fractional
  probability mass.
