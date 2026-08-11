# Classify plausibility of the prior predictive distribution

Unlike posterior predictive checking, the goal here is NOT to make the
simulated distribution match the observed data closely – it is to flag
obviously implausible/extreme implications of the CHOSEN priors, before
conditioning on any outcomes. Computes, per prior state and averaged
across states: the fraction of event-class probabilities
\\\Phi(\mu\_{ed})\\ below 0.001 or above 0.999 (near-deterministic
implied resistance), the fraction of generated events that are
degenerate (all classes resistant or all susceptible), and the spread of
per-"hospital" (the model's primary declared random-effect grouping)
mean implied probability – a very large spread implies an implausibly
extreme facility-to-facility prior.

## Usage

``` r
compute_prior_predictive_status(prior_draws, thresholds = list())
```

## Arguments

- prior_draws:

  An `"amr_prior_predictive_draws"` object from
  [`simulate_probit_prior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_prior_predictive.md).

- thresholds:

  Named list overriding any of the defaults:
  `max_fraction_extreme_probability = 0.10`,
  `severe_fraction_extreme_probability = 0.50`,
  `max_fraction_degenerate_profiles = 0.50`,
  `severe_fraction_degenerate_profiles = 0.90`,
  `max_hospital_spread_sd = 0.35`. All thresholds are project decisions,
  not universal truths – override freely and document why.

## Value

List with `status` (one of `"pass"`,
`"warning_extreme_prior_predictions"`,
`"warning_excessive_hospital_heterogeneity"`,
`"warning_degenerate_profiles"`, `"fail_implausible_prior_predictive"`,
`"insufficient_prior_check"`), `reasons`, `thresholds_used`, and
`summary` (the averaged plausibility fractions).
