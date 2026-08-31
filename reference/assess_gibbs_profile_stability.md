# Assess Numerical Stability of Correlated Profile Completion

Repeats the existing conditional truncated-MVN Gibbs profile completion
for fixed posterior draws and a fixed, reproducibly selected subset of
incomplete events. It does not refit Stan or alter the fitted model. The
comparison therefore isolates finite Gibbs-chain length. The conditional
completion algorithm is the correct profile-completion algorithm; this
diagnostic checks whether its finite Monte Carlo run length is adequate
for a particular fit, especially when posterior residual correlations
are high.

## Usage

``` r
assess_gibbs_profile_stability(
  fitted_model,
  schedules = list(baseline = list(burnin = 30L, kept = 50L), medium = list(burnin = 60L,
    kept = 100L), long = list(burnin = 100L, kept = 250L)),
  n_posterior_draws = 200L,
  max_events = 250L,
  seed = 123L,
  tolerance = 0.01,
  event_selection = c("stratified", "all")
)
```

## Arguments

- fitted_model:

  Object returned by \[fit_bayesian_multivariate_probit()\].

- schedules:

  Named list of \`list(burnin=, kept=)\` schedules. It must contain
  \`baseline\`, \`medium\`, and \`long\`.

- n_posterior_draws:

  Number of fixed posterior draws to use.

- max_events:

  Maximum number of incomplete events to assess.

- seed:

  Seed used once to select draws/events and then reset before each
  schedule. Thus schedules share posterior draws, events, and RNG
  streams.

- tolerance:

  Maximum aggregate difference for status \`"pass"\`; twice this value
  is the \`"warning"\` upper bound.

- event_selection:

  Currently \`"stratified"\` or \`"all"\`.

## Value

A structured list with summary, event/aggregate comparisons, selected
events, schedules, posterior draw indices, and a separate stability
status.

## Details

The default 0.01/0.02 thresholds are numerical-stability conventions for
this diagnostic, not universal statistical thresholds.
