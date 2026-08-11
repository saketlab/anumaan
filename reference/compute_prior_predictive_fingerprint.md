# Fingerprint the prior-generative configuration of a fitted model

Two experiments that share the same prior-generative structure (same
residual structure, same beta/tau/LKJ priors, same declared
random-effect block structure, same fixed-effect design-column
definitions, same `D`) can safely REUSE a prior predictive check's
results rather than rerunning
[`simulate_probit_prior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_prior_predictive.md)
for every experiment – prior predictive checking depends only on the
generative structure, never on the fitted posterior. If ANY component
changes (a prior, a scaling, the declared block structure), the
fingerprint changes and any cached result must be invalidated.

## Usage

``` r
compute_prior_predictive_fingerprint(fitted_model)
```

## Arguments

- fitted_model:

  A fitted model object as returned by
  [`fit_bayesian_multivariate_probit`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

## Value

List with `fingerprint` (a single character hash via
[`rlang::hash()`](https://rlang.r-lib.org/reference/hash.html), already
an anumaan dependency – no new package dependency introduced) and
`components` (the exact named list that was hashed, for
audit/debugging).
