# Fixed-effect coefficient diagnostics, split by contrast type and faceted by antimicrobial class

Returns a named list of up to two ggplot objects: `other` (every
non-hospital-dummy covariate – intercept, age, gender, etc.) and
`hospital` (hospital contrasts only, if the fit used hospital as a fixed
effect rather than a random effect). Both facet by antimicrobial class
rather than putting every class's coefficients on one axis.

## Usage

``` r
plot_probit_fixed_effect_diagnostics(fit, title_base = "", ci_level = 0.95)
```

## Arguments

- fit:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- title_base:

  Character. Prefixed to plot titles.

- ci_level:

  Numeric. Credible interval level. Default `0.95`.

## Value

Named list of ggplot objects (`other`, `hospital`), either of which may
be absent if that group has no covariates in this fit. `NULL` if
`fit$draws` or `fit$X_event` is missing.
