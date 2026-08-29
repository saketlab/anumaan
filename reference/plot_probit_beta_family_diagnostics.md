# Fixed-effect convergence broken out by source of variation

Splits the fixed-effect ("beta") parameter group by which covariate
family each coefficient belongs to (Intercept / Hospital / Age / Sex /
Specimen / Ward-ICU / Other, by covariate-name pattern), so a reviewer
can see directly whether e.g. hospital contrasts specifically are
destabilising the fit, rather than only a single aggregate Rhat/ESS for
all of "beta".

## Usage

``` r
plot_probit_beta_family_diagnostics(fit, title_base = "")
```

## Arguments

- fit:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- title_base:

  Character. Prefixed to plot titles.

## Value

Named list with `table` (a formatted summary page) and `bar` (a
percent-Rhat-above-1.01-by-family bar chart), or `NULL` if fixed-effect
parameter diagnostics are unavailable.
