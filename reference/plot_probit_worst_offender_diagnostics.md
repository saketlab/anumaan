# Trace and rank plots for the worst-converging structural parameters

Auto-selects the union of the top `n_worst` parameters by Rhat and the
top `n_worst` by lowest bulk ESS (from
`fit$diagnostics_detail$monitored_parameters`), and shows their chain
trace and rank-overlay plots with human-readable facet labels (covariate
x class, or correlation pair) instead of raw Stan parameter names –
enough to tell slow-mixing, chain separation, a sticky region, and
multimodality apart, which a single Rhat number cannot.

## Usage

``` r
plot_probit_worst_offender_diagnostics(
  fit,
  draws,
  title_base = "",
  n_worst = 6L
)
```

## Arguments

- fit:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- draws:

  A
  [`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
  (typically `fit$draws`).

- title_base:

  Character. Prefixed to plot titles.

- n_worst:

  Integer. How many parameters to select by each criterion (Rhat, bulk
  ESS) before deduplicating. Default `6`.

## Value

Named list with `trace` and `rank` ggplot objects, or `NULL` if no
structural parameter diagnostics are available.
