# Human-readable worst-converging structural parameters

Replaces raw Stan parameter names (e.g. `"beta[12,2]"`) with the
substantive covariate/class (or correlation pair) they refer to, and
shows only parameters that actually exceed the Rhat threshold – not an
arbitrary fixed top-N – alongside their bulk and tail ESS: a high Rhat
with healthy ESS is a different diagnostic situation than a high Rhat
with very low ESS, and showing only one of the two obscures that
distinction.

## Usage

``` r
plot_probit_worst_parameters(
  fit,
  title_base = "",
  rhat_threshold = 1.01,
  max_rows = 40L
)
```

## Arguments

- fit:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- title_base:

  Character. Prefixed to plot titles.

- rhat_threshold:

  Numeric. Only parameters with Rhat above this are shown. Default
  `1.01`.

- max_rows:

  Integer cap on rows shown, in case very many parameters fail. Default
  `40`.

## Value

A ggplot object (table page), or `NULL` if no structural parameter
exceeds the threshold, or diagnostics are unavailable.
