# Near-degeneracy diagnostic for the posterior Omega correlation matrix

A multivariate-probit correlated-residual fit can push Omega toward a
near-singular boundary (very high pairwise correlations across the
board) without that being visible from any single Omega\\i,j\\ entry.
For each posterior draw this computes the smallest eigenvalue of the
full D x D Omega matrix – a value approaching 0 means that specific
draw's correlation matrix is nearly singular – and reports the
distribution across draws, plus the condition number (largest / smallest
eigenvalue) as a supplementary summary.

## Usage

``` r
plot_omega_degeneracy_diagnostic(
  fit,
  class_cols,
  title_base = "",
  degenerate_threshold = 0.05
)
```

## Arguments

- fit:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- class_cols:

  Character vector of all class names, in the fit's canonical order.

- title_base:

  Character. Prefixed to the plot title.

- degenerate_threshold:

  Numeric. Draws with smallest eigenvalue below this are counted as
  "near-degenerate" in the subtitle. Default `0.05`.

## Value

A ggplot object, or `NULL` if this is not a correlated-residual fit, or
Omega draws are unavailable.
