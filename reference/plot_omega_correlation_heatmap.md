# Posterior-median latent cross-class correlation heatmap (Omega)

Rows/columns are antimicrobial classes; each cell is the posterior
median of the corresponding `Omega[i,j]` (or `R_block[r,i,j]`) entry.
Only meaningful for correlated-residual fits – pass `NULL` (e.g. an
identity-residual fit's
[`summarize_fit_correlation_matrix()`](https://saketlab.github.io/anumaan/reference/summarize_fit_correlation_matrix.md)
result) to get `NULL` back rather than an empty/misleading plot.

## Usage

``` r
plot_omega_correlation_heatmap(
  corr_summary,
  class_cols,
  title_base = "",
  rhat_flag_threshold = 1.01
)
```

## Arguments

- corr_summary:

  Tibble from
  [`summarize_fit_correlation_matrix()`](https://saketlab.github.io/anumaan/reference/summarize_fit_correlation_matrix.md)
  (with `class_1, class_2, correlation_median, rhat` columns), or
  `NULL`.

- class_cols:

  Character vector of all class names, in the fit's canonical order (for
  consistent row/column ordering).

- title_base:

  Character. Prefixed to the plot title.

- rhat_flag_threshold:

  Numeric. Cells with `rhat` above this are marked with an asterisk,
  cross-referencing the companion convergence heatmap rather than
  repeating full diagnostics on this plot.

## Value

A ggplot object, or `NULL` if `corr_summary` is `NULL`/empty.
