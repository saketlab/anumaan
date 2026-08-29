# Full per-pair Omega summary table (rho, 95% CrI, Rhat, ESS bulk/tail)

The two heatmaps show the correlation estimate and its Rhat as separate
pictures; this page puts every number for every class pair – median rho,
credible interval, Rhat, bulk and tail ESS – in one table, for a reader
who wants the exact figures rather than reading them off a colour scale.

## Usage

``` r
plot_omega_summary_table(corr_summary, class_cols, title_base = "")
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

## Value

A ggplot object (table page), or `NULL` if `corr_summary` is
`NULL`/empty.
