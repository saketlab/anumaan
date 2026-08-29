# Sampling-diagnostic (Rhat) heatmap for Omega

Companion to
[`plot_omega_correlation_heatmap()`](https://saketlab.github.io/anumaan/reference/plot_omega_correlation_heatmap.md)
– shows Rhat for each `Omega[i,j]` entry, so a reader cannot look at an
attractive correlation heatmap without also seeing whether those
specific parameters actually converged. A uniformly high correlation
heatmap alongside a uniformly poor Rhat heatmap here is the signature of
a near-degenerate, weakly-identified correlation matrix, not a settled
finding.

## Usage

``` r
plot_omega_convergence_heatmap(corr_summary, class_cols, title_base = "")
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

A ggplot object, or `NULL` if `corr_summary` is `NULL`/empty or has no
`rhat` column (e.g. computed by an older
[`summarize_fit_correlation_matrix()`](https://saketlab.github.io/anumaan/reference/summarize_fit_correlation_matrix.md)
call).
