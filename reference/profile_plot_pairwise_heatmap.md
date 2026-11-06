# Pairwise Co-Resistance Heatmap, One Hospital

Class x class heatmap of the model's pairwise P(R,R) for one hospital,
mirroring
[`plot_omega_correlation_heatmap`](https://saketlab.github.io/anumaan/reference/plot_omega_correlation_heatmap.md)'s
exact visual style. This is a joint-Bernoulli pairwise probability, not
a latent residual correlation matrix – a different quantity that happens
to render the same way (classes x classes, symmetric, one cell per
pair).

## Usage

``` r
profile_plot_pairwise_heatmap(ppc_pairwise, group, class_cols, title_base = "")
```

## Arguments

- ppc_pairwise:

  [`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md)'s
  `$pairwise` table.

- group:

  Character. Which hospital to plot.

- class_cols:

  Character vector, all classes in the fit's canonical order (for
  consistent row/column ordering).

- title_base:

  Character. Prefixed to the plot title.

## Value

A ggplot object, or `NULL` if this group has no evaluated pairwise
cells.
