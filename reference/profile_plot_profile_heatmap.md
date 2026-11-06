# Hospital x Profile-State Heatmap

Visualizes the estimated or observed resistance-profile distribution
across hospitals and states in one grid – columns are hospitals, rows
are states (profiles), fill is the selected probability column. Mirrors
[`plot_omega_correlation_heatmap`](https://saketlab.github.io/anumaan/reference/plot_omega_correlation_heatmap.md)'s
visual style.

## Usage

``` r
profile_plot_profile_heatmap(
  profile_table,
  value_col,
  title_base = "",
  max_states = 16L
)
```

## Arguments

- profile_table:

  Either
  [`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md)'s
  `$complete_profile` table or
  [`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)'s
  output.

- value_col:

  Character. Which column to plot as the fill – `"model_mean"` or
  `"empirical_frequency"` for a `profile_ppc()$complete_profile` table,
  or `"R_mean"` for a
  [`profile_aggregate_for_daly`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)
  table.

- title_base:

  Character. Prefixed to the plot title.

- max_states:

  Integer or `NULL`. Keeps only the `max_states` states with the largest
  mean value across hospitals (for readability when `D` is large).
  Default `16L`.

## Value

A ggplot object, or `NULL` if nothing to plot.
