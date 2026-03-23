# GBD Figure 4 – Deaths by Organism Group (Overlapping Bars)

Overlapping bar chart showing deaths associated with vs attributable to
AMR, grouped by organism, replicating GBD Lancet 2022 Figure 4.

## Usage

``` r
plot_gbd_fig4(deaths_df, total_deaths, title = NULL)
```

## Arguments

- deaths_df:

  Data frame with columns: `org_group`, `deaths_associated`,
  `deaths_attributable`.

- total_deaths:

  Numeric. Total deaths for subtitle.

- title:

  Character. Plot title. Default generates a standard title.

## Value

A `ggplot` object.
