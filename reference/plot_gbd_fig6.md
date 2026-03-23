# GBD Figure 6 – Attributable Deaths Heatmap (Pathogen x Drug Class)

Heatmap of deaths attributable to AMR by pathogen and antibiotic class,
replicating GBD Lancet 2022 Figure 6. Includes an "All pathogens"
summary row.

## Usage

``` r
plot_gbd_fig6(profile_df, total_deaths, total_attributable, title = NULL)
```

## Arguments

- profile_df:

  Data frame with columns: `pathogen`, `profile` (drug class),
  `deaths_attributable`.

- total_deaths:

  Numeric. Total deaths for subtitle.

- total_attributable:

  Numeric. Total attributable deaths for subtitle.

- title:

  Character. Plot title.

## Value

A `ggplot` object.
