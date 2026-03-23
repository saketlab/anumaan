# Create Interactive Leaflet Map

Creates an interactive web map of AMR data.

## Usage

``` r
create_interactive_map(
  spatial_obj,
  variable = "resistance_rate",
  popup_vars = NULL,
  title = "AMR Interactive Map",
  palette = "YlOrRd"
)
```

## Arguments

- spatial_obj:

  sf object with metrics

- variable:

  Character. Variable to map

- popup_vars:

  Character vector. Variables to show in popup

- title:

  Character. Map title

- palette:

  Character. Color palette function name (e.g., "YlOrRd")

## Value

leaflet map object
