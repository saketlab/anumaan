# Create Choropleth Map

Creates a thematic map showing AMR metrics by geographic area.

## Usage

``` r
create_choropleth_map(
  spatial_obj,
  variable = "resistance_rate",
  title = "AMR Resistance Rate by District",
  palette = "sequential",
  breaks = NULL,
  labels = NULL
)
```

## Arguments

- spatial_obj:

  sf object with metrics

- variable:

  Character. Variable to map

- title:

  Character. Map title

- palette:

  Character. Color palette: "sequential", "diverging", "hotcold".
  Default "sequential".

- breaks:

  Numeric vector. Custom breaks for classification

- labels:

  Character vector. Labels for breaks

## Value

tmap object
