# Calculate Spatial Autocorrelation (Moran's I)

Tests for spatial clustering of AMR rates using Moran's I statistic.
Requires sf object with spatial geometries.

## Usage

``` r
calculate_spatial_autocorrelation(
  spatial_obj,
  variable = "resistance_rate",
  neighbors_k = 4
)
```

## Arguments

- spatial_obj:

  sf object from create_spatial_object()

- variable:

  Character. Variable name to test for autocorrelation

- neighbors_k:

  Integer. Number of nearest neighbors. Default 4.

## Value

List with Moran's I statistic, p-value, and interpretation
