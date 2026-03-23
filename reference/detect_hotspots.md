# Detect Spatial Hotspots (Getis-Ord Gi\*)

Identifies statistically significant hotspots and coldspots of high/low
AMR.

## Usage

``` r
detect_hotspots(
  spatial_obj,
  variable = "resistance_rate",
  neighbors_k = 4,
  alpha = 0.05
)
```

## Arguments

- spatial_obj:

  sf object from create_spatial_object()

- variable:

  Character. Variable name to analyze

- neighbors_k:

  Integer. Number of nearest neighbors. Default 4.

- alpha:

  Numeric. Significance level. Default 0.05.

## Value

sf object with hotspot classification
