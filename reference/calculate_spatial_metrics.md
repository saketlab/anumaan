# Calculate AMR Metrics by Geographic Unit

Aggregates AMR data by district or facility and calculates key metrics
including resistance rates, isolate counts, and diversity indices.

## Usage

``` r
calculate_spatial_metrics(
  data,
  geo_level = "district",
  stratify_by = NULL,
  min_isolates = 30
)
```

## Arguments

- data:

  Data frame with processed AMR data

- geo_level:

  Character. "district" or "facility"

- stratify_by:

  Character vector. Additional stratification variables (e.g.,
  "organism_group")

- min_isolates:

  Numeric. Minimum isolates required to calculate rate. Default 30.

## Value

Data frame with geographic AMR metrics
