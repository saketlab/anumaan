# Create Spatial Object from AMR Data

Converts AMR metrics to spatial (sf) object using coordinates. Requires
geocoded locations or joins with spatial reference data.

## Usage

``` r
create_spatial_object(metrics_data, geo_level = "district", coords_data = NULL)
```

## Arguments

- metrics_data:

  Data frame from calculate_spatial_metrics()

- geo_level:

  Character. "district" or "facility"

- coords_data:

  Data frame with columns: name, longitude, latitude

## Value

sf object with spatial geometries
