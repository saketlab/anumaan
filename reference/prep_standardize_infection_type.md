# Standardize Infection Type Column

Normalizes HAI/CAI variants to a controlled vocabulary: "HAI", "CAI", or
NA.

## Usage

``` r
prep_standardize_infection_type(data, col = "infection_type")
```

## Arguments

- data:

  Data frame.

- col:

  Character. Infection type column. Default "infection_type".

## Value

Data frame with standardized `infection_type` column.
