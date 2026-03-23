# Groom Optional Columns

Cleans and standardizes all optional columns in a single pass. Handles
value normalization, whitespace trimming, and consistency checks.

## Usage

``` r
groom_optional_columns(data, optional_cols = NULL)
```

## Arguments

- data:

  Data frame

- optional_cols:

  Character vector. Optional column names to groom. If NULL, grooms all
  known optional columns. Default NULL.

## Value

Data frame with groomed optional columns

## Examples

``` r
if (FALSE) { # \dontrun{
data_groomed <- groom_optional_columns(data)

# Groom specific columns
data_groomed <- groom_optional_columns(
  data,
  optional_cols = c("hospital_department", "unit_type", "comorbidities")
)
} # }
```
