# Apply Column Map to Rename Columns

Renames data frame columns using a standard_name = raw_name map. Warns
on unmapped columns; does not drop them.

## Usage

``` r
prep_apply_column_map(data, column_map)
```

## Arguments

- data:

  Data frame.

- column_map:

  Named character vector from
  [`prep_build_column_map()`](https://saketlab.github.io/anumaan/reference/prep_build_column_map.md).

## Value

Data frame with standard column names where mapping existed.
