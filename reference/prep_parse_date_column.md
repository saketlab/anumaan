# Detect and decode encrypted or non-standard date columns

Inspects a column for date-like content and automatically selects the
correct parsing strategy (Excel serial, Unix timestamp ms, ISO/DMY/MDY
strings, reversed YYYYDDMM).

## Usage

``` r
prep_parse_date_column(x, col_name = "date", table_label = "table")
```

## Arguments

- x:

  Vector (character, numeric, or Date-like) to parse.

- col_name:

  Character. Column name used in messages.

- table_label:

  Character. Table label used in messages.

## Value

A Date vector of the same length as `x`.
