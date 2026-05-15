# Detect and convert all date-like columns in a table

Auto-detects columns with date/time-like names (or uses a supplied
list), then calls
[`prep_parse_date_column()`](https://saketlab.github.io/anumaan/reference/prep_parse_date_column.md)
on each.

## Usage

``` r
prep_coerce_dates(data, cols = NULL, table_label = "table")
```

## Arguments

- data:

  Data frame.

- cols:

  Character vector of column names to convert. When `NULL`, columns
  whose names match a date/time pattern are detected automatically.

- table_label:

  Character. Label used in messages.

## Value

Data frame with date columns converted to `Date`.
