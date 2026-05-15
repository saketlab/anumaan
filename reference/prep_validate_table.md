# Run all pre-join sanity checks for one table

Convenience wrapper that calls
[`prep_check_columns()`](https://saketlab.github.io/anumaan/reference/prep_check_columns.md),
[`prep_check_keys()`](https://saketlab.github.io/anumaan/reference/prep_check_keys.md),
and
[`prep_coerce_dates()`](https://saketlab.github.io/anumaan/reference/prep_coerce_dates.md)
in sequence. Returns the date-coerced data along with a check report.

## Usage

``` r
prep_validate_table(
  data,
  required_cols = character(),
  key_col = NULL,
  expected_types = character(),
  date_cols = NULL,
  table_label = "table",
  stop_on_missing = TRUE
)
```

## Arguments

- data:

  Data frame.

- required_cols:

  Character vector of required column names.

- key_col:

  Character. Primary/foreign key column to quality-check.

- expected_types:

  Named character vector of expected column classes.

- date_cols:

  Character vector of date column names to coerce. NULL triggers
  auto-detection.

- table_label:

  Character. Label used in all messages.

- stop_on_missing:

  Logical. Passed to
  [`prep_check_columns()`](https://saketlab.github.io/anumaan/reference/prep_check_columns.md).

## Value

A list with `data` (date-coerced) and `report` (check summary).
