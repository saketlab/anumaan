# Validate Data Quality

Runs quality checks: minimum row count, required column presence, and
maximum missing-value percentage per column.

## Usage

``` r
validate_data_quality(
  data,
  min_rows = 10,
  max_missing_pct = 50,
  required_cols = c("patient_id", "organism_normalized"),
  stop_on_failure = FALSE
)
```

## Arguments

- data:

  Data frame to validate.

- min_rows:

  Integer. Minimum acceptable number of rows. Default 10.

- max_missing_pct:

  Numeric. Maximum acceptable percent missing per column (0-100).
  Default 50.

- required_cols:

  Character vector. Columns that must be present.

- stop_on_failure:

  Logical. Stop with error on failure. Default FALSE.

## Value

List with quality assessment.
