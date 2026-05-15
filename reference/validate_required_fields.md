# Validate required fields exist and meet completeness threshold

Validate required fields exist and meet completeness threshold

## Usage

``` r
validate_required_fields(
  data,
  required_cols,
  stop_on_failure = TRUE,
  allow_na = FALSE,
  min_completeness = 0.8
)
```

## Arguments

- data:

  Data frame.

- required_cols:

  Character vector of required column names.

- stop_on_failure:

  Logical. Stop if validation fails. Default TRUE.

- allow_na:

  Logical. Skip completeness checks if TRUE. Default FALSE.

- min_completeness:

  Numeric (0-1). Minimum completeness required. Default 0.8.

## Value

List with validation result details.
