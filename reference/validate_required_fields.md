# Validate Required Fields

Checks that all required columns are present and contain valid data.
Returns validation results and optionally stops on failure.

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

  Data frame to validate

- required_cols:

  Character vector. Required column names.

- stop_on_failure:

  Logical. If TRUE, stops execution on validation failure. If FALSE,
  returns validation report. Default TRUE.

- allow_na:

  Logical. If TRUE, allows NA values in required columns. Default FALSE.

- min_completeness:

  Numeric. Minimum proportion of non-NA values required (0-1). Default
  0.8 (80 percent completeness).

## Value

List with validation results:

- `valid`: Logical. Overall validation status

- `missing_cols`: Character vector of missing columns

- `incomplete_cols`: Data frame of columns below min_completeness

- `messages`: Character vector of validation messages

## Examples

``` r
if (FALSE) { # \dontrun{
# Stop on failure (default)
validate_required_fields(
  data,
  required_cols = c("patient_id", "organism_normalized", "antibiotic_normalized")
)

# Get validation report without stopping
validation <- validate_required_fields(
  data,
  required_cols = c("patient_id", "organism_normalized"),
  stop_on_failure = FALSE
)
} # }
```
