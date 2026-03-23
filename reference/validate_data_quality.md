# Validate Data Quality

Assesses overall data quality metrics and flags datasets that don't meet
minimum quality thresholds.

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

  Data frame

- min_rows:

  Numeric. Minimum number of rows required. Default 10.

- max_missing_pct:

  Numeric. Maximum percent of missing values allowed per column (0-100).
  Default 50.

- required_cols:

  Character vector. Columns that must be present. Default
  c("patient_id", "organism_normalized").

- stop_on_failure:

  Logical. If TRUE, stops on quality failure. Default FALSE.

## Value

List with quality assessment: - passes_quality: Logical - n_rows: Number
of rows - n_cols: Number of columns - overall_completeness: Overall
proportion non-missing (0-1) - column_completeness: Data frame of
per-column completeness - quality_issues: Character vector of identified
issues

## Examples

``` r
if (FALSE) { # \dontrun{
quality <- validate_data_quality(data, min_rows = 100, max_missing_pct = 30)

if (!quality$passes_quality) {
  print(quality$quality_issues)
}
} # }
```
