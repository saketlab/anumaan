# Check required columns exist and report types

Runs before any table is merged. Reports missing columns, type
mismatches, and blank/duplicate column names. Stops on hard failures;
warns on soft ones.

## Usage

``` r
prep_check_columns(
  data,
  required = character(),
  expected_types = character(),
  table_label = "table",
  stop_on_missing = TRUE
)
```

## Arguments

- data:

  Data frame to inspect.

- required:

  Character vector of column names that must be present.

- expected_types:

  Named character vector mapping column names to expected R classes
  (e.g. `c(PatientInformation_id = "character")`). Only checked when the
  column exists.

- table_label:

  Character. Label used in messages (e.g. "patient_KGMU").

- stop_on_missing:

  Logical. Stop if required columns are absent (default TRUE).

## Value

Invisibly returns a tibble summarising every checked column.
