# Convert Wide Format to Long Format

Converts wide format data (where each antibiotic is a column) to long
format (one row per organism-antibiotic combination). This is the first
step before normalization and analysis.

## Usage

``` r
pivot_wide_to_long(
  data,
  antibiotic_cols = NULL,
  pattern = NULL,
  id_cols = c("patient_id", "event_id", "organism_name", "date_of_culture"),
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  remove_missing = TRUE,
  create_event_id = FALSE
)
```

## Arguments

- data:

  Data frame in wide format (antibiotics as columns)

- antibiotic_cols:

  Character vector. Names of antibiotic columns to pivot. If NULL, will
  auto-detect based on pattern. Default NULL.

- pattern:

  Character. Regex pattern to identify antibiotic columns if
  antibiotic_cols not provided. Default NULL (no auto-detection).

- id_cols:

  Character vector. Columns to keep as identifiers (not pivoted).
  Default c("patient_id", "event_id", "organism_name",
  "date_of_culture").

- antibiotic_name_col:

  Character. Name for the new column containing antibiotic names.
  Default "antibiotic_name".

- antibiotic_value_col:

  Character. Name for the new column containing susceptibility results.
  Default "antibiotic_value".

- remove_missing:

  Logical. Remove rows where antibiotic_value is NA, empty, or "-".
  Default TRUE.

- create_event_id:

  Logical. Create event_id column if it doesn't exist (uses row
  numbers). Default FALSE.

## Value

Data frame in long format

## Examples

``` r
if (FALSE) { # \dontrun{
# Specify antibiotic columns explicitly
long_data <- pivot_wide_to_long(
  data = raw_data,
  antibiotic_cols = c("AMIKACIN", "GENTAMICIN", "CIPROFLOXACIN")
)

# Auto-detect columns by pattern (columns 12-53)
long_data <- pivot_wide_to_long(
  data = raw_data,
  antibiotic_cols = names(raw_data)[12:53]
)

# Auto-detect uppercase antibiotic names
long_data <- pivot_wide_to_long(
  data = raw_data,
  pattern = "^[A-Z]+$"
)
} # }
```
