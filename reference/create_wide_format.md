# Create Wide Format Dataset

Converts long format (one row per organism-antibiotic) to wide format
(one row per event with antibiotic columns). Useful for analysis and
machine learning applications.

## Usage

``` r
create_wide_format(
  data,
  event_col = "event_id",
  antibiotic_col = "antibiotic_normalized",
  susceptibility_col = "antibiotic_value",
  prefix = "abx_",
  keep_cols = c("patient_id", "organism_normalized", "date_of_culture")
)
```

## Arguments

- data:

  Data frame in long format

- event_col:

  Character. Event ID column. Default "event_id".

- antibiotic_col:

  Character. Antibiotic/class column to pivot. Default
  "antibiotic_normalized".

- susceptibility_col:

  Character. Susceptibility column. Default "antibiotic_value".

- prefix:

  Character. Prefix for pivoted columns. Default "abx\_".

- keep_cols:

  Character vector. Additional columns to keep from original data.
  Default c("patient_id", "organism_normalized", "date_of_culture").

## Value

Wide format data frame (one row per event)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic wide format
wide_data <- create_wide_format(data)

# Class-level wide format
wide_data <- create_wide_format(
  data,
  antibiotic_col = "antibiotic_class",
  prefix = "class_"
)
} # }
```
