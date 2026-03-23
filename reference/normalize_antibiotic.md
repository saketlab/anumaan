# Normalize Antibiotic Names

Standardizes antibiotic names using fuzzy matching against WHO
reference. Similar approach to organism normalization.

## Usage

``` r
normalize_antibiotic(
  data,
  antibiotic_col = "antibiotic_name",
  who_table = NULL,
  add_class = TRUE,
  add_aware = TRUE
)
```

## Arguments

- data:

  Data frame with antibiotic data

- antibiotic_col:

  Character. Column name with antibiotic names. Default
  "antibiotic_name".

- who_table:

  Data frame. WHO AWaRe classification table. If NULL, loads from
  inst/extdata/WHO_aware_class.csv.

- add_class:

  Logical. Add antibiotic_class column. Default TRUE.

- add_aware:

  Logical. Add aware_category column. Default TRUE.

## Value

Data frame with antibiotic_normalized, antibiotic_class, aware_category

## Examples

``` r
if (FALSE) { # \dontrun{
data <- normalize_antibiotic(data, antibiotic_col = "antibiotic_name")
} # }
```
