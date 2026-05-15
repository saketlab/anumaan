# Normalize Specimen/Sample Type

Normalizes specimen/sample type names and adds sample_category and
sterile_classification from the reference CSV file. Includes rule-based
text cleaning and fuzzy matching for common shorthand and minor
misspellings.

## Usage

``` r
prep_standardize_specimens(
  data,
  specimen_col = "specimen_type",
  add_categories = TRUE
)
```

## Arguments

- data:

  Data frame with specimen column.

- specimen_col:

  Character. Specimen column name. Default "specimen_type".

- add_categories:

  Logical. Add sample_category and sterile_classification. Default TRUE.

## Value

Data frame with specimen_normalized, sample_category, and
sterile_classification columns.
