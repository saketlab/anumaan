# Prepare Diagnosis Text

Combines a primary diagnosis column with an optional fallback column,
trims whitespace, and collapses repeated separators. Optionally appends
the organism name to the text (as the notebook does for context-aware
embedding).

## Usage

``` r
prep_diagnosis_text(
  data,
  diagnosis_col,
  fallback_col = NULL,
  output_col = "diagnosis_text",
  include_organism = FALSE,
  organism_col = "organism_normalized",
  keep_original = TRUE
)
```

## Arguments

- data:

  Data frame.

- diagnosis_col:

  Character. Primary diagnosis column.

- fallback_col:

  Character or NULL. Used when `diagnosis_col` is NA. Default NULL.

- output_col:

  Character. Name for the output text column. Default
  `"diagnosis_text"`.

- include_organism:

  Logical. Append organism name to the text. Default FALSE.

- organism_col:

  Character. Organism column. Used only when `include_organism = TRUE`.
  Default `"organism_normalized"`.

- keep_original:

  Logical. Keep the original diagnosis columns. Default TRUE.

## Value

Data frame with `output_col` added.
