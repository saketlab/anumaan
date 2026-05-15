# Map ICD Codes to Infectious Syndromes

Joins an ICD-to-syndrome reference table onto the data. The reference
must have at minimum two columns: an ICD code/description column and a
syndrome column. Typically this is a manually curated mapping CSV.

## Usage

``` r
prep_map_icd_to_syndrome(
  data,
  icd_col = "icd_code",
  icd_to_syndrome_ref,
  ref_icd_col = "icd_code",
  ref_syndrome_col = "infectious_syndrome",
  output_col = "syndrome",
  unmatched_label = NA_character_
)
```

## Arguments

- data:

  Data frame (output of
  [`prep_map_diagnosis_to_icd()`](https://saketlab.github.io/anumaan/reference/prep_map_diagnosis_to_icd.md)
  or any data frame with an ICD column).

- icd_col:

  Character. Column in `data` holding ICD codes or descriptions to match
  on. Default `"icd_code"`.

- icd_to_syndrome_ref:

  Data frame. Reference table with at minimum `icd_col` values and a
  syndrome column.

- ref_icd_col:

  Character. The ICD column name in `icd_to_syndrome_ref`. Default
  `"icd_code"`.

- ref_syndrome_col:

  Character. The syndrome column name in the reference. Default
  `"infectious_syndrome"`.

- output_col:

  Character. Name for the syndrome column added to `data`. Default
  `"syndrome"`.

- unmatched_label:

  Character. Label for rows with no syndrome match. Default
  `NA_character_`.

## Value

Data frame with `output_col` added.
