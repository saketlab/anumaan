# Assign Syndrome from Wide-Format Syndrome Flags

Converts a wide-format data frame (one column per syndrome, values
indicating presence) to long format and assigns the highest-priority
syndrome per patient using the infectious syndrome hierarchy.

## Usage

``` r
infer_patient_syndrome_wide(
  data,
  patient_col = "patient_id",
  syndrome_cols = NULL,
  positive_values = c(1, "1", TRUE, "TRUE", "True", "true", "Yes", "YES", "yes"),
  collapse_unspecified_respiratory = TRUE,
  keep_only_burden_syndromes = FALSE
)
```

## Arguments

- data:

  Data frame in wide format where syndrome columns contain presence
  indicators (`1`/`TRUE`/`"Yes"` etc.).

- patient_col:

  Character. Patient ID column. Default `"patient_id"`.

- syndrome_cols:

  Character vector. Columns to treat as syndrome flags. If NULL, all
  columns except `patient_col` are used. Default NULL.

- positive_values:

  Character/logical/numeric vector. Values treated as "syndrome
  present". Default covers common representations of TRUE/1/Yes.

- collapse_unspecified_respiratory:

  Logical. Passed to
  [`prep_assign_patient_syndrome()`](https://saketlab.github.io/anumaan/reference/prep_assign_patient_syndrome.md).
  Default TRUE.

- keep_only_burden_syndromes:

  Logical. Passed to
  [`prep_assign_patient_syndrome()`](https://saketlab.github.io/anumaan/reference/prep_assign_patient_syndrome.md).
  Default FALSE.

## Value

Data frame with one row per patient and the selected syndrome.

## Details

This is a convenience wrapper around
[`prep_assign_patient_syndrome()`](https://saketlab.github.io/anumaan/reference/prep_assign_patient_syndrome.md).
All syndrome-selection logic and parameters (hierarchy, burden filter,
respiratory collapsing) are handled there.
