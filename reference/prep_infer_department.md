# Infer Hospital Department

Heuristically infers hospital department from specimen type, diagnosis,
and patient age. Low-confidence inference; used only to fill missing
values.

## Usage

``` r
prep_infer_department(
  data,
  department_col = "hospital_department",
  specimen_col = "specimen_type",
  diagnosis_col = "diagnosis_1",
  age_col = "Age",
  overwrite = FALSE
)
```

## Arguments

- data:

  Data frame.

- department_col:

  Character. Department column. Default "hospital_department".

- specimen_col:

  Character. Specimen type column. Default "specimen_type".

- diagnosis_col:

  Character. Diagnosis column. Default "diagnosis_1".

- age_col:

  Character. Age column. Default "Age".

- overwrite:

  Logical. Overwrite existing values. Default FALSE.

## Value

Data frame with `hospital_department` enriched.
