# Enrich Hospital Department

Infers hospital department from contextual data (specimen type,
diagnosis, patient demographics).

## Usage

``` r
enrich_hospital_department(
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

  Data frame

- department_col:

  Character. Department column. Default "hospital_department".

- specimen_col:

  Character. Specimen type column. Default "specimen_type".

- diagnosis_col:

  Character. Diagnosis column. Default "diagnosis_1".

- age_col:

  Character. Age column. Default "Age".

- overwrite:

  Logical. Recalculate even if present. Default FALSE.

## Value

Data frame with hospital_department enriched

## Examples

``` r
if (FALSE) { # \dontrun{
data_enriched <- enrich_hospital_department(data)
} # }
```
