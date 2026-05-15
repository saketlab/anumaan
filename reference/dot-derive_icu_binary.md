# Derive ICU Binary Flag per Patient

Collapses unit-type data (one row per drug test per patient) to a
patient-level binary ICU indicator using the **"ever in ICU"** rule: if
a patient has *any* row with an ICU-type location during that admission,
`ICU = 1`; otherwise `ICU = 0`.

## Usage

``` r
.derive_icu_binary(
  data,
  patient_id_col = "PatientInformation_id",
  unit_type_col = "unit_type",
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  missing_threshold = 0.1
)
```

## Arguments

- data:

  Data frame at drug-test level (before patient collapse).

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- unit_type_col:

  Character. Column containing ward/ICU labels. Default `"unit_type"`.

- icu_values:

  Character vector. Values (matched case-insensitively) that indicate
  ICU admission. Default
  `c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU")`.

- missing_threshold:

  Numeric in \[0, 1\]. Proportion of patients with entirely missing
  unit-type above which a 3-level factor is returned. Default `0.10`.

## Value

Patient-level data frame: `patient_id_col` and `ICU` (integer 0/1, or
factor `"ICU"/"Ward"/"Unknown"` when missing is high).

## Details

Missing unit-type values are coded as `NA`. When the proportion of
patients with *all* unit-type rows missing exceeds `missing_threshold`,
the function returns a three-level ordered factor (`"ICU"` / `"Ward"` /
`"Unknown"`) and emits a message recommending that the `"Unknown"` level
be included in the model.
