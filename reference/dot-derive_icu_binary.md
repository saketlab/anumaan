# Derive ICU Binary Flag per Patient

Collapses unit-type data (one row per drug test per patient) to a
patient-level binary ICU indicator using the **"ever in ICU"** rule: if
a patient has *any* row with an ICU-type location during that admission,
`ICU = 1`; otherwise `ICU = 0`.

Collapses unit-type data to a patient-level binary ICU indicator using
the "ever in ICU" rule. Returns integer 0/1 or a 3-level factor when
missing data is high.

## Usage

``` r
.derive_icu_binary(
  data,
  patient_id_col = "PatientInformation_id",
  unit_type_col = "unit_type",
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  missing_threshold = 0.1
)

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

  Data frame at drug-test level.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- unit_type_col:

  Character. Default `"unit_type"`.

- icu_values:

  Character vector. ICU location labels.

- missing_threshold:

  Numeric. Proportion above which a 3-level factor is returned. Default
  `0.10`.

## Value

Patient-level data frame: `patient_id_col` and `ICU` (integer 0/1, or
factor `"ICU"/"Ward"/"Unknown"` when missing is high).

Patient-level data frame with `patient_id_col` and `ICU`.

## Details

Missing unit-type values are coded as `NA`. When the proportion of
patients with *all* unit-type rows missing exceeds `missing_threshold`,
the function returns a three-level ordered factor (`"ICU"` / `"Ward"` /
`"Unknown"`) and emits a message recommending that the `"Unknown"` level
be included in the model.
