# Prepare LOS Dataset

Calculates length of stay from admission and outcome dates, filters to
discharged patients with valid LOS within a plausible range.

## Usage

``` r
prepare_los_data(
  data,
  admission_col = "date_of_admission",
  outcome_date_col = "final_outcome_date",
  outcome_col = "final_outcome",
  patient_id_col = "PatientInformation_id",
  max_los = 200
)
```

## Arguments

- data:

  Data frame with patient-level records.

- admission_col:

  Character. Column name for admission date. Default
  `"date_of_admission"`.

- outcome_date_col:

  Character. Column name for outcome/discharge date. Default
  `"final_outcome_date"`.

- outcome_col:

  Character. Column name for outcome status. Default `"final_outcome"`.

- patient_id_col:

  Character. Column name for patient identifier. Default
  `"PatientInformation_id"`.

- max_los:

  Numeric. Maximum plausible LOS in days. Default 200.

## Value

Data frame with one row per patient-admission, including `LOS_days`.
