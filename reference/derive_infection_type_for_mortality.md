# Derive Infection Type (HAI / CAI) for Mortality RR Model

Variant of
[`derive_infection_type()`](https://saketlab.github.io/anumaan/reference/derive_infection_type.md)
tailored for the mortality relative-risk model. Key differences from the
LOS version:

1.  Processes **all** patients (not just discharged), because the
    mortality model requires both dead and surviving patients.

2.  Recognises explicit HAI/CAI labels first (and common synonyms such
    as "hospital-acquired infection" / "community acquired infection").
    Date-gap derivation is applied only when the label is absent or
    ambiguous.

3.  Performs a dedicated death-date data-quality check: patients whose
    `final_outcome` equals `death_value` but whose outcome date is
    missing – yet admission or culture date is present – are listed.
    HAI/CAI derivation is unaffected (it uses admission vs culture
    date), but the missing death date is surfaced for review.

Variant of
[`derive_infection_type()`](https://saketlab.github.io/anumaan/reference/derive_infection_type.md)
tailored for the mortality relative-risk model. Processes all patients
(dead and surviving), recognises explicit HAI/CAI labels first, and
falls back to date-gap derivation when labels are ambiguous. Also
performs a death-date data-quality check.

## Usage

``` r
derive_infection_type_for_mortality(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  death_value = "Death",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
)

derive_infection_type_for_mortality(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  death_value = "Death",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
)
```

## Arguments

- data:

  Data frame.

- infection_type_col:

  Character. Default `"type_of_infection"`.

- date_admission_col:

  Character. Default `"date_of_admission"`.

- date_culture_col:

  Character. Default `"date_of_first_positive_culture"`.

- final_outcome_col:

  Character. Default `"final_outcome"`.

- final_outcome_date_col:

  Character. Default `"final_outcome_date"`.

- death_value:

  Character. Default `"Death"`.

- hai_threshold_hours:

  Numeric. Default `48`.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

## Value

`data` with column `infection_type_derived` (`"HAI"` / `"CAI"` /
`"Not Known"`). Attribute `"missing_death_date_patients"` is a data
frame of rows where death is confirmed but the outcome date is absent.

`data` with column `infection_type_derived` (`"HAI"` / `"CAI"` /
`"Not Known"`).
