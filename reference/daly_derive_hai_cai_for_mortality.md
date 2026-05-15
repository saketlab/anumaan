# Classify HAI/CAI Infection Type for the Mortality Cohort

Derives the healthcare-associated infection (HAI) vs community-acquired
infection (CAI) classification for each patient record in the mortality
cohort, using admission-to-culture date differences. Also performs a
data-quality check for deceased patients missing an outcome date.

## Usage

``` r
daly_derive_hai_cai_for_mortality(
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

  Data frame of patient-level records.

- infection_type_col:

  Character. Column containing the raw infection type label. Default
  `"type_of_infection"`.

- date_admission_col:

  Character. Admission date column. Default `"date_of_admission"`.

- date_culture_col:

  Character. First positive culture date column. Default
  `"date_of_first_positive_culture"`.

- final_outcome_col:

  Character. Final outcome column. Default `"final_outcome"`.

- final_outcome_date_col:

  Character. Final outcome date column. Default `"final_outcome_date"`.

- death_value:

  Character. Value in `final_outcome_col` indicating a fatal outcome.
  Default `"Death"`.

- hai_threshold_hours:

  Numeric. Hours after admission at which a positive culture is
  classified as HAI. Default `48`.

- patient_id_col:

  Character. Patient identifier column. Default
  `"PatientInformation_id"`.

## Value

The input data frame with a derived `infection_type` column (values
`"HAI"` or `"CAI"`).
