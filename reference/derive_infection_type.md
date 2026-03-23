# Derive Infection Type (HAI / CAI) per Patient

Classifies each row as HAI or CAI. Uses `infection_type_col` when it
contains a valid value. For rows where that column is `NA`,
`"Not known"`, or `"NULL"`, derives the classification from the gap
between `date_culture_col` and `date_admission_col`:

- gap \<= `hai_threshold_hours` -\> **CAI**

- gap \> `hai_threshold_hours` -\> **HAI**

## Usage

``` r
derive_infection_type(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
)
```

## Arguments

- data:

  Data frame.

- infection_type_col:

  Character. Raw infection type column. Default `"type_of_infection"`.

- date_admission_col:

  Character. Default `"date_of_admission"`.

- date_culture_col:

  Character. Date of first positive culture. Default
  `"date_of_first_positive_culture"`.

- hai_threshold_hours:

  Numeric. Gap threshold in hours. Default `48`.

- patient_id_col:

  Character. Unique patient identifier column. Default
  `"PatientInformation_id"`.

## Value

`data` with column `infection_type_derived` (`"HAI"` / `"CAI"` /
`"Unknown"`).
