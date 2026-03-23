# Compute Patient-Level Post-Infection LOS

Computes LOS with infection-type-specific clock start:

- **CAI**: LOS = date_discharge - date_admission

- **HAI**: LOS = date_discharge - date_culture

- **Unknown / Not known**: LOS = date_discharge - date_culture (if
  culture date present), else date_discharge - date_admission. If
  neither reference date is available the patient is excluded.

Only discharged patients are retained. A checkpoint reports how many are
missing a discharge date; those patients fall back to `los_col` when
provided, otherwise they are excluded. Rows with LOS \<= 0 or LOS \>
`max_los` are dropped. Returns one row per patient.

## Usage

``` r
compute_patient_los(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  syndrome_name = NULL,
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  infection_type_derived_col = "infection_type_derived",
  los_col = NULL,
  max_los = 200
)
```

## Arguments

- data:

  Data frame (after
  [`derive_infection_type()`](https://saketlab.github.io/anumaan/reference/derive_infection_type.md)
  has been run).

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- facility_col:

  Character. Default `"center_name"`.

- facility_name:

  Character or `NULL`. If provided, filters data to the specified
  facility before computing LOS. Default `NULL`.

- organism_col:

  Character. Column containing organism/pathogen names. Used for episode
  grouping: a new episode is only created when the culture date gap
  exceeds 14 days AND the organism differs from the episode's first
  organism. Default `"organism_name"`.

- syndrome_col:

  Character. Syndrome column name. Only used when `syndrome_name` is not
  `NULL`. Default `"syndrome"`.

- syndrome_name:

  Character or `NULL`. If provided, only patients with this syndrome are
  retained before LOS computation. Default `NULL`.

- date_admission_col:

  Character. Default `"date_of_admission"`.

- date_discharge_col:

  Character. Default `"final_outcome_date"`.

- date_culture_col:

  Character. Default `"date_of_first_positive_culture"`.

- final_outcome_col:

  Character. Default `"final_outcome"`.

- final_outcome_value:

  Character. Default `"Discharged"`.

- infection_type_derived_col:

  Character. Column from
  [`derive_infection_type()`](https://saketlab.github.io/anumaan/reference/derive_infection_type.md).
  Default `"infection_type_derived"`.

- los_col:

  Character or `NULL`. Optional name of a pre-computed LOS column (e.g.
  recorded directly in the data). Used only for patients whose discharge
  date is missing. Values that are `NA`, zero, or negative are treated
  as invalid and those patients are excluded. Default `NULL`.

- max_los:

  Numeric. Upper cap on LOS in days; patients exceeding this are
  excluded. Default `200`.

## Value

Data frame: one row per patient with `patient_id_col`, `facility_col`,
`infection_type_derived_col`, `LOS_days`.
