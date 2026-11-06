# Compute Patient-Level Length of Stay

Computes hospital length of stay with a single definition for every
patient, regardless of infection type (HAI/CAI) or culture date: \$\$LOS
= final\\ outcome\\ (discharge)\\ date - admission\\ date\$\$ Only
patients whose final outcome equals `final_outcome_value` (default
`"Discharged"`) are retained. LOS is `NA` when either date is missing;
such patients are excluded. A checkpoint reports how many discharged
patients are missing a discharge date. Rows with LOS \<= 0 or LOS \>
`max_los` are dropped.

## Usage

``` r
daly_compute_patient_los(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
  syndrome_col = "syndrome",
  syndrome_name = NULL,
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  max_los = 200
)
```

## Arguments

- data:

  Data frame.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- facility_col:

  Character. Default `"center_name"`.

- facility_name:

  Character or `NULL`. If provided, filters data to the specified
  facility before computing LOS. Default `NULL`.

- syndrome_col:

  Character. Syndrome column name. Only used when `syndrome_name` is not
  `NULL`. Default `"syndrome"`.

- syndrome_name:

  Character or `NULL`. If provided, only patients with this syndrome are
  retained before LOS computation. Default `NULL`.

- date_admission_col:

  Character. Default `"date_of_admission"`.

- date_discharge_col:

  Character. Final outcome / discharge date column. Default
  `"final_outcome_date"`.

- final_outcome_col:

  Character. Default `"final_outcome"`.

- final_outcome_value:

  Character. Default `"Discharged"`.

- max_los:

  Numeric. Upper cap on LOS in days; patients exceeding this are
  excluded. Default `200`.

## Value

Data frame with one row per distinct patient x facility x LOS
(`patient_id_col`, `facility_col`, `LOS_days`).
