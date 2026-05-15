# Flag and Classify Readmissions

For each patient, classifies admissions into:

- index:

  First admission for this patient.

- linked_readmission:

  Readmission within `gap_linked_days` – treated as the same episode.
  Excluded from HAI incidence counts.

- new_readmission:

  Readmission between `gap_linked_days` and `gap_new_days`.

- late_readmission:

  Readmission after `gap_new_days` – fully new event.

## Usage

``` r
prep_flag_readmission(
  data,
  patient_col = "patient_id",
  admission_col = "admission_date",
  gap_linked_days = 30,
  gap_new_days = 90
)
```

## Arguments

- data:

  Data frame.

- patient_col:

  Character. Patient ID column. Default "patient_id".

- admission_col:

  Character. Admission date column. Default "admission_date".

- gap_linked_days:

  Numeric. Gap threshold for linked readmissions. Default 30.

- gap_new_days:

  Numeric. Gap threshold for new vs late readmissions. Default 90.

## Value

Data frame with `readmission_class` column added.
