# Flag Polymicrobial Infections

Flags polymicrobial status by counting distinct organisms per patient
(optionally scoped within facility and/or syndrome).

## Usage

``` r
prep_flag_polymicrobial(
  data,
  patient_col = "patient_id",
  organism_col = "organism_normalized",
  facility_col = NULL,
  facility_name = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)

flag_polymicrobial(
  data,
  patient_col = "patient_id",
  organism_col = "organism_normalized",
  facility_col = NULL,
  facility_name = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL
)
```

## Arguments

- data:

  Data frame.

- patient_col:

  Character. Patient ID column. Default "patient_id".

- organism_col:

  Character. Organism column. Default "organism_normalized".

- facility_col:

  Character or NULL. Facility column for scoped counting.

- facility_name:

  Character or NULL. Filter to this facility before flagging.

- syndrome_col:

  Character or NULL. Syndrome column for scoped counting.

- syndrome_name:

  Character or NULL. Filter to this syndrome before flagging.

## Value

Data frame with `n_organisms` and `is_polymicrobial` (0/1).
