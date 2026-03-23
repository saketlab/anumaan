# Flag Polymicrobial Infections (No Specimen Type, No Episode ID)

Flags polymicrobial status by counting distinct organisms per patient
(optionally within facility/syndrome scope). No specimen type and no
episode_id are used.

## Usage

``` r
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

  Data frame with patient and organism columns.

- patient_col:

  Character. Patient ID column. Default `"patient_id"`.

- organism_col:

  Character. Organism column. Default `"organism_normalized"`.

- facility_col:

  Character or `NULL`. Facility column. When supplied, polymicrobial
  status is scoped within each facility.

- facility_name:

  Character or `NULL`. When supplied with `facility_col`, data are first
  filtered to that facility.

- syndrome_col:

  Character or `NULL`. Syndrome column. When supplied, polymicrobial
  status is scoped within each syndrome.

- syndrome_name:

  Character or `NULL`. When supplied with `syndrome_col`, data are first
  filtered to that syndrome.

## Value

Data frame with `n_organisms` and `is_polymicrobial` (0/1).
