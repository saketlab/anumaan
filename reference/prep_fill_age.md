# Fill Missing Age Values

Fills missing Age by computing from DOB and a reference date. Tracks
derivation method and confidence for each record.

## Usage

``` r
prep_fill_age(
  data,
  age_col = "Age",
  dob_col = "DOB",
  date_col = "date_of_culture",
  overwrite = FALSE
)
```

## Arguments

- data:

  Data frame.

- age_col:

  Character. Age column. Default "Age".

- dob_col:

  Character. Date of birth column. Default "DOB".

- date_col:

  Character. Reference date column. Default "date_of_culture".

- overwrite:

  Logical. Recalculate even when Age is present. Default FALSE.

## Value

Data frame with Age enriched plus `age_method` and `age_confidence`
columns.
