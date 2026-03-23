# Derive Age from Date of Birth

Calculates age in years from date of birth and a reference date.
Implements proxy logic when Age is already present.

## Usage

``` r
derive_age(
  data,
  dob_col = "DOB",
  reference_date_col = "date_of_culture",
  force = FALSE
)
```

## Arguments

- data:

  Data frame

- dob_col:

  Character. Name of DOB column. Default "DOB".

- reference_date_col:

  Character. Reference date column (usually culture date). Default
  "date_of_culture".

- force:

  Logical. If TRUE, recalculates even if Age present. Default FALSE.

## Value

Data frame with Age column added/updated and Age_derived flag
