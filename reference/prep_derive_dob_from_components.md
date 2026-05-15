# Derive Date of Birth from Age Components

Stewardship-specific helper: DOB is sometimes stored as separate
year/month/day components or as a decimal age in years. This function
reconstructs an approximate DOB or derives age directly when DOB is
absent.

## Usage

``` r
prep_derive_dob_from_components(
  data,
  dob_year_col = "dob_year",
  dob_month_col = "dob_month",
  dob_day_col = "dob_day",
  age_years_col = "age_years",
  reference_date_col = "admission_date",
  dob_output_col = "dob"
)
```

## Arguments

- data:

  Data frame.

- dob_year_col:

  Character. Year component column. Default "dob_year".

- dob_month_col:

  Character. Month component column. Default "dob_month".

- dob_day_col:

  Character. Day component column. Default "dob_day".

- age_years_col:

  Character. Decimal age in years. Default "age_years".

- reference_date_col:

  Character. Reference date for age-based DOB estimation. Default
  "admission_date".

- dob_output_col:

  Character. Output DOB column. Default "dob".

## Value

Data frame with `dob` column added/populated.

## Details

Strategies (applied in order):

1.  If `dob_year_col`, `dob_month_col`, and `dob_day_col` are present,
    assemble into a Date.

2.  If `age_years_col` is present and a reference date exists, estimate
    DOB as `reference_date - age_years * 365.25`.
