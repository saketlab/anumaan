# Assign Age Bins

Categorizes age into bins for stratification using GBD standard or
custom break points. Supports ages stored as a single column in years,
months, or days, as well as compound datasets where years, months, and
days are in separate columns (e.g. neonatal records).

## Usage

``` r
prep_assign_age_bins(
  data,
  age_col = "Age",
  bins = "GBD_standard",
  age_unit = "years",
  age_months_col = NULL,
  age_days_col = NULL,
  negative_age_strategy = c("fallback", "na"),
  fallback_years_col = "year",
  fallback_months_col = "months",
  fallback_days_col = "age_days",
  fallback_dob_col = "dob",
  fallback_admission_col = "admission_date"
)
```

## Arguments

- data:

  Data frame with an age column.

- age_col:

  Character. Primary age column. Default `"Age"`.

- bins:

  Character or character vector. Preset name (`"GBD_standard"`,
  `"pediatric"`, `"geriatric"`, `"neonatal"`) or a custom vector of bin
  labels. Default `"GBD_standard"`.

- age_unit:

  Character. Unit of `age_col`: `"years"` (default), `"months"`, or
  `"days"`. The value is converted to decimal years before binning.

- age_months_col:

  Character or NULL. Optional column holding the months component of age
  (0-11). Added to `age_col` after unit conversion.

- age_days_col:

  Character or NULL. Optional column holding the days component of age
  (0-30). Added to `age_col` after unit conversion.

- negative_age_strategy:

  Character. How to handle negative ages in `age_col`. `"fallback"`
  (default) tries alternate columns / DOB logic and falls back to 0
  years when unresolved; `"na"` sets negative ages to `NA`.

- fallback_years_col:

  Character or NULL. Optional years column used when `age_col` is
  negative. Default `"year"`.

- fallback_months_col:

  Character or NULL. Optional months column used when `age_col` is
  negative. Default `"months"`.

- fallback_days_col:

  Character or NULL. Optional days column used when `age_col` is
  negative. Default `"age_days"`.

- fallback_dob_col:

  Character or NULL. Optional DOB column used with
  `fallback_admission_col` to derive age when `age_col` is negative.
  Default `"dob"`.

- fallback_admission_col:

  Character or NULL. Optional admission/reference date column used with
  `fallback_dob_col`. Default `"admission_date"`.

## Value

Data frame with `Age_bin` factor column added.

## Details

When `age_months_col` or `age_days_col` are supplied they are treated as
the fractional remainder on top of `age_col` (e.g. age_years = 0,
age_months = 3, age_days = 5 -\> 0.27 years). Missing values in the
component columns are treated as zero so that a partially-recorded age
is still binnable.
