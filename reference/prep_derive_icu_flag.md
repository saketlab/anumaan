# Derive ICU Flag

Sets `icu_flag = TRUE` where ward/department information indicates the
patient was in an ICU at the time of culture.

## Usage

``` r
prep_derive_icu_flag(
  data,
  ward_col = "ward_icu",
  dept_col = "hospital_department"
)
```

## Arguments

- data:

  Data frame.

- ward_col:

  Character. Ward/unit column. Default "ward_icu".

- dept_col:

  Character. Department column. Default "hospital_department".

## Value

Data frame with `icu_flag` logical column added.
