# Validate Date Logic

Checks that date sequences and derived values are internally consistent.
Consolidates all cross-column date/age checks in one place: 1.
admission_date \<= culture_date \<= outcome_date 2. dob \<=
admission_date 3. age in \[0, 120\] 4. age vs DOB-computed age within
2-year tolerance 5. "Died" rows must have an outcome date

## Usage

``` r
prep_validate_date_logic(
  data,
  admission_col = "admission_date",
  culture_col = "culture_date",
  outcome_col = "outcome_date",
  dob_col = "dob",
  age_col = "age",
  outcome_value_col = "final_outcome"
)
```

## Arguments

- data:

  Data frame with standard column names.

- admission_col:

  Character. Admission date column. Default "admission_date".

- culture_col:

  Character. Culture/event date column. Default "culture_date".

- outcome_col:

  Character. Outcome date column. Default "outcome_date".

- dob_col:

  Character. Date of birth column. Default "dob".

- age_col:

  Character. Age column. Default "age".

- outcome_value_col:

  Character. Outcome value column for death check. Default
  "final_outcome".

## Value

Invisibly returns data. Prints violation summary.

## Details

Reports violations as warnings; does not drop or modify rows.
