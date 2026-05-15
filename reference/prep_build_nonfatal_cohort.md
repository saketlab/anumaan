# Build Non-Fatal Cohort

Extracts records for patients who survived, for use in YLD calculation.

## Usage

``` r
prep_build_nonfatal_cohort(
  data,
  outcome_col = "final_outcome",
  survived_value = "Survived"
)
```

## Arguments

- data:

  Analysis-ready data frame.

- outcome_col:

  Character. Final outcome column. Default "final_outcome".

- survived_value:

  Character. Value indicating survival. Default "Survived".

## Value

Data frame with only non-fatal episodes.
