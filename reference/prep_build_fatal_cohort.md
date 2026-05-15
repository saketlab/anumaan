# Build Fatal Cohort

Extracts records for patients who died, for use in YLL calculation.

## Usage

``` r
prep_build_fatal_cohort(
  data,
  outcome_col = "final_outcome",
  died_value = "Died"
)
```

## Arguments

- data:

  Analysis-ready data frame.

- outcome_col:

  Character. Final outcome column. Default "final_outcome".

- died_value:

  Character. Value indicating death. Default "Died".

## Value

Data frame with only fatal episodes.
