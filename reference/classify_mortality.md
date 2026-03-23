# Classify Infection-Related Mortality

Determines if death was related to infection using date window logic.
Implements proxy logic when dates are missing.

## Usage

``` r
classify_mortality(
  data,
  outcome_col = "final_outcome",
  event_date_col = "date_of_culture",
  outcome_date_col = "date_of_final_outcome",
  window = 14
)
```

## Arguments

- data:

  Data frame

- outcome_col:

  Character. Outcome column. Default "final_outcome".

- event_date_col:

  Character. Event/culture date. Default "date_of_culture".

- outcome_date_col:

  Character. Outcome date. Default "date_of_final_outcome".

- window:

  Numeric. Days after event to classify death as infection-related.
  Default 14.

## Value

Data frame with mortality_infection, mortality_method,
mortality_confidence columns
