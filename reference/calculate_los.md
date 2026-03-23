# Calculate Length of Stay

Calculates hospital length of stay in days.

## Usage

``` r
calculate_los(
  data,
  admission_col = "date_of_admission",
  outcome_col = "date_of_final_outcome"
)
```

## Arguments

- data:

  Data frame

- admission_col:

  Character. Admission date column.

- outcome_col:

  Character. Outcome/discharge date column.

## Value

Data frame with Length_of_stay column
