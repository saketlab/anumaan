# Parse Dates Safely

Attempts to parse date columns using multiple common formats via
lubridate. Returns Date objects or NA for unparseable values.

## Usage

``` r
parse_dates(
  data,
  date_columns = c("date_of_admission", "date_of_culture", "date_of_final_outcome",
    "DOB")
)
```

## Arguments

- data:

  Data frame containing date columns

- date_columns:

  Character vector of column names to parse as dates. Default includes
  common date fields.

## Value

Data frame with date columns converted to Date class
