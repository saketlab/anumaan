# Standardize Outcome Values

Maps various outcome representations to standard values: "Died",
"Discharged", "LAMA", "Unknown"

## Usage

``` r
standardize_outcome(data, col = "final_outcome")
```

## Arguments

- data:

  Data frame containing outcome column

- col:

  Character. Name of outcome column. Default "final_outcome".

## Value

Data frame with standardized outcome values
