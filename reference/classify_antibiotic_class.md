# Classify Antibiotic to WHO Class

Maps antibiotic names to WHO antibiotic classes.

## Usage

``` r
classify_antibiotic_class(
  data,
  antibiotic_col = "antibiotic_normalized",
  who_table = NULL
)
```

## Arguments

- data:

  Data frame

- antibiotic_col:

  Character. Normalized antibiotic column. Default
  "antibiotic_normalized".

- who_table:

  Data frame. WHO classification table. If NULL, uses built-in mapping.

## Value

Data frame with antibiotic_class column added
