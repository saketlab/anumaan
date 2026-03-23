# Classify AWaRe Category

Assigns WHO AWaRe (Access, Watch, Reserve) categories to antibiotics.

## Usage

``` r
classify_aware(
  data,
  antibiotic_col = "antibiotic_normalized",
  who_table = NULL
)
```

## Arguments

- data:

  Data frame

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_normalized".

- who_table:

  Data frame. WHO AWaRe table. If NULL, uses built-in.

## Value

Data frame with aware_category column added
