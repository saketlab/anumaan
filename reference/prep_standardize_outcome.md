# Standardize Outcome Values (deprecated)

Deprecated. Use
[`prep_standardize_final_outcome()`](https://saketlab.github.io/anumaan/reference/prep_standardize_final_outcome.md)
instead, which adds `outcome_std` ("Survived"/"Died"/NA) without
overwriting the original column.

## Usage

``` r
prep_standardize_outcome(data, col = "final_outcome")
```

## Arguments

- data:

  Data frame containing outcome column.

- col:

  Character. Name of outcome column. Default "final_outcome".

## Value

Data frame with `outcome_std` column added.
