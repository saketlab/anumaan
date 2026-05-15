# Missingness Report

Generates a per-column summary of missing values. Flags columns where
missingness exceeds a threshold.

## Usage

``` r
prep_missingness_report(data, threshold = 20, cols = NULL)
```

## Arguments

- data:

  Data frame.

- threshold:

  Numeric. Percentage threshold (0-100) above which a column is flagged
  as high-missingness. Default 20.

- cols:

  Character vector. Subset of columns to report on. NULL = all.

## Value

Data frame with columns: col_name, n_total, n_missing, pct_missing,
is_high_missing.
