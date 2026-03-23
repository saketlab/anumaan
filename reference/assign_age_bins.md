# Assign Age Bins

Categorizes age into bins for stratification.

## Usage

``` r
assign_age_bins(data, age_col = "Age", bins = "GBD_standard")
```

## Arguments

- data:

  Data frame with Age column

- age_col:

  Character. Name of age column. Default "Age".

- bins:

  Character or numeric vector. Either "GBD_standard", "pediatric",
  "geriatric", or custom bin breaks. Default "GBD_standard".

## Value

Data frame with Age_bin column (factor)
