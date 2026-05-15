# Get Age Bin Labels

Returns predefined age bin label sets for use with
[`prep_assign_age_bins()`](https://saketlab.github.io/anumaan/reference/prep_assign_age_bins.md).

## Usage

``` r
get_age_bins(type = "GBD_standard")
```

## Arguments

- type:

  Character. One of `"GBD_standard"` (5-year bins), `"pediatric"`,
  `"geriatric"`, or `"neonatal"`. The `"neonatal"` preset uses
  fractional-year labels corresponding to: `<0.02` (0-7 days),
  `0.02-0.08` (7-28 days), `0.08-0.25` (28-90 days), `0.25-1` (3
  months-1 year), then standard pediatric bands through adulthood.

## Value

Character vector of age bin labels.
