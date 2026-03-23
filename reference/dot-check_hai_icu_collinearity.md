# Check Collinearity Between HAI and ICU Covariates

Computes the phi (Pearson) correlation coefficient for the 2x2
contingency table of HAI x ICU. When ICU is hospital-acquired, the two
covariates can be highly correlated, making regression coefficients
unstable.

Computes the phi correlation for the HAI x ICU 2x2 table and warns when
it exceeds `phi_threshold`.

## Usage

``` r
.check_hai_icu_collinearity(
  df,
  hai_col = "HAI",
  icu_col = "ICU",
  phi_threshold = 0.7
)

.check_hai_icu_collinearity(
  df,
  hai_col = "HAI",
  icu_col = "ICU",
  phi_threshold = 0.7
)
```

## Arguments

- df:

  Patient-level data frame with integer 0/1 columns.

- hai_col:

  Character. Default `"HAI"`.

- icu_col:

  Character. Default `"ICU"`.

- phi_threshold:

  Numeric. Default `0.7`.

## Value

Named list: `phi`, `tbl` (2x2 table), `warning_issued` (logical).

Named list: `phi`, `tbl`, `warning_issued`.

## Details

Perfect separation (any 2x2 marginal = 0) is detected and reported
separately, as it guarantees model non-convergence.
