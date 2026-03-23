# Compare Distribution Fits by AIC

Fits Weibull, Lognormal, and Gamma distributions to a numeric vector and
returns their AIC values for comparison.

## Usage

``` r
compare_distribution_aic(x, fits = NULL)
```

## Arguments

- x:

  Numeric vector of positive values (e.g., LOS in days).

- fits:

  Optional. Pre-computed fits from
  [`fit_distributions()`](https://saketlab.github.io/anumaan/reference/fit_distributions.md).

## Value

Data frame with columns `Weibull_AIC`, `Lognormal_AIC`, `Gamma_AIC`.

## Examples

``` r
compare_distribution_aic(rlnorm(200, meanlog = 2, sdlog = 0.5))
#>   Weibull_AIC Lognormal_AIC Gamma_AIC
#> 1    1150.919      1123.855   1129.85
```
