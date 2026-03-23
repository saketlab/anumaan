# Fit Multiple Distributions

Fits Weibull, Lognormal, and Gamma distributions to a numeric vector,
returning all fit objects. Failed fits are `NULL`.

## Usage

``` r
fit_distributions(x)
```

## Arguments

- x:

  Numeric vector of positive values (e.g., LOS in days).

## Value

Named list with elements `weibull`, `lnorm`, `gamma`, each a `fitdist`
object or `NULL`.

## Examples

``` r
fits <- fit_distributions(rlnorm(200, 2, 0.5))
fits$lnorm$estimate
#>   meanlog     sdlog 
#> 2.0396379 0.5113832 
```
