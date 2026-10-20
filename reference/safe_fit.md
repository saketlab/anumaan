# Safely Fit a Distribution

Wrapper around
[`fitdistrplus::fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
that returns `NULL` instead of throwing an error when fitting fails.

## Usage

``` r
safe_fit(x, dist)
```

## Arguments

- x:

  Numeric vector of positive values (e.g., LOS in days).

- dist:

  Character. Distribution name: `"weibull"`, `"lnorm"`, or `"gamma"`.

## Value

A `fitdist` object, or `NULL` on failure.

## Examples

``` r
safe_fit(rlnorm(100), "lnorm")
#> Fitting of the distribution ' lnorm ' by maximum likelihood 
#> Parameters:
#>             estimate Std. Error
#> meanlog -0.004476225 0.09817305
#> sdlog    0.981730542 0.06941851
```
