# Summarise a Fitted Distribution

Extracts mean, median, SD, and parameter values from a `fitdist` object
for Weibull, Lognormal, or Gamma distributions.

## Usage

``` r
summarise_distribution(fit, dist)
```

## Arguments

- fit:

  A `fitdist` object from
  [`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).

- dist:

  Character. One of `"weibull"`, `"lnorm"`, or `"gamma"`.

## Value

Data frame with columns `Mean_LOS`, `Median_LOS`, `SD_LOS`,
`Parameters`.

## Examples

``` r
fit <- fitdistrplus::fitdist(rlnorm(200), "lnorm")
summarise_distribution(fit, "lnorm")
#>         Mean_LOS Median_LOS SD_LOS                  Parameters
#> meanlog     1.55       0.98    1.9 meanlog=-0.023, sdlog=0.959
```
