# Plot LOS Distribution with Fitted Overlays

Creates a histogram of LOS values with Weibull, Lognormal, and Gamma
density curves overlaid.

## Usage

``` r
plot_los_distributions(los_vec, title, bins = 35, fits = NULL)
```

## Arguments

- los_vec:

  Numeric vector of LOS values (days).

- title:

  Character. Plot title.

- bins:

  Integer. Number of histogram bins. Default 35.

- fits:

  Optional. Pre-computed fits from
  [`fit_distributions()`](https://saketlab.github.io/anumaan/reference/fit_distributions.md).

## Value

A `ggplot` object.

## Examples

``` r
plot_los_distributions(rlnorm(200, 2, 0.5), "Example LOS Distribution")
```
