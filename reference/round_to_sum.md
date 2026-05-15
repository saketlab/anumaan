# Largest-Remainder Rounding

Rounds a numeric vector so that the individual rounded values sum
exactly to a target total. Uses the largest-remainder method (Hamilton
method). Input is coerced with
[`as.numeric()`](https://rdrr.io/r/base/numeric.html); non-numeric
values become `NA`.

## Usage

``` r
round_to_sum(x, target = round(sum(x, na.rm = TRUE)))
```

## Arguments

- x:

  Numeric vector to round.

- target:

  Integer target sum. Default is `round(sum(x, na.rm = TRUE))`.

## Value

Integer vector of the same length as `x` whose sum equals `target`
across non-missing entries; missing values are preserved as `NA`.

## Examples

``` r
round_to_sum(c(3.3, 3.3, 3.4), target = 10)
#> [1] 3 3 4
```
