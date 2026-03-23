# Largest-Remainder Rounding

Rounds a numeric vector so that the individual rounded values sum
exactly to a target total. Uses the largest-remainder method (Hamilton
method).

## Usage

``` r
round_to_sum(x, target = round(sum(x)))
```

## Arguments

- x:

  Numeric vector to round.

- target:

  Integer target sum. Default is `round(sum(x))`.

## Value

Integer vector of the same length as `x` whose sum equals `target`.

## Examples

``` r
round_to_sum(c(3.3, 3.3, 3.4), target = 10)
#> [1] 3 3 4
```
