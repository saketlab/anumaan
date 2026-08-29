# Human-readable class-pair label

Human-readable class-pair label

## Usage

``` r
class_pair_label(class_1, class_2, short = TRUE)
```

## Arguments

- class_1, class_2:

  Character vectors of raw class names.

- short:

  Logical; use
  [`class_short_label()`](https://saketlab.github.io/anumaan/reference/class_short_label.md)
  instead of
  [`class_display_label()`](https://saketlab.github.io/anumaan/reference/class_display_label.md).
  Default `TRUE` (pairwise facet strips get crowded fast with full
  names).

## Value

Character vector, e.g. `"AMG x CARB"`.
