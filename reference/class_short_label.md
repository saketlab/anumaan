# Short antimicrobial class label (for dense facet strips)

Same lookup/fallback behaviour as
[`class_display_label()`](https://saketlab.github.io/anumaan/reference/class_display_label.md),
but returns the abbreviated form (e.g. `"3GC"`). Falls back to the first
4 characters of the underscore-cleaned name, upper-cased, for any class
not in the reference table.

## Usage

``` r
class_short_label(x)
```

## Arguments

- x:

  Character vector of raw class names.

## Value

Character vector of short labels, same length as `x`.
