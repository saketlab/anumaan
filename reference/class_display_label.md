# Human-readable antimicrobial class label

Looks up `inst/extdata/class_display_labels.csv`; falls back to
replacing underscores with spaces for any class not in the table (e.g. a
pathogen-specific class not yet added there), so this never errors on an
unmapped class – it just degrades to a readable-if-generic label.

## Usage

``` r
class_display_label(x)
```

## Arguments

- x:

  Character vector of raw class names (e.g.
  `"Third_generation_cephalosporins"`).

## Value

Character vector of display labels, same length as `x`.
