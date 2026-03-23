# Resistance Pattern Heatmap

Creates a heatmap showing resistance patterns across isolates and
antibiotic classes. Shows S (susceptible), R (resistant), Partial
(mixed), or NT (not tested) for each isolate-class combination.

## Usage

``` r
plot_resistance_heatmap(
  data,
  isolate_col,
  class_col,
  result_col,
  top_n = NULL,
  title = NULL,
  show_labels = FALSE
)
```

## Arguments

- data:

  Data frame with isolate/event data

- isolate_col:

  Character. Column name for isolate/event ID

- class_col:

  Character. Column name for antibiotic class

- result_col:

  Character. Column name for resistance result (S/R/I/NT)

- top_n:

  Numeric. Show only top N isolates. Default NULL (show all).

- title:

  Character. Plot title

- show_labels:

  Logical. Show result labels in tiles. Default FALSE.

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
plot_resistance_heatmap(data,
  isolate_col = "event_id",
  class_col = "antibiotic_class",
  result_col = "resistance_status",
  top_n = 20
)
} # }
```
