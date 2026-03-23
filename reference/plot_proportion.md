# Proportion Bar Plot

Creates a 100

## Usage

``` r
plot_proportion(
  data,
  x,
  fill,
  facet = NULL,
  title = NULL,
  xlab = NULL,
  ylab = "Proportion",
  palette = "default",
  show_counts = TRUE,
  show_totals = TRUE,
  flip_coords = FALSE,
  order = "desc",
  text_size = NULL,
  text_angle = NULL,
  bar_width = 0.7
)
```

## Arguments

- data:

  Data frame

- x:

  Character. Column name for x-axis

- fill:

  Character. Column name for fill/stacking variable

- facet:

  Character. Optional column name for faceting

- title:

  Character. Plot title

- xlab:

  Character. X-axis label

- ylab:

  Character. Y-axis label

- palette:

  Character or named vector. Color palette

- show_counts:

  Logical. Show counts alongside percentages. Default TRUE.

- show_totals:

  Logical. Show total counts above bars. Default TRUE.

- flip_coords:

  Logical. Flip coordinates. Default FALSE.

- order:

  Character. Order bars by: "desc" (descending), "asc" (ascending), or
  "none". Default "desc".

- text_size:

  Numeric. Label text size. NULL for auto-adjustment. Default NULL.

- text_angle:

  Numeric. X-axis text angle (0-90). NULL for auto-adjustment. Default
  NULL.

- bar_width:

  Numeric. Bar width (0-1). Default 0.7.

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
plot_proportion(data,
  x = "organism", fill = "resistance_status",
  palette = "resistance"
)
} # }
```
