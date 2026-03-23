# Grouped (Dodged) Bar Plot

Creates a grouped bar plot for side-by-side comparisons.

## Usage

``` r
plot_grouped_bar(
  data,
  x,
  fill,
  y = NULL,
  facet = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  palette = "default",
  show_labels = TRUE,
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

  Character. Column name for grouping variable

- y:

  Character. Optional column name for y-axis. If NULL, uses count.

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

- show_labels:

  Logical. Show value labels. Default TRUE.

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
plot_grouped_bar(data, x = "hospital", fill = "mrsa_status")
} # }
```
