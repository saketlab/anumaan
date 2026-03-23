# Stacked Bar Plot

Creates a stacked bar plot showing composition.

## Usage

``` r
plot_stacked_bar(
  data,
  x,
  fill,
  y = NULL,
  position = "stack",
  facet = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  palette = "default",
  show_labels = TRUE,
  show_percentages = TRUE,
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

  Character. Column name for stacking/fill variable

- y:

  Character. Optional column name for y-axis. If NULL, uses count.

- position:

  Character. "stack" (default), "fill" (100 percent stacked), or "dodge"
  (grouped)

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

- show_percentages:

  Logical. Show percentages in labels. Default TRUE.

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
