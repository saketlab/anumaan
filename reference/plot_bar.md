# Generic Bar Plot

Creates a bar plot with flexible column mapping.

## Usage

``` r
plot_bar(
  data,
  x,
  y = NULL,
  fill = NULL,
  facet = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  palette = "default",
  show_labels = TRUE,
  flip_coords = FALSE,
  top_n = NULL,
  order = "desc",
  text_size = NULL,
  text_angle = NULL,
  bar_width = 0.8
)
```

## Arguments

- data:

  Data frame

- x:

  Character. Column name for x-axis (categorical variable)

- y:

  Character. Column name for y-axis (numeric variable). If NULL, uses
  count.

- fill:

  Character. Optional column name for fill color

- facet:

  Character. Optional column name for faceting

- title:

  Character. Plot title

- xlab:

  Character. X-axis label

- ylab:

  Character. Y-axis label

- palette:

  Character or named vector. Color palette name or custom colors

- show_labels:

  Logical. Show value labels on bars. Default TRUE.

- flip_coords:

  Logical. Flip coordinates (horizontal bars). Default FALSE.

- top_n:

  Numeric. Show only top N categories. Default NULL (show all).

- order:

  Character. Order bars by: "desc" (descending), "asc" (ascending), or
  "none". Default "desc".

- text_size:

  Numeric. Label text size. NULL for auto-adjustment. Default NULL.

- text_angle:

  Numeric. X-axis text angle (0-90). NULL for auto-adjustment. Default
  NULL.

- bar_width:

  Numeric. Bar width (0-1). Default 0.8.

## Value

A ggplot object
