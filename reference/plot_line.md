# Line Chart

Creates a line chart for trends over time or continuous variables.

## Usage

``` r
plot_line(
  data,
  x,
  y,
  group = NULL,
  color = NULL,
  facet = NULL,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  palette = "default",
  show_points = TRUE,
  show_labels = FALSE
)
```

## Arguments

- data:

  Data frame

- x:

  Character. Column name for x-axis (usually time/date)

- y:

  Character. Column name for y-axis

- group:

  Character. Optional column name for grouping lines

- color:

  Character. Optional column name for line color

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

- show_points:

  Logical. Show points on lines. Default TRUE.

- show_labels:

  Logical. Show value labels. Default FALSE.

## Value

A ggplot object
