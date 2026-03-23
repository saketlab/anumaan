# Heatmap

Creates a heatmap (tile plot) for matrix data.

## Usage

``` r
plot_heatmap(
  data,
  x,
  y,
  fill,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  palette = "green",
  show_labels = TRUE
)
```

## Arguments

- data:

  Data frame

- x:

  Character. Column name for x-axis

- y:

  Character. Column name for y-axis

- fill:

  Character. Column name for fill value

- title:

  Character. Plot title

- xlab:

  Character. X-axis label

- ylab:

  Character. Y-axis label

- palette:

  Character or vector. Color palette. Default "green".

- show_labels:

  Logical. Show value labels in tiles. Default TRUE.

## Value

A ggplot object
