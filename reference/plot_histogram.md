# Histogram

Creates a histogram for continuous variables.

## Usage

``` r
plot_histogram(
  data,
  x,
  fill = NULL,
  facet = NULL,
  bins = 30,
  title = NULL,
  xlab = NULL,
  ylab = "Count",
  palette = "default"
)
```

## Arguments

- data:

  Data frame

- x:

  Character. Column name for continuous variable

- fill:

  Character. Optional column name for fill color

- facet:

  Character. Optional column name for faceting

- bins:

  Numeric. Number of bins. Default 30.

- title:

  Character. Plot title

- xlab:

  Character. X-axis label

- ylab:

  Character. Y-axis label

- palette:

  Character or named vector. Color palette

## Value

A ggplot object
