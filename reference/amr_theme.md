# AMR Theme for ggplot2

Universal theme for all AMR plots with consistent styling.

## Usage

``` r
amr_theme(base_size = 16, legend_position = "top")
```

## Arguments

- base_size:

  Base font size. Default 12.

- legend_position:

  Legend position. Default "top".

## Value

A ggplot2 theme object

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  amr_theme()
```
