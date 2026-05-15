# EDA ggplot2 Theme

Consistent theme used across all EDA plot functions in this file. Built
on top of
[`ggpubr::theme_pubr()`](https://rpkgs.datanovia.com/ggpubr/reference/theme_pubr.html).

## Usage

``` r
eda_theme(base_size = 14, legend_position = "top")
```

## Arguments

- base_size:

  Numeric. Base font size. Default 14.

- legend_position:

  Character. Legend position ("top", "bottom", "left", "right", "none").
  Default "top".

## Value

A ggplot2 theme object.
