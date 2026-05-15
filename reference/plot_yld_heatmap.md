# Heatmap of YLD per 1 000 Admissions by Organism Group

Produces a two-column tile heatmap where the x-axis shows YLD Associated
and YLD Attributable, the y-axis shows organism groups (or individual
pathogens), and each cell's fill reflects the YLD burden per 1 000
admissions.

## Usage

``` r
plot_yld_heatmap(
  data,
  n_admissions = NULL,
  pathogen_col = "pathogen",
  assoc_col = "YLD_associated",
  attr_col = "YLD_attributable",
  group_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  show_values = TRUE,
  value_digits = 2,
  palette = "Blues",
  coord_fixed = TRUE,
  base_size = 14,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- data:

  Data frame. Must contain the columns named by `pathogen_col`,
  `assoc_col`, and `attr_col`.

- n_admissions:

  Numeric or `NULL`. Total admissions used as the denominator for per-1
  000 normalisation. If `NULL` (default), values are plotted as supplied
  (assumed pre-normalised).

- pathogen_col:

  Character. Column containing pathogen names. Default `"pathogen"`.

- assoc_col:

  Character. Column containing YLD associated values. Default
  `"YLD_associated"`.

- attr_col:

  Character. Column containing YLD attributable values. Default
  `"YLD_attributable"`.

- group_col:

  Character or `NULL`. Column containing organism group labels (e.g.
  `"org_group"`). If `NULL` (default) or not present in `data`, pathogen
  names are used on the y-axis.

- syndrome_col:

  Character or `NULL`. Syndrome filter column. Default `NULL`.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`. Default `NULL`.

- show_values:

  Logical. Print rounded values inside each cell. Default `TRUE`.

- value_digits:

  Integer. Decimal places for cell labels. Default 2.

- palette:

  Character. RColorBrewer palette for the fill scale. Default `"Blues"`.

- coord_fixed:

  Logical. Whether to apply
  [`ggplot2::coord_fixed()`](https://ggplot2.tidyverse.org/reference/coord_fixed.html)
  to produce square tiles. Default `TRUE`.

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom plot title. Auto-generated if `NULL`.

- subtitle:

  Character. Custom subtitle. Auto-generated if `NULL`.

## Value

A `ggplot` object.

## Details

The input `data` should contain one row per pathogen with
`YLD_associated` and `YLD_attributable` columns – typically the output
of
[`daly_calc_yld_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_attributable.md)
(which augments the data frame with both columns). If an organism-group
column has been pre-joined (e.g. `org_group`), supply its name via
`group_col` to aggregate cells to the group level.
