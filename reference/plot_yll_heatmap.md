# Heatmap of YLL per 1 000 Admissions by Resistance Class and Pathogen Group

Produces a tile heatmap where the x-axis shows resistance classes
(dominant AMR profiles), the y-axis shows pathogen groups (or individual
pathogens), and each cell's fill reflects the associated or attributable
YLL burden per 1 000 admissions.

## Usage

``` r
plot_yll_heatmap(
  data,
  type = c("associated", "attributable"),
  value_col,
  n_admissions = NULL,
  pathogen_col = "pathogen",
  class_col = "profile",
  group_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  show_values = TRUE,
  value_digits = 2,
  palette = "Spectral",
  base_size = 14,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- data:

  Data frame. Must contain the columns named by `pathogen_col`,
  `class_col`, and `value_col`.

- type:

  Character. `"associated"` (default) or `"attributable"`. Used only for
  axis and title labels.

- value_col:

  Character. Required. Column containing the YLL value to plot (e.g.
  `"YLL_class"` after distributing per-pathogen YLL across resistance
  classes, or `"YLL_attributable_k"` if your data already has one row
  per pathogen x class).

- n_admissions:

  Numeric or `NULL`. Total admissions used as the denominator for per-1
  000 normalisation. If `NULL` (default), values are plotted as supplied
  (assumed pre-normalised).

- pathogen_col:

  Character. Column containing pathogen names. Default `"pathogen"`.

- class_col:

  Character. Column containing the resistance class / dominant profile.
  Default `"profile"`.

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

  Character. RColorBrewer palette for the fill scale. Default
  `"Spectral"`.

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom plot title. Auto-generated if `NULL`.

- subtitle:

  Character. Custom subtitle. Auto-generated if `NULL`.

## Value

A `ggplot` object.

## Details

The input `data` should be the `$by_class` slot (for associated) or the
`$by_profile` slot (for attributable) returned by
[`daly_calc_yll_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_associated.md)
/
[`daly_calc_yll_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_attributable.md).
Both slots contain columns `pathogen`, `profile`, and `YLL_Kdelta`.

If an organism-group column has been pre-joined onto `data` (e.g.
`org_group`), supply its name via `group_col` and the function will
aggregate to the group x class level. When `group_col` is `NULL` or
absent, individual pathogen names are used on the y-axis.
