# Plot YLL or YLD Associated vs Attributable by Organism

Produces a horizontal grouped bar chart showing the top `n` organisms
ranked by their associated burden (YLL or YLD), with a paired
attributable bar alongside each. Organisms are ranked by the associated
value so the highest-burden pathogen appears at the top.

## Usage

``` r
plot_burden_by_organism(
  data,
  metric = c("YLL", "YLD"),
  n = 10,
  organism_col = NULL,
  assoc_col = NULL,
  attr_col = NULL,
  n_admissions = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  label_digits = 1,
  colours = NULL,
  base_size = 14,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- data:

  Data frame. One row per organism. Must contain the columns named by
  `organism_col`, `assoc_col`, and `attr_col`.

- metric:

  Character. `"YLL"` or `"YLD"`. Controls default column names, colours,
  and axis labels. Ignored when `assoc_col` and `attr_col` are supplied
  explicitly. Default `"YLL"`.

- n:

  Integer. Number of top organisms to display. Default 10.

- organism_col:

  Character or `NULL`. Column containing organism / pathogen names. If
  `NULL` (default), auto-set to `"organism_name"` for YLL (matching
  `daly_calc_yll_*` examples) and `"pathogen"` for YLD (matching
  `daly_calc_yld_*` defaults).

- assoc_col:

  Character or `NULL`. Column for the associated burden value. If `NULL`
  (default), auto-set to `"YLL_associated_k"` (YLL) or
  `"YLD_associated"` (YLD) – matching the output column names of
  `compute_yll_associated()` and `compute_yld_associated()`
  respectively.

- attr_col:

  Character or `NULL`. Column for the attributable burden value. If
  `NULL` (default), auto-set to `"YLL_attributable_k"` (YLL) or
  `"YLD_attributable"` (YLD) – matching the output column names of
  `compute_yll_attributable()` and `compute_yld_attributable()`
  respectively.

- n_admissions:

  Numeric or `NULL`. When supplied, the associated and attributable
  values are divided by `n_admissions` and multiplied by 1 000 before
  plotting. Default `NULL` (no normalisation; values are plotted in
  thousands as returned by the YLL/YLD functions).

- syndrome_col:

  Character or `NULL`. Column name containing syndrome labels. If `NULL`
  (default), no filtering is applied.

- syndrome_name:

  Character or `NULL`. Syndrome value to retain. Requires
  `syndrome_col`. Default `NULL`.

- label_digits:

  Integer. Decimal places in bar value labels. Default 1.

- colours:

  Named character vector with elements `"associated"` and
  `"attributable"`. If `NULL` (default), auto-chosen based on `metric`.

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom plot title. Auto-generated if `NULL`.

- subtitle:

  Character. Custom subtitle. Auto-generated if `NULL`.

## Value

A `ggplot` object.

## Details

The input `data` should be a pathogen-level summary data frame with one
row per organism, already containing both the associated and
attributable burden columns. For YLD this is typically the output of
`compute_yld_associated()` left-joined with
`compute_yld_attributable()`. For YLL it is the `by_pathogen` slot of
`compute_yll_associated()` left-joined with the equivalent slot from
`compute_yll_attributable()`.
