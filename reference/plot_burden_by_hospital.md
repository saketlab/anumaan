# Plot YLL or YLD Associated vs Attributable per Hospital

Produces a horizontal grouped bar chart comparing the associated and
attributable burden (YLL or YLD) per 1 000 patients across hospitals.
Bars are sorted by the associated value (largest at top by default).
Each bar is labelled with its rounded value.

## Usage

``` r
plot_burden_by_hospital(
  data,
  metric = c("YLL", "YLD"),
  facility_col = "center_name",
  center_col = NULL,
  assoc_col = NULL,
  attr_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  sort_by = c("associated", "attributable", "none"),
  label_digits = 1,
  colours = NULL,
  base_size = 14,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- data:

  Data frame. One row per hospital (or one row per hospital x syndrome).
  Must contain the columns named by `facility_col`, `assoc_col`, and
  `attr_col`.

- metric:

  Character. `"YLL"` or `"YLD"`. Controls which default column names and
  axis labels are used. Ignored when `assoc_col` and `attr_col` are
  supplied explicitly. Default `"YLL"`.

- facility_col:

  Character. Column containing hospital / facility names. Default
  `"center_name"`.

- center_col:

  Deprecated alias for `facility_col`.

- assoc_col:

  Character or `NULL`. Column for the associated per-1 000 value. If
  `NULL` (default), auto-set to `"YLL_assoc_per_1000"` or
  `"YLD_assoc_per_1000"` based on `metric`.

- attr_col:

  Character or `NULL`. Column for the attributable per-1 000 value. If
  `NULL` (default), auto-set to `"YLL_attr_per_1000"` or
  `"YLD_attr_per_1000"` based on `metric`.

- syndrome_col:

  Character or `NULL`. Column name containing syndrome labels (e.g.
  `"infectious_syndrome"`). If `NULL` (default), no syndrome filtering
  is applied.

- syndrome_name:

  Character or `NULL`. The syndrome value to retain (e.g.
  `"Bloodstream infection"`). Requires `syndrome_col` to be set. If
  `NULL` (default), all rows are used.

- sort_by:

  Character. `"associated"` (default), `"attributable"`, or `"none"`.
  Controls hospital ordering on the y-axis.

- label_digits:

  Integer. Decimal places in bar value labels. Default 1.

- colours:

  Named character vector with exactly two elements: `"associated"` and
  `"attributable"`. If `NULL` (default), a colour pair is chosen
  automatically based on `metric` (red tones for YLL, blue tones for
  YLD).

- base_size:

  Numeric. Base font size. Default 14.

- title:

  Character. Custom plot title. Auto-generated if `NULL`.

- subtitle:

  Character. Custom subtitle. Auto-generated if `NULL`.

## Value

A `ggplot` object.

## Details

The input `data` should be a hospital-level summary data frame (e.g. the
output of
[`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md))
that already contains the per-1 000 columns. If the summary was computed
across multiple syndromes, use `syndrome_col` and `syndrome_name` to
restrict the plot to one syndrome before pivoting.
