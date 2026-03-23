# Burden Comparison Bar Plot

Horizontal grouped bar chart comparing associated vs attributable burden
(YLL, YLD, or DALY) across hospitals.

## Usage

``` r
plot_burden_comparison(
  hospital_daly,
  metric = "DALY",
  label_map = NULL,
  syndrome_label = "cases"
)
```

## Arguments

- hospital_daly:

  Data frame from
  [`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md).

- metric:

  Character. One of `"YLL"`, `"YLD"`, or `"DALY"`.

- label_map:

  Optional named character vector mapping center names to display
  labels. Default `NULL` (use center names as-is).

- syndrome_label:

  Character. Label for the infectious syndrome used in axis/legend text
  (e.g., `"BSI patients"`, `"UTI cases"`). Default `"cases"`.

## Value

A `ggplot` object.
