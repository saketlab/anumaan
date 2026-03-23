# DALY Lollipop Plot

Lollipop chart ranking hospitals by DALY (associated or attributable)
per 1000 cases, with point size proportional to case count.

## Usage

``` r
plot_daly_lollipop(
  hospital_daly,
  type = "attributable",
  label_map = NULL,
  syndrome_label = "cases"
)
```

## Arguments

- hospital_daly:

  Data frame from
  [`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md).

- type:

  Character. `"attributable"` or `"associated"`. Default
  `"attributable"`.

- label_map:

  Optional named character vector for center name labels.

- syndrome_label:

  Character. Label for infectious syndrome in axis/legend text. Default
  `"cases"`.

## Value

A `ggplot` object.
