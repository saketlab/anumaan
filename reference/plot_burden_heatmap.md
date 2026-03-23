# Burden Heatmap Across Hospitals

Heatmap showing all burden metrics (YLL, YLD, DALY x associated,
attributable) per 1000 cases across hospitals.

## Usage

``` r
plot_burden_heatmap(hospital_daly, label_map = NULL, syndrome_label = "cases")
```

## Arguments

- hospital_daly:

  Data frame from
  [`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md).

- label_map:

  Optional named character vector for center name labels.

- syndrome_label:

  Character. Label for infectious syndrome in legends. Default
  `"cases"`.

## Value

A `ggplot` object.
