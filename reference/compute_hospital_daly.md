# Compute Hospital-Level DALY Breakdown

Distributes pooled YLL/YLD totals to individual hospitals proportionally
by deaths and discharges, and computes per-1000 rates.

## Usage

``` r
compute_hospital_daly(
  hospital_counts,
  total_deaths,
  total_discharged,
  yll_base,
  yll_associated,
  yll_attributable,
  yld_base,
  yld_associated,
  yld_attributable
)
```

## Arguments

- hospital_counts:

  Data frame with columns: `center_name`, `deaths_h`, `discharged_h`,
  `cases_h`.

- total_deaths:

  Numeric. Pooled total deaths across all hospitals.

- total_discharged:

  Numeric. Pooled total discharges.

- yll_base:

  Numeric. Total baseline YLL.

- yll_associated:

  Numeric. Total YLL associated with resistance.

- yll_attributable:

  Numeric. Total YLL attributable to resistance.

- yld_base:

  Numeric. Total baseline YLD.

- yld_associated:

  Numeric. Total YLD associated.

- yld_attributable:

  Numeric. Total YLD attributable.

## Value

Data frame with hospital-level YLL, YLD, DALY (absolute and per 1000
cases) columns.
