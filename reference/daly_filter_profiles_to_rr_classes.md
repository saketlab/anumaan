# Filter Resistance Profiles to Classes with RR Estimates

Drops profiles whose resistant classes have no RR estimate in `rr_table`
and re-normalises the remaining probabilities. The all-susceptible
reference profile is always kept.

## Usage

``` r
daly_filter_profiles_to_rr_classes(
  profiles_with_rr,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  probability_col = "probability",
  fallback_rr = 1
)
```

## Arguments

- profiles_with_rr:

  Named list returned by
  [`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md).

- rr_table:

  Data frame of RR estimates.

- pathogen_col:

  Character. Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- probability_col:

  Character. Default `"probability"`.

- fallback_rr:

  Numeric. Default `1`.

## Value

Named list (one entry per pathogen) with filtered and re-normalised
profiles.
