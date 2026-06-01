# Validate Pre-computed Aggregate Marginal Inputs

Validates a pre-computed marginal resistance table (e.g. from GBD or a
surveillance network) before it is passed to
[`estimate_profiles_convex()`](https://saketlab.github.io/anumaan/reference/estimate_profiles_convex.md).
Checks prevalence bounds, duplicate rows, zero-denominator cells,
impossible pairwise values, and minimum class counts per pathogen.

## Usage

``` r
validate_aggregate_inputs(
  marginals,
  pairwise = NULL,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rate_col = "marginal_resistance",
  n_tested_col = "n_tested",
  n_resistant_col = "n_resistant",
  class1_col = "antibiotic_class_1",
  class2_col = "antibiotic_class_2",
  pairwise_rate_col = "pairwise_resistance_prevalence",
  min_classes_per_pathogen = 2L
)
```

## Arguments

- marginals:

  Data frame with at minimum columns `pathogen_col`, `class_col`, and
  `rate_col`. Optionally also `n_tested_col` and `n_resistant_col`.

- pairwise:

  Data frame or `NULL`. Pairwise co-resistance table with columns
  `pathogen_col`, `class1_col`, `class2_col`, and `pairwise_rate_col`.
  Default `NULL`.

- pathogen_col:

  Character. Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- rate_col:

  Character. Marginal resistance rate column. Default
  `"marginal_resistance"`.

- n_tested_col:

  Character or `NULL`. Default `"n_tested"`.

- n_resistant_col:

  Character or `NULL`. Default `"n_resistant"`.

- class1_col:

  Character. Default `"antibiotic_class_1"`.

- class2_col:

  Character. Default `"antibiotic_class_2"`.

- pairwise_rate_col:

  Character. Default `"pairwise_resistance_prevalence"`.

- min_classes_per_pathogen:

  Integer. Pathogens with fewer classes are flagged. Default `2L`.

## Value

Tibble of check results (same schema as
[`validate_profile_inputs()`](https://saketlab.github.io/anumaan/reference/validate_profile_inputs.md)).
Stops on any `"fail"`.
