# Filter Profiles to Classes with Actual RR Estimates

After
[`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md),
some profiles may be resistant only to antibiotic classes that have no
RR estimate in the 7-centres LOS data (e.g. a class present in the
susceptibility dataset but not tested in the LOS cohort). Those profiles
receive `fallback_rr = 1` and contribute zero to the PAF numerator while
their probability mass still enters the denominator, silently distorting
the estimate.

Drops profiles where every resistant class is absent from the RR table
(they would receive fallback_rr = 1 and silently distort the PAF), then
re-normalises the surviving profile probabilities to sum to 1.

## Usage

``` r
filter_profiles_to_rr_classes(
  profiles_with_rr,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  probability_col = "probability",
  fallback_rr = 1
)

filter_profiles_to_rr_classes(
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

  Named list from
  [`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md).

- rr_table:

  Data frame. Must have `pathogen_col` and `class_col`.

- pathogen_col:

  Character. Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- probability_col:

  Character. Default `"probability"`.

- fallback_rr:

  Numeric. Default `1`.

## Value

Named list with the same structure as `profiles_with_rr` but with
unmatched profiles removed and probabilities re-normalised.

Named list with unmatched profiles removed and probabilities
re-normalised.

## Details

This function drops any profile where **every** resistant class
(\\\delta_c = 1\\) is absent from the per-pathogen RR table, then
re-normalises the surviving profile probabilities to sum to 1.
