# Compute Marginal Resistance per Pathogen and Antibiotic Class

Collapses drug-level susceptibility data to antibiotic-class level per
isolate (\\R\_{e,k,c} = 1\\ if resistant to **any** drug in class
\\c\\), then computes marginal resistance for every pathogen x class
combination found in the data.

## Usage

``` r
compute_marginal_resistance(
  data,
  pathogen_col = "organism_name",
  org_group_col = "org_group",
  isolate_col = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_value_col = "antibiotic_value",
  zero_threshold = 0,
  min_n_tested = 30,
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL
)
```

## Arguments

- data:

  Data frame. Pre-processed AMR data at isolate x antibiotic level (one
  row per isolate-antibiotic combination). Results must already be
  binary (`"S"` / `"R"`); no reclassification is applied.

- pathogen_col:

  Character. Column with pathogen names. Default `"organism_name"`.

- org_group_col:

  Character. Column with organism group labels. Default `"org_group"`.

- isolate_col:

  Character. Column uniquely identifying each isolate. Default
  `"isolate_id"`.

- antibiotic_class_col:

  Character. Column with the antibiotic class for each drug. Default
  `"antibiotic_class"`.

- antibiotic_value_col:

  Character. Column with susceptibility result (`"S"` or `"R"`). Default
  `"antibiotic_value"`.

- zero_threshold:

  Numeric. Classes with `marginal_resistance <= zero_threshold` are
  listed in `$near_zero`. Default `0`.

- min_n_tested:

  Integer or `NULL`. Minimum number of isolates that must have been
  tested for a pathogen-class combination to be retained. Combinations
  with `n_tested < min_n_tested` are dropped from `$marginal` **and**
  from `$class_long`, so the exclusion propagates automatically into
  [`compute_pairwise_coresistance()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_coresistance.md)
  and
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).
  Set to `NULL` or `0` to disable the filter. Default `30`.

- facility_col:

  Character or `NULL`. Name of the column identifying the facility/site.
  When provided together with `facility_name`, data are filtered to the
  specified facility **before** any computation. Both `facility_col` and
  `facility_name` must be supplied together or both left `NULL`. When
  provided, `facility_col` is also retained in `$class_long` and
  `$marginal` so that downstream steps can apply the same filter.
  Default `NULL`.

- facility_name:

  Character or `NULL`. The facility value to retain (matched via `==`
  against `facility_col`). Default `NULL`.

- outcome_col:

  Character or `NULL`. Name of the column containing patient outcomes
  (e.g. `"final_outcome"`). When provided together with `outcome_value`,
  data are filtered to isolates with the specified outcome **before**
  any computation. Both `outcome_col` and `outcome_value` must be
  supplied together or both left `NULL`. When provided, `outcome_col` is
  retained in `$class_long` and `$marginal` so that downstream steps can
  apply the same filter. Default `NULL`.

- outcome_value:

  Character or `NULL`. The outcome value to retain (e.g. `"discharged"`,
  `"dead"`; matched via `==` against `outcome_col`). Default `NULL`.

## Value

Named list:

- `marginal`:

  Data frame with columns: `pathogen_col`, `org_group_col`,
  `antibiotic_class_col`, `n_tested`, `n_resistant`,
  `marginal_resistance`. Sorted descending by `marginal_resistance`
  within each pathogen.

- `near_zero`:

  Subset of `marginal` where `marginal_resistance <= zero_threshold`.
  These classes are candidates for exclusion in downstream profiling.

- `class_long`:

  Collapsed isolate x pathogen x class data frame (columns:
  `isolate_col`, `pathogen_col`, `org_group_col`,
  `antibiotic_class_col`, `class_result`). Pass this directly to
  [`compute_pairwise_coresistance()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_coresistance.md).

## Details

Classes whose marginal resistance is at or below `zero_threshold` are
listed in `$near_zero` as a flag for downstream use – they are **not**
removed here.

## Examples

``` r
if (FALSE) { # \dontrun{
marg <- compute_marginal_resistance(
  data                 = amr_clean,
  pathogen_col         = "organism_name",
  org_group_col        = "org_group",
  isolate_col          = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_value_col = "antibiotic_value"
)

marg$marginal # full marginal resistance table
marg$near_zero # classes flagged as near-zero
} # }
```
