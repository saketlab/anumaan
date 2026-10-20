# Compute YLL Attributable to AMR

Takes the per-patient YLL data produced by
[`daly_calc_yll_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_associated.md)
(which already contains `life_expectancy`, `death_weight`, and
`yll_contribution`) and multiplies by the mortality PAF from
[`daly_calc_paf_mortality()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_mortality.md)
to produce AMR-attributable YLL.

## Usage

``` r
daly_calc_yll_attributable(
  yll_patient_data,
  paf_mort,
  pathogen_col,
  patient_col,
  age_bin_col,
  sex_col = ".sex_norm",
  resistance_profile_col = NULL,
  profile_col = "profile",
  facility_col = NULL,
  syndrome_col = NULL,
  stratify_by = NULL
)
```

## Arguments

- yll_patient_data:

  Data frame: the `patient_data` element from
  [`daly_calc_yll_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_associated.md),
  with `life_expectancy`, `death_weight`, and `yll_contribution`
  columns.

- paf_mort:

  Named list returned by
  [`daly_calc_paf_mortality()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_mortality.md).

- pathogen_col:

  Character. Pathogen column.

- patient_col:

  Character. Patient identifier column.

- age_bin_col:

  Character. Age bin column.

- sex_col:

  Character. Normalised sex column. Default `".sex_norm"`.

- resistance_profile_col:

  Character or `NULL`. Column holding the resistance profile identifier
  per patient row. When supplied, per-profile PAF is applied instead of
  a scalar PAF_k. Default `NULL`.

- profile_col:

  Character. Profile identifier column in the `paf_mort` per-profile
  data frames. Default `"profile"`.

- facility_col:

  Character or `NULL`. Facility identifier column. When supplied results
  include a `by_facility` breakdown. Default `NULL`.

- syndrome_col:

  Character or `NULL`. Syndrome column.

- stratify_by:

  Character vector or `NULL`. Additional columns to aggregate results
  by. Default `NULL`.

## Value

A list:

- total:

  Scalar: total AMR-attributable YLL.

- per_pathogen:

  One row per pathogen: `pathogen_col`, `n_patients`,
  `YLL_associated_k`, `PAF_k_mort`, `YLL_attributable_k`.

- by_age_sex:

  YLL attributable stratified by age_bin x sex.

- by_pathogen_age_sex:

  YLL attributable by pathogen x age_bin x sex.

- by_facility:

  Per-facility breakdown (when `facility_col` supplied).

- by_syndrome:

  Per-syndrome breakdown (when `syndrome_col` supplied).

- by_syndrome_pathogen:

  Per-syndrome x pathogen (when `syndrome_col` supplied).

- stratified:

  User-defined stratification (when `stratify_by` supplied).

- patient_data:

  Patient-level data augmented with `PAF_kd` and
  `YLL_attributable_contribution`.

## Details

**Two PAF modes:**

- PAF_k scalar mode (default):

  When `resistance_profile_col` is `NULL`, the overall mortality PAF per
  pathogen k (`PAF_k_mort`) from `paf_mort` is used as a scalar
  multiplier.

- Per-profile mode:

  When `resistance_profile_col` is supplied, each patient row is matched
  to its resistance profile delta and the profile-specific
  `PAF_mortality` is applied. Patients whose profile does not appear in
  `paf_mort` receive `NA` and a warning is issued.

\$\$\text{YLL}^{\text{attr}}\_{i,k} = \text{yll\\contribution}\_{i,k}
\times \text{PAF}\_{k(,\delta)}\$\$

## Examples

``` r
cohort <- data.frame(
  patient_id    = sprintf("P%02d", 1:6),
  final_outcome = c("Death", "Death", "Discharged", "Death", "Discharged", "Death"),
  organism_name = rep(c("Escherichia coli", "Klebsiella pneumoniae"), 3),
  age_bin       = c("40-45", "60-65", "20-25", "0-1", "70-75", "85+"),
  sex           = c("Male", "Female", "Male", "Female", "Male", "Female")
)
yll <- daly_calc_yll_associated(cohort,
  outcome_col = "final_outcome", pathogen_col = "organism_name",
  patient_col = "patient_id", age_bin_col = "age_bin", sex_col = "sex"
)
#> Computing polymicrobial weights (global)...
#> Death cohort: 4 rows | 4 patients | 2 pathogens
#> YLL associated: total = 131.60 years | 2 pathogens | 4 patients
paf_mort <- list(
  "Escherichia coli"      = list(PAF_k_mort = 0.20),
  "Klebsiella pneumoniae" = list(PAF_k_mort = 0.35)
)
yll_attr <- daly_calc_yll_attributable(yll$patient_data, paf_mort,
  pathogen_col = "organism_name", patient_col = "patient_id",
  age_bin_col = "age_bin"
)
#> PAF mode: PAF_k scalar
#> YLL attributable: total = 41.09 years | associated = 131.60 years | 2 pathogens | 4 patients
yll_attr$per_pathogen
#> # A tibble: 2 × 5
#>   organism_name        n_patients YLL_associated_k YLL_attributable_k PAF_k_mort
#>   <chr>                     <int>            <dbl>              <dbl>      <dbl>
#> 1 Escherichia coli              1             33.1               6.62       0.2 
#> 2 Klebsiella pneumoni…          3             98.5              34.5        0.35
```
