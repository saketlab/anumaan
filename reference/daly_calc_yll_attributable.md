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
if (FALSE) { # \dontrun{
yll_assoc <- daly_calc_yll_associated(data = cohort, ...)
mort_or <- fit_mortality_rr_logistic(data = rr_data, ...)
profiles_or <- assign_rr_to_profiles(profiles_out,
  rr_table = mort_or,
  rr_col = "OR_death"
)
paf_mort <- daly_calc_paf_mortality(profiles_or)

yll_attr <- daly_calc_yll_attributable(
  yll_patient_data = yll_assoc$patient_data,
  paf_mort         = paf_mort,
  pathogen_col     = "organism_name",
  patient_col      = "PatientInformation_id",
  age_bin_col      = "Age_bin",
  sex_col          = ".sex_norm",
  facility_col     = "center_name",
  syndrome_col     = "infectious_syndrome",
  stratify_by      = c("location", "infectious_syndrome")
)
} # }
```
