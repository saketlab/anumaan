# Compute YLL Attributable to AMR

Takes the per-patient YLL data produced by
[`compute_yll_associated()`](https://saketlab.github.io/anumaan/reference/compute_yll_associated.md)
(which already contains `life_expectancy`, `death_weight`, and
`yll_contribution`) and multiplies by the mortality PAF from
[`compute_paf_rr_mortality()`](https://saketlab.github.io/anumaan/reference/compute_paf_rr_mortality.md)
to produce AMR-attributable YLL.

Implements YLL – Attributable (Bhaswati Ganguli, DALY Methodology,
2026):

## Usage

``` r
compute_yll_attributable(
  profiles_with_rr,
  base_yll,
  P_Lk,
  pathogen_col = "pathogen",
  p_lk_col = "P_Lk",
  probability_col = "probability",
  rr_profile_col = "RR_death_profile",
  dominant_class_col = "dominant_class",
  patient_data = NULL,
  patient_col = NULL,
  outcome_col = NULL,
  death_value = "Death",
  age_bin_col = NULL,
  sex_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  patient_pathogen_col = NULL,
  le_path = NULL,
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c(`<1` = "0-1")
)

compute_yll_attributable(
  profiles_with_rr,
  base_yll,
  P_Lk,
  pathogen_col = "pathogen",
  p_lk_col = "P_Lk",
  probability_col = "probability",
  rr_profile_col = "RR_death_profile",
  dominant_class_col = "dominant_class",
  patient_data = NULL,
  patient_col = NULL,
  outcome_col = NULL,
  death_value = "Death",
  age_bin_col = NULL,
  sex_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  patient_pathogen_col = NULL,
  le_path = NULL,
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c(`<1` = "0-1")
)
```

## Arguments

- profiles_with_rr:

  Named list from `assign_rr_to_profiles(rr_col = "RR_death")`, one
  element per pathogen. Each element must have `probability_col`,
  `rr_profile_col`, and `dominant_class_col`.

- base_yll:

  Numeric scalar or list from
  [`compute_base_yll_from_dl()`](https://saketlab.github.io/anumaan/reference/compute_base_yll_from_dl.md)
  (uses `$total`).

- P_Lk:

  Data frame (pathogen column + P_Lk column) or list from
  [`calculate_P_Lk()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk.md).
  `NULL` applies equal weighting 1/K.

- pathogen_col:

  Character. Pathogen column in `P_Lk`. Default `"pathogen"`.

- p_lk_col:

  Character. P_Lk value column. Default `"P_Lk"`.

- probability_col:

  Character. Profile prevalence R'\_Kdelta column. Default
  `"probability"`.

- rr_profile_col:

  Character. Dominant-class mortality RR column. Default
  `"RR_death_profile"`.

- dominant_class_col:

  Character. Dominant class column. Default `"dominant_class"`.

- patient_data:

  Data frame or `NULL`. Patient-level records for per-patient
  attributable YLL output.

- patient_col:

  Character or `NULL`. Patient identifier column in `patient_data`.

- outcome_col:

  Character or `NULL`. Final outcome column in `patient_data`.

- death_value:

  Character. Value(s) indicating a fatal outcome. Default `"Death"`.

- age_bin_col:

  Character or `NULL`. Age bin column in `patient_data`.

- sex_col:

  Character or `NULL`. Sex column in `patient_data`.

- syndrome_col:

  Character or `NULL`. Syndrome column.

- syndrome_name:

  Character or `NULL`. If supplied, filters `patient_data` to this
  syndrome.

- patient_pathogen_col:

  Character or `NULL`. Pathogen column name in `patient_data` (e.g.
  `"organism_name"`).

- le_path:

  Character or `NULL`. Path to the life expectancy file.

- male_value:

  Character. Value in `sex_col` for males. Default `"Male"`.

- female_value:

  Character. Value in `sex_col` for females. Default `"Female"`.

- age_bin_map:

  Named character vector remapping non-standard age bin labels. Default
  `c("<1" = "0-1")`.

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

Named list:

- `by_profile`:

  Data frame: one row per pathogen x resistant profile with `R_prime`,
  `RR_dominant`, `excess_risk`, `numerator_delta`, `denominator`,
  `MortPAF`, `P_Lk`, `base_yll`, `YLL_Kdelta`.

- `by_pathogen`:

  Data frame: one row per pathogen with `P_Lk`, `denominator`,
  `sum_MortPAF`, `YLL_K`.

- `total_yll`:

  Scalar: \\\sum_K YLL_K\\.

- `base_yll_used`:

  Scalar base YLL used.

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

**Eq. 15 – Mortality PAF per profile delta:** \$\$
\text{MortPAF}\_{K\delta} = \frac{R'\_{K\delta} \cdot
(RR\_{Kd^\*}(\delta) - 1)} {1 + \sum\_{\delta'} R'\_{K\delta'} \cdot
(RR\_{Kd^\*}(\delta') - 1)} \$\$

Numerator is **per profile** delta; denominator is **shared** across all
resistant profiles of pathogen K. The sum \\\sum\_{\delta'}\\ runs over
resistant profiles only (`dominant_class != "all_susceptible"`).

**Eq. 14 – Attributable YLL per profile:** \$\$ YLL\_{K\delta} =
\text{BaseYLL} \times P\_{KL} \times \text{MortPAF}\_{K\delta} \$\$

**Attributable YLL per pathogen:** \$\$YLL_K = \sum\_\delta
YLL\_{K\delta}\$\$

## References

Bhaswati Ganguli. DALY Methodology for AMR. March 2026. Eq. 14-15.

## Examples

``` r
if (FALSE) { # \dontrun{
yll_assoc <- compute_yll_associated(data = cohort, ...)
mort_or <- fit_mortality_rr_logistic(data = rr_data, ...)
profiles_or <- assign_rr_to_profiles(profiles_out,
  rr_table = mort_or,
  rr_col = "OR_death"
)
paf_mort <- compute_paf_rr_mortality(profiles_or)

yll_attr <- compute_yll_attributable(
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
