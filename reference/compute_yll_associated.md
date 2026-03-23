# Compute YLL Associated with AMR (Patient-Level, Facility-Direct)

Computes years of life lost (YLL) associated with AMR directly from
patient-level facility records, without requiring population-level P_LK
(syndrome fractions) or R_kd (resistance profile scalars).

Implements YLL – Associated (Eq. 12):

## Usage

``` r
compute_yll_associated(
  base_yll,
  P_Lk,
  fatal_resistance_prevalence,
  pathogen_col = "pathogen",
  p_lk_col = "P_Lk",
  profiles_with_rr = NULL,
  rr_profile_col = "RR_death_profile",
  probability_col = "probability",
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

compute_yll_associated(
  base_yll,
  P_Lk,
  fatal_resistance_prevalence,
  pathogen_col = "pathogen",
  p_lk_col = "P_Lk",
  profiles_with_rr = NULL,
  rr_profile_col = "RR_death_profile",
  probability_col = "probability",
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

- base_yll:

  Numeric scalar or list from
  [`compute_base_yll_from_dl()`](https://saketlab.github.io/anumaan/reference/compute_base_yll_from_dl.md)
  (uses `$total`).

- P_Lk:

  Named list from
  [`calculate_P_Lk()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk.md)
  or the P_Lk data frame directly (column `pathogen_col` and
  `p_lk_col`). `NULL` applies equal weight 1/K.

- fatal_resistance_prevalence:

  Named list from
  [`compute_fatal_resistance_prevalence()`](https://saketlab.github.io/anumaan/reference/compute_fatal_resistance_prevalence.md).

- pathogen_col:

  Character. Pathogen column in P_Lk. Default `"pathogen"`.

- p_lk_col:

  Character. P_Lk value column. Default `"P_Lk"`.

- profiles_with_rr:

  Named list or `NULL`. Output from
  [`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
  for per-profile / per-class breakdowns.

- rr_profile_col:

  Character. Dominant-class RR column in profile data. Default
  `"RR_death_profile"`.

- probability_col:

  Character. Profile prevalence column. Default `"probability"`.

- dominant_class_col:

  Character. Column identifying the dominant antibiotic class. Default
  `"dominant_class"`.

- patient_data:

  Data frame or `NULL`. Patient-level records for per-patient YLL
  output.

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

A named list:

- `total`:

  Scalar: total YLL associated across all pathogens and facilities.

- `per_pathogen`:

  Data frame: YLL summed per pathogen k, pooled across facilities.
  Columns: `pathogen_col`, `n_patients`, `YLL_associated_k`.

- `by_age_sex`:

  Data frame: YLL by `age_bin_col` x sex.

- `by_pathogen_age_sex`:

  Data frame: YLL by pathogen x `age_bin_col` x sex.

- `by_facility`:

  Data frame (only when `facility_col` is supplied): per-facility YLL,
  one row per facility x pathogen.

- `by_syndrome`:

  Data frame (only when `syndrome_col` is supplied and `syndrome_name`
  is `NULL`): YLL by syndrome.

- `by_syndrome_pathogen`:

  Data frame (only when `syndrome_col` supplied): YLL by syndrome x
  pathogen.

- `stratified`:

  Data frame (only when `stratify_by` is supplied): YLL aggregated by
  the requested columns.

- `patient_data`:

  The death-cohort data frame used for computation, with
  `polymicrobial_weight`, `life_expectancy`, and `yll_contribution`
  columns attached.

Named list:

- `by_pathogen`:

  Data frame: `pathogen`, `P_Lk`, `sum_r_prime`, `susceptible_fraction`,
  `denominator`, `R_K_star`, `base_yll`, `YLL_K`.

- `by_pathogen_profile`:

  Data frame: one row per pathogen x resistant profile with `P_Lk`,
  `R_star`, `YLL_Kdelta`.

- `total_yll`:

  Scalar: sum of YLL_K.

- `base_yll_used`:

  Scalar base YLL.

## Details

For every fatal patient with pathogen k (and optionally syndrome L) the
individual YLL contribution is: \$\$\text{YLL}\_{r,k} =
\text{LE}(\text{age\\bin}\_r, \text{sex}\_r) \times w\_{r,k}\$\$ where
\\w\_{r,k}\\ is the polymicrobial death weight for patient r and
pathogen k (= 1 for monomicrobial, 0-1 for polymicrobial episodes).
Total YLL associated: \$\$\text{YLL}\_{\text{associated}} = \sum\_{r,k}
\text{YLL}\_{r,k}\$\$

**Polymicrobial weights** are computed via
[`flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/flag_polymicrobial.md) +
[`compute_polymicrobial_weight()`](https://saketlab.github.io/anumaan/reference/compute_polymicrobial_weight.md)
from `weight.R`. When `facility_col` is provided the weights are derived
per facility (reflecting local organism distributions), with automatic
fallback to globally-pooled proportions for facilities whose
monomicrobial reference pool is smaller than `min_mono_per_facility`. If
`date_culture_col` is `NULL` polymicrobial flagging is skipped and all
weights default to 1.

\$\$ YLL_K = \left\[\sum_L \left(\sum_x D_L^x e_x^\*\right) \times
P\_{KL}\right\] \times R_K^\* \$\$

where \\R_K^\* = \sum\_\delta R^\*\_{K\delta}\\ (sum of per-profile
fatal resistance prevalences from
`compute_fatal_resistance_prevalence`). For a single syndrome (BSI):
YLL_K = base_YLL x P_KL x R_K\*.

Also returns per-profile YLL_Kdelta = base_YLL x P_KL x R\*\_Kdelta.

## References

Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026. Eq.
12.

## Examples

``` r
if (FALSE) { # \dontrun{
yll <- compute_yll_associated(
  data             = cohort_df,
  outcome_col      = "final_outcome",
  death_value      = "Death",
  pathogen_col     = "organism_name",
  patient_col      = "PatientInformation_id",
  age_bin_col      = "Age_bin",
  sex_col          = "gender",
  facility_col     = "center_name",
  syndrome_col     = "infectious_syndrome",
  date_culture_col = "culture_date",
  specimen_col     = "sample_type",
  stratify_by      = c("location", "infectious_syndrome")
)
yll$total
yll$per_pathogen
yll$by_age_sex
yll$stratified
} # }
```
