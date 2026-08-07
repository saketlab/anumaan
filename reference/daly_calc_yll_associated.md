# Compute YLL Associated with AMR (Patient-Level, Facility-Direct)

Computes years of life lost (YLL) associated with AMR directly from
patient-level facility records, without requiring population-level P_LK
(syndrome fractions) or R_kd (resistance profile scalars).

## Usage

``` r
daly_calc_yll_associated(
  data,
  outcome_col,
  death_value = "Death",
  pathogen_col,
  patient_col,
  age_bin_col,
  sex_col,
  facility_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  date_culture_col = NULL,
  specimen_col = NULL,
  poly_weight_method = "monomicrobial_proportion",
  min_mono_per_facility = 30L,
  gap_days = 14,
  le_path = system.file("extdata", "life_expectancy_all.xlsx", package = "anumaan"),
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c(`<1` = "0-1"),
  stratify_by = NULL
)
```

## Arguments

- data:

  Data frame of patient-level facility records.

- outcome_col:

  Character. Final outcome column.

- death_value:

  Character. Value(s) indicating a fatal outcome. Default `"Death"`.
  Pass a vector to match multiple labels (e.g. `c("Death","Died")`).

- pathogen_col:

  Character. Pathogen / organism column (k).

- patient_col:

  Character. Unique patient identifier column.

- age_bin_col:

  Character. Column containing GBD-standard age bin labels (e.g.
  `"0-1"`, `"1-5"`, ..., `"85+"`). Use `age_bin_map` to recode
  non-standard labels.

- sex_col:

  Character. Column containing patient sex.

- facility_col:

  Character or `NULL`. Facility identifier column. When supplied,
  polymicrobial weights are computed per facility and results include a
  `by_facility` breakdown. Default `NULL`.

- syndrome_col:

  Character or `NULL`. Syndrome column. When `NULL` all syndromes are
  pooled; when supplied the `by_syndrome` and `by_syndrome_pathogen`
  outputs are populated.

- syndrome_name:

  Character or `NULL`. If supplied, data are filtered to this syndrome
  before computation. `NULL` = all.

- date_culture_col:

  Character or `NULL`. Culture date column used for polymicrobial
  episode detection. When `NULL` polymicrobial flagging is skipped and
  all weights default to 1.

- specimen_col:

  Character or `NULL`. Specimen type column used alongside
  `date_culture_col` for polymicrobial detection.

- poly_weight_method:

  Character. Method for computing polymicrobial death weights. Default
  `"monomicrobial_proportion"`.

- min_mono_per_facility:

  Integer. Minimum monomicrobial records per facility to use
  facility-specific weights; smaller facilities fall back to global
  weights. Default `30L`.

- gap_days:

  Integer. Window (days) used for polymicrobial episode detection.
  Default `14`.

- le_path:

  Character. Path to the India life expectancy xlsx file. Defaults to
  the bundled `inst/extdata` copy.

- male_value:

  Character. Value in `sex_col` for males. Default `"Male"`.

- female_value:

  Character. Value in `sex_col` for females. Default `"Female"`. All
  other values use the combined LE.

- age_bin_map:

  Named character vector remapping non-standard age bin labels to
  LE-table labels. Default `c("<1" = "0-1")`.

- stratify_by:

  Character vector or `NULL`. Additional columns to aggregate results
  by. Default `NULL`.

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

## Details

For every fatal patient with pathogen k (and optionally syndrome L) the
individual YLL contribution is: \$\$\text{YLL}\_{r,k} =
\text{LE}(\text{age\\bin}\_r, \text{sex}\_r) \times w\_{r,k}\$\$ where
\\w\_{r,k}\\ is the polymicrobial death weight for patient r and
pathogen k (= 1 for monomicrobial, 0-1 for polymicrobial episodes).
Total YLL associated: \$\$\text{YLL}\_{\text{associated}} = \sum\_{r,k}
\text{YLL}\_{r,k}\$\$

**Polymicrobial weights** are computed via
[`prep_flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/prep_flag_polymicrobial.md) +
[`prep_compute_poly_weights()`](https://saketlab.github.io/anumaan/reference/prep_compute_poly_weights.md).
When `facility_col` is provided the weights are derived per facility
(reflecting local organism distributions), with automatic fallback to
globally-pooled proportions for facilities whose monomicrobial reference
pool is smaller than `min_mono_per_facility`. If `date_culture_col` is
`NULL` polymicrobial flagging is skipped and all weights default to 1.

## Examples

``` r
if (FALSE) { # \dontrun{
yll <- daly_calc_yll_associated(
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
