# Compute Base YLL Block from a Scalar D_L (Steps 1-3)

Implements the first part of the YLL pipeline. Takes \\D_L\\ (deaths by
infectious syndrome, a single scalar), disaggregates it into age x sex
strata using observed proportions from the patient death cohort, and
multiplies by the India life expectancy for each stratum to yield the
base YLL block \\\sum_x D_L^x \\ e\_\*^x\\ (Equation 10).

## Usage

``` r
compute_base_yll_from_dl(
  dl,
  patient_data,
  patient_col,
  outcome_col,
  death_value = "Death",
  age_bin_col,
  sex_col,
  syndrome_col = NULL,
  syndrome_name = NULL,
  pathogen_col = NULL,
  facility_col = NULL,
  facility_name = NULL,
  use_sex = TRUE,
  le_path = system.file("extdata", "life_table_india.csv", package = "anumaan"),
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c(`<1` = "0-1"),
  stratify_by = NULL
)
```

## Arguments

- dl:

  Numeric scalar. \\D_L\\, the number of deaths for syndrome L (e.g.
  `45.2` deaths per 1000 incidences).

- patient_data:

  Data frame. Patient-level records (one row per patient x antibiotic
  test or per patient x pathogen). Used to derive observed age (x sex)
  proportions from the death cohort.

- patient_col:

  Character. Unique patient identifier column in `patient_data`.

- outcome_col:

  Character. Final-outcome column in `patient_data`.

- death_value:

  Character (scalar or vector). Value(s) in `outcome_col` indicating a
  fatal outcome. Default `"Death"`.

- age_bin_col:

  Character. Column in `patient_data` containing GBD-standard age bin
  labels (e.g. `"0-1"`, `"1-5"`, ..., `"85+"`). Use `age_bin_map` to
  recode non-standard labels (default remaps `"<1"` -\> `"0-1"`).

- sex_col:

  Character. Sex column in `patient_data`.

- syndrome_col:

  Character or `NULL`. Syndrome column in `patient_data` (e.g.
  `"infectious_syndrome"`). Required when `syndrome_name` is not `NULL`.

- syndrome_name:

  Character or `NULL`. If supplied, `patient_data` is filtered to rows
  where `syndrome_col == syndrome_name` before computing proportions and
  outputs (e.g. `"Bloodstream infection"`). `NULL` uses all syndromes.

- pathogen_col:

  Character or `NULL`. Pathogen (organism) column. When supplied,
  per-pathogen age proportions are computed and the `per_pathogen` and
  `by_pathogen_age_sex` outputs are populated.

- facility_col:

  Character or `NULL`. Facility identifier column. When supplied,
  per-facility proportions are computed and `by_facility` is populated.

- facility_name:

  Character or `NULL`. If provided, filters `patient_data` to the
  specified facility before computation.

- use_sex:

  Logical. Whether to disaggregate by sex as well as age (uses
  sex-specific life expectancy). Default `TRUE`.

- le_path:

  Character. Path to the India life expectancy xlsx file. Defaults to
  the bundled `inst/extdata` copy.

- male_value:

  Character. Value in `sex_col` for males. Default `"Male"`.

- female_value:

  Character. Value in `sex_col` for females. Default `"Female"`. All
  other values map to `"Combined"`.

- age_bin_map:

  Named character vector. Remaps non-standard age bin labels before
  joining to the life table. Default `c("<1" = "0-1")`.

- stratify_by:

  Character vector or `NULL`. Additional columns from `patient_data` to
  include in the `stratified` output.

## Value

A named list:

- `total`:

  Scalar: total base YLL = \\\sum_x D_L^x e\_\*^x\\.

- `by_age_sex`:

  Data frame: base YLL by age bin x sex.

- `per_pathogen`:

  Data frame: base YLL per pathogen K, computed using each pathogen's
  own death age distribution (only when `pathogen_col` is supplied).

- `by_pathogen_age_sex`:

  Data frame: base YLL by pathogen x age bin x sex (only when
  `pathogen_col` is supplied).

- `by_facility`:

  Data frame: base YLL by facility (only when `facility_col` is
  supplied).

- `by_syndrome_pathogen`:

  Data frame: base YLL by syndrome x pathogen (only when both
  `syndrome_col` and `pathogen_col` are supplied).

- `stratified`:

  Data frame: base YLL aggregated by `stratify_by` (only when
  `stratify_by` is supplied).

- `disaggregated_dl`:

  Data frame: the full expanded table with columns `D_x_L` (\\D_L^x\\),
  `proportion` (\\\hat{p}\_x\\), `life_expectancy` (\\e\_\*^x\\), and
  `yll_contribution` (\\D_L^x \times e\_\*^x\\) for every age x sex
  stratum. This is the row-level audit trail.

## Details

**Formula:** \$\$ \text{base\\YLL}\_L = \sum_x D_L^x \\ e\_\*^x \$\$
where \$\$D_L^x = D_L \times \hat{p}\_x\$\$ and \\\hat{p}\_x\\ is the
observed proportion of deaths in age bin \\x\\ (within sex \\s\\ when
`use_sex = TRUE`), estimated from `patient_data` (filtered to deaths,
and optionally to `syndrome_name`).

**D_L note:** For 1000-incidence normalisation: \$\$D_L = 1000 \times
\text{death\\rate}\_L,\quad \text{death\\rate}\_L =
\frac{\\\text{deaths}\mid L}{\\\text{incidence}\_L}\$\$ The caller
supplies this scalar directly via `dl`.

**Per-pathogen / per-facility YLL:** age (x sex) proportions are
re-computed within each subgroup so that each subgroup's YLL reflects
its own age structure, scaled by the shared `dl`.

## References

Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026.

Antimicrobial Resistance Collaborators. Global burden of bacterial
antimicrobial resistance in 2019. Lancet. 2022.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- compute_base_yll_from_dl(
  dl = 45.2,
  patient_data = bsi_data,
  outcome_col = "final_outcome",
  death_value = "Death",
  age_bin_col = "Age_bin",
  sex_col = "gender",
  syndrome_col = "infectious_syndrome",
  syndrome_name = "Bloodstream infection",
  pathogen_col = "organism_name",
  facility_col = "center_name",
  le_path = here::here(
    "anumaan", "inst", "extdata",
    "life_expectancy_all.xlsx"
  )
)

result$total
result$by_age_sex
result$per_pathogen
result$by_syndrome_pathogen
result$disaggregated_dl
} # }
```
