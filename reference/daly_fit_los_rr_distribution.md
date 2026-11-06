# Estimate Per-Class LOS Relative Risk via Parametric Distribution Fitting

For each pathogen x antibiotic class, splits patients into resistant (R)
and susceptible (S) groups, fits the candidate parametric
distribution(s) (default: gamma) to each group's LOS, selects the best
fit by AIC, and takes the **theoretical median** of the selected
distribution: RR_LOS = median_LOS(R) / median_LOS(S)

## Usage

``` r
daly_fit_los_rr_distribution(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "infection_type",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  syndrome_name = NULL,
  organism_name = NULL,
  facility_name = NULL,
  hai_threshold_hours = 48,
  distributions = "gamma",
  max_los = 365,
  min_n = 10L
)
```

## Arguments

- data:

  Data frame of isolate/patient rows (long format).

- patient_id_col:

  Character. Patient identifier column.

- facility_col:

  Character. Facility/centre column.

- organism_col:

  Character. Organism name column.

- syndrome_col:

  Character. Syndrome column (used only when `syndrome_name` is not
  `NULL`).

- infection_type_col:

  Character. Column already holding the final HAI/CAI classification
  (e.g.
  [`prep_derive_hai_cai()`](https://saketlab.github.io/anumaan/reference/prep_derive_hai_cai.md)'s
  output). Used as-is – **not** re-derived from admission/culture dates.
  Default `"infection_type"`.

- antibiotic_class_col:

  Character. Antibiotic class column.

- antibiotic_name_col:

  Character. Antibiotic name column.

- antibiotic_value_col:

  Character. Antibiotic value column (S/I/R).

- date_admission_col:

  Character. Date of admission column.

- date_discharge_col:

  Character. Date of discharge column.

- date_culture_col:

  Character. Unused – retained only for backward API compatibility. LOS
  is computed as final outcome date - admission date and does not depend
  on the culture date.

- final_outcome_col:

  Character. Final outcome column.

- final_outcome_value:

  Character. Value indicating discharge.

- syndrome_name:

  Character or `NULL`. Restrict to this syndrome.

- organism_name:

  Character vector or `NULL`. Restrict to these pathogen(s); otherwise
  all pathogens in the filtered data.

- facility_name:

  Character or `NULL`. If provided, filters data to the specified
  facility before fitting. Default `NULL`.

- hai_threshold_hours:

  Numeric. Unused – retained only for backward API compatibility.
  HAI/CAI is read directly from `infection_type_col` rather than
  re-derived from a date gap.

- distributions:

  Character vector. Candidate distributions to fit. Default `"gamma"`.

- max_los:

  Numeric. Maximum plausible LOS in days. Default `200`.

- min_n:

  Integer. Minimum patients required in both R and S groups. Default
  `10L`.

## Value

A data frame with one row per pathogen x class:

- pathogen:

  Pathogen name.

- antibiotic_class:

  Antibiotic class name.

- RR_LOS:

  Ratio of fitted medians: median_LOS(R) / median_LOS(S).

- median_LOS_R:

  Fitted median LOS (days) for resistant patients.

- median_LOS_S:

  Fitted median LOS (days) for susceptible patients.

- median_LOS_R_years, median_LOS_S_years:

  The same medians in years (days / 365).

- best_dist_R:

  Distribution fitted to the R group.

- best_dist_S:

  Distribution fitted to the S group.

- n_R:

  Number of resistant patients.

- n_S:

  Number of susceptible patients.

- syndrome_scope:

  Syndrome filter applied, or "all".

## Details

Theoretical medians: Weibull (shape k, scale lambda) \\\lambda (\ln
2)^{1/k}\\; log-normal (meanlog mu) \\e^{\mu}\\; gamma (shape alpha,
rate beta) `qgamma(0.5, alpha, beta)` (no closed form). A class is
skipped when either group has fewer than `min_n` patients or its
distribution fit fails.

This is an alternative to
[`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md)
that avoids regression and models the full LOS distribution per group.
The output flat data frame is fully compatible with
`assign_rr_to_profiles()` and `filter_profiles_to_rr_classes()`.

## Two options for profile-level LOS

This function is **option 2** (per-class with max rule): its output is
passed to
[`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
(with `los_r_col = "median_LOS_R"`, `los_s_col = "median_LOS_S"`), which
gives each resistance profile the estimate of its dominant resistant
class. The S group is "susceptible to that class", not fully
susceptible. **Option 1** is
[`daly_fit_los_distn_by_profile()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_distn_by_profile.md),
which follows Section 5.3 directly by fitting LOS per patient-level
resistance profile against the fully susceptible profile; use it when
patient-level profiles are available and adequately populated.

## See also

[`daly_fit_los_rr`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md),
[`daly_assign_rr_to_profiles`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md),
[`daly_fit_los_distn_by_profile`](https://saketlab.github.io/anumaan/reference/daly_fit_los_distn_by_profile.md)
(option 1)
