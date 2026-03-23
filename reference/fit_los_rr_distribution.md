# Estimate Per-Class LOS Relative Risk via Parametric Distribution Fitting

For each pathogen x antibiotic class, splits patients into resistant (R)
and susceptible (S) groups, fits a parametric distribution (default:
gamma) to each group's LOS, and returns: RR_LOS = mean_LOS(R) /
mean_LOS(S)

## Usage

``` r
fit_los_rr_distribution(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
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

  Character. Column used to derive HAI/CAI.

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

  Character. Date of culture column.

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

  Numeric. Hours after admission before a culture is classified as HAI.
  Default `48`.

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

  Fitted mean ratio: mean_LOS(R) / mean_LOS(S).

- mean_LOS_R:

  Fitted mean LOS for resistant patients.

- mean_LOS_S:

  Fitted mean LOS for susceptible patients.

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

This is an alternative to
[`fit_los_rr_poisson()`](https://saketlab.github.io/anumaan/reference/fit_los_rr_poisson.md)
that avoids regression and models the full LOS distribution per group.
The output flat data frame is fully compatible with
[`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
and
[`filter_profiles_to_rr_classes()`](https://saketlab.github.io/anumaan/reference/filter_profiles_to_rr_classes.md).

## See also

[`fit_los_rr_poisson`](https://saketlab.github.io/anumaan/reference/fit_los_rr_poisson.md),
[`assign_rr_to_profiles`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
