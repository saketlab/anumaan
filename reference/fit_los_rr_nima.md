# Estimate Per-Profile LOS Relative Risk via Distribution Fitting (Nima Procedure)

For each pathogen k and each resistance profile delta from
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md),
collects the LOS vector of patients whose class-level resistance pattern
matches delta, fits Weibull, Lognormal, and Gamma distributions (best by
AIC), and computes: RR_LOS_k_delta = E\[LOS \| delta\] / E\[LOS \|
delta_0\] where delta_0 is the all-susceptible profile ("SSS...S").

## Usage

``` r
fit_los_rr_nima(
  data,
  resistance_profiles,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
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
  distributions = c("weibull", "lnorm", "gamma"),
  los_summary = c("mean", "median"),
  hai_threshold_hours = 48,
  max_los = 365,
  min_n = 10L
)
```

## Arguments

- data:

  Data frame. Merged AMR dataset at isolate x antibiotic level.

- resistance_profiles:

  Named list returned by
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).
  One entry per pathogen.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- facility_col:

  Character. Default `"center_name"`.

- facility_name:

  Character or `NULL`. When supplied, data are filtered to this facility
  before computation. Default `NULL`.

- organism_col:

  Character. Default `"organism_name"`.

- syndrome_col:

  Character. Default `"syndrome"`.

- infection_type_col:

  Character. Raw infection-type column used for HAI/CAI derivation.
  Default `"type_of_infection"`.

- antibiotic_class_col:

  Character. Default `"antibiotic_class"`.

- antibiotic_name_col:

  Character. Default `"antibiotic_name"`.

- antibiotic_value_col:

  Character. Default `"antibiotic_value"`.

- date_admission_col:

  Character. Default `"date_of_admission"`.

- date_discharge_col:

  Character. Default `"final_outcome_date"`.

- date_culture_col:

  Character. Default `"date_of_first_positive_culture"`.

- final_outcome_col:

  Character. Default `"final_outcome"`.

- final_outcome_value:

  Character. Default `"Discharged"`.

- syndrome_name:

  Character or NULL. Filter to syndrome. NULL = all.

- organism_name:

  Character or NULL. Pathogen(s) to process. NULL = all.

- distributions:

  Character vector. Distributions to try. Default
  `c("weibull","lnorm","gamma")`. Ignored when `los_summary = "median"`.

- los_summary:

  Character. How to summarise each profile's LOS vector before computing
  the RR. `"mean"` (default) fits the best parametric distribution
  (Weibull / Lognormal / Gamma by AIC) and returns the analytical mean.
  `"median"` uses the empirical median directly – no distribution
  fitting, more robust to outliers but ignores the shape of the LOS
  distribution.

- hai_threshold_hours:

  Numeric. Gap threshold (hours) for HAI/CAI derivation. Default `48`.

- max_los:

  Numeric. Cap on LOS (days). Default `200`.

- min_n:

  Integer. Minimum patients required per profile to fit a distribution.
  Profiles below this threshold return `NA` for RR. Default `10`.

## Value

Named list (one entry per pathogen k):

- `RR_k_delta`:

  Data frame with columns: `profile`, `probability`, `n_patients`,
  `best_dist`, `mean_LOS`, `RR_LOS` (relative to all-S profile).

- `classes`:

  Character vector of antibiotic class names (from resistance_profiles).

- `reference_profile`:

  The all-susceptible profile label (e.g. `"SSSS"`).

- `mean_LOS_S`:

  Fitted mean LOS for the all-S reference profile.

- `best_dist_S`:

  Best-fit distribution for the all-S profile.

- `syndrome_scope`:

  Syndrome filter applied, or `"all"`.

- `facility_scope`:

  Facility filter applied, or `"all"`.

## Details

LOS uses the HAI/CAI-specific clock (via `compute_patient_los`): HAI:
date_discharge - date_of_first_positive_culture CAI: date_discharge -
date_of_admission HAI/CAI is derived from `infection_type_col`; rows
where it is NA / "Not known" are classified by the culture-admission
gap.

Resistance is classified at antibiotic class level (class = R if any
drug in that class is R), matching the profile classes from
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).
Patients not tested for all profile classes are excluded from profile
assignment.

When `facility_name` is supplied the data are pre-filtered to that
facility; the `resistance_profiles` argument should correspondingly be
the facility-specific output of
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).
When `facility_name = NULL` all facilities are pooled.
