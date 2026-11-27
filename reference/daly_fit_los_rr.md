# Fit Relative Length of Stay with a Gamma GLM (Log Link)

For each pathogen and antibiotic class, fits a Gamma generalised linear
model with log link to the length of stay (LOS) of discharged
admissions:

## Usage

``` r
daly_fit_los_rr(
  data,
  admission_id_col = "PatientInformation_id",
  cluster_id_col = NULL,
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
  age_col = "Age",
  sex_col = "Gender",
  unit_type_col = NULL,
  comorbidity_col = NULL,
  syndrome_name = NULL,
  organism_name = NULL,
  classes = NULL,
  hospital_effect = c("fixed", "none", "stratified"),
  reference_hospital = NULL,
  intermediate_as_resistant = TRUE,
  hai_value = "HAI",
  cai_value = "CAI",
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  max_los = 365,
  min_n = 20L,
  min_resistant = 5L,
  min_susceptible = 5L,
  flag_min = 10L,
  ci_method = c("bootstrap", "delta", "none"),
  n_boot = 1000L,
  conf_level = 0.95,
  seed = NULL,
  return_draws = FALSE,
  patient_id_col = NULL,
  hai_threshold_hours = NULL
)
```

## Arguments

- data:

  Data frame with one row per antibiotic test result.

- admission_id_col:

  Character. Admission identifier; the unit of analysis. Default
  `"PatientInformation_id"`.

- cluster_id_col:

  Character or `NULL`. Identifier used to cluster repeated admissions of
  the same patient in the bootstrap and the sandwich covariance. `NULL`
  clusters by admission.

- facility_col:

  Character or `NULL`. Hospital column. `NULL` pools everything and
  forces `hospital_effect = "none"`.

- organism_col:

  Character. Pathogen column.

- syndrome_col:

  Character. Syndrome column (used with `syndrome_name`).

- infection_type_col:

  Character. Column holding the final HAI/CAI classification. Default
  `"infection_type"`.

- antibiotic_class_col, antibiotic_name_col, antibiotic_value_col:

  Character. Antibiotic class, drug name and S/I/R result columns.

- date_admission_col, date_discharge_col:

  Character. Admission and final outcome date columns.

- date_culture_col:

  Character. Culture date column, used to take covariates from each
  admission's earliest culture of the pathogen. When absent, the first
  row is used.

- final_outcome_col, final_outcome_value:

  Character. Final outcome column and the value marking discharged
  admissions.

- age_col, sex_col:

  Character or `NULL`. Age and sex columns; `NULL` leaves the covariate
  out.

- unit_type_col:

  Character or `NULL`. Ward / ICU column.

- comorbidity_col:

  Character or `NULL`. Comorbidity column.

- syndrome_name:

  Character or `NULL`. Restrict to one syndrome.

- organism_name:

  Character vector or `NULL`. Pathogens to model; `NULL` models all.

- classes:

  `NULL` (all classes), a character vector (the same classes for every
  pathogen), or a named list of character vectors (a panel per pathogen,
  names matching `organism_name` values).

- hospital_effect:

  One of `"fixed"`, `"none"`, `"stratified"`. See "Hospital effect".

- reference_hospital:

  Character or `NULL`. Reference level of the hospital fixed effect;
  `NULL` uses the hospital with most admissions in each model. Does not
  change any estimate.

- intermediate_as_resistant:

  Logical. See "Resistance".

- hai_value, cai_value:

  Character. Values of `infection_type_col` for HAI and CAI; other
  values are treated as missing.

- icu_values:

  Character vector. Values of `unit_type_col` counted as ICU (matched
  case-insensitively as substrings).

- max_los:

  Numeric. Maximum LOS in days. Default 365.

- min_n, min_resistant, min_susceptible:

  Integer. Minimum admissions, resistant admissions and susceptible
  admissions for a model to be fitted.

- flag_min:

  Integer. Hospitals with fewer resistant or susceptible admissions than
  this are flagged `low_support` in the hospital table (they are still
  estimated).

- ci_method:

  One of `"bootstrap"`, `"delta"`, `"none"`.

- n_boot:

  Integer. Bootstrap replicates. Default 1000.

- conf_level:

  Numeric. Confidence level. Default 0.95.

- seed:

  Integer or `NULL`. Random seed for the bootstrap.

- return_draws:

  Logical. Return the bootstrap replicates.

- patient_id_col:

  Deprecated. Use `admission_id_col`.

- hai_threshold_hours:

  Deprecated and ignored. HAI/CAI is read from `infection_type_col`.

## Value

A list with:

- `rr`:

  One row per pathogen x class (x hospital for `"stratified"`): `RR_LOS`
  with `CI_lower`, `CI_upper`, `beta_log_rr`, `se_log_rr_cluster`;
  counts; observed mean LOS of R and S admissions; all-hospital `LOSR`,
  `LOSS` and `excess_LOS` with intervals; `dispersion`; `RxHospital_p`;
  covariates used and settings.

- `hospital`:

  One row per pathogen x class x hospital: counts, observed and
  predicted mean LOS of R and S admissions (a calibration check),
  `LOSR`, `LOSS`, `excess_LOS` with intervals, `RR_LOS`, `low_support`
  and `fallback_all_admissions`.

- `draws`:

  When `return_draws = TRUE` and `ci_method = "bootstrap"`: a list named
  `"pathogen | class"` (`"pathogen | class | hospital"` for
  `"stratified"`), each with `RR_LOS` (vector), `LOSR`, `LOSS`
  (replicate x hospital matrices, days) and `LOSR_all`, `LOSS_all`
  (all-hospital values, vectors). Pass it to
  [`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md)
  /
  [`daly_calc_yld_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_attributable.md)
  as `los_draws`.

- `skipped`:

  Pathogen x class (x hospital) combinations not fitted, with the
  reason.

## Details

\$\$\log E(LOS_i) = \beta_0 + \beta_R R_i + \beta_1 HAI_i + \beta_2
Age_i + \beta_3 Sex_i + \beta_4 ICU_i + \beta_5 Comorbidity_i +
\gamma\_{h(i)}\$\$

and returns the relative LOS, \\RR\_{LOS} = \exp(\beta_R)\\, together
with the adjusted LOS of resistant admissions (`LOSR`) and of the same
admissions had they been susceptible (`LOSS`), per hospital and for all
hospitals combined.

## Unit of analysis

One row per admission x pathogen x antibiotic class. LOS is computed
once per admission by
[`daly_compute_patient_los()`](https://saketlab.github.io/anumaan/reference/daly_compute_patient_los.md)
(final outcome date minus admission date, discharged admissions only, 0
\< LOS \<= `max_los`). An admission with several pathogens enters each
pathogen's models.

## Resistance (worst phenotype)

For an admission and pathogen, a class is resistant (`R = 1`) when any
result of that pathogen for any drug in the class is R, and susceptible
(`R = 0`) when it was tested and no result is R. Results of other
pathogens in the same admission are not used. With
`intermediate_as_resistant = TRUE`, I counts as R except for colistin,
where I counts as S (the convention of
[`prep_recode_intermediate_ast()`](https://saketlab.github.io/anumaan/reference/prep_recode_intermediate_ast.md));
with `FALSE`, I results are dropped. Admissions not tested for a class
are excluded from that class's model.

## Covariates

Taken from each admission's earliest culture of the pathogen
(`date_culture_col`): HAI (1 if `infection_type_col` equals `hai_value`,
0 if it equals `cai_value`, otherwise missing), age (linear), sex, ICU
(when `unit_type_col` is given) and comorbidity (when `comorbidity_col`
is given). Admissions with a missing value in a requested covariate are
dropped from that model and counted in `n_dropped_covariates`. A
covariate with no variation in a model is left out of that model.

## Hospital effect

- `"fixed"` (default):

  One model per pathogen x class across all hospitals, with hospital as
  a fixed effect. Each hospital has its own LOS level; the resistance
  effect (`RR_LOS`) is shared. Small hospitals are kept and borrow the
  resistance effect.

- `"none"`:

  One model with no hospital term.

- `"stratified"`:

  A separate model in each hospital, so each hospital has its own
  `RR_LOS`. Hospitals below the minimum counts are skipped.

For `"fixed"`, `RxHospital_p` is the F-test p-value for adding a
resistance x hospital interaction, i.e. evidence that the resistance
effect differs between hospitals.

## LOSR and LOSS

Model predictions are averaged over a hospital's own resistant
admissions: `LOSR` with `R = 1` and `LOSS` with `R = 0` for the same
admissions. `LOSS` is therefore the counterfactual LOS of the resistant
admissions had they been susceptible, and `LOSR / LOSS = RR_LOS`. A
hospital with no resistant admissions in a class is averaged over all
its tested admissions instead (`fallback_all_admissions = TRUE`). The
all-hospital values in `rr` are averaged over every resistant admission.

## Confidence intervals

- `"bootstrap"` (default):

  Clusters (`cluster_id_col`, e.g. a patient ID; the admission when
  `NULL`) are resampled with replacement within each hospital, the model
  is refitted, and `RR_LOS`, `LOSR`, `LOSS` and the excess
  (`LOSR - LOSS`, computed within each replicate) are recomputed.
  Intervals are percentiles of the `n_boot` replicates.

- `"delta"`:

  Delta method with a cluster-robust (sandwich) covariance matrix.
  Treats each hospital's admissions as fixed.

- `"none"`:

  No intervals.

`se_log_rr_cluster` (the cluster-robust standard error of \\\beta_R\\)
is always returned. With `return_draws = TRUE` the bootstrap replicates
are returned so they can be combined with profile draws downstream.
