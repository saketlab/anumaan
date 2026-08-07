# Fit relative LOS using Gamma GLM with log link

Estimates class-specific relative length of stay (RR_LOS) for each
pathogen using a Gamma generalized linear model with log link:

## Usage

``` r
daly_fit_los_rr(
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
  age_col = "Age",
  sex_col = "Gender",
  unit_type_col = NULL,
  comorbidity_col = NULL,
  syndrome_name = NULL,
  organism_name = NULL,
  hai_threshold_hours = 48,
  max_los = 365,
  min_n = 20L,
  min_resistant = 5L,
  min_susceptible = 5L
)
```

## Arguments

- data:

  Data frame.

- patient_id_col:

  Character. Patient identifier column.

- facility_col:

  Character or NULL. Facility / centre column. When `NULL`, no hospital
  effect is included in the model and a single pooled RR is returned
  (hospital column will be `NA`). When provided, the model is fit
  separately per hospital and the result contains one row per
  pathogen-class-hospital combination (hospital-specific RR).

- organism_col:

  Character. Pathogen column.

- syndrome_col:

  Character. Syndrome column.

- infection_type_col:

  Character. Raw infection type column.

- antibiotic_class_col:

  Character. Antibiotic class column.

- antibiotic_name_col:

  Character. Antibiotic name column.

- antibiotic_value_col:

  Character. Antibiotic susceptibility column.

- date_admission_col:

  Character. Admission date column.

- date_discharge_col:

  Character. Discharge / final outcome date column.

- date_culture_col:

  Character. First positive culture date column.

- final_outcome_col:

  Character. Final outcome column.

- final_outcome_value:

  Character. Value indicating discharged patients.

- age_col:

  Character. Age column.

- sex_col:

  Character. Sex column.

- unit_type_col:

  Character or NULL. ICU / ward location column.

- comorbidity_col:

  Character or NULL. Comorbidity column.

- syndrome_name:

  Character or NULL. Restrict to one syndrome.

- organism_name:

  Character vector or NULL. Restrict to these pathogens.

- hai_threshold_hours:

  Numeric. HAI derivation threshold.

- max_los:

  Numeric. Maximum retained LOS in days.

- min_n:

  Integer. Minimum patients required in a class model.

- min_resistant:

  Integer. Minimum resistant patients required.

- min_susceptible:

  Integer. Minimum susceptible patients required.

## Value

Data frame with pathogen, antibiotic_class, hospital (NA when
`facility_col` is `NULL`), RR_LOS, CI_lower, CI_upper, beta_log_rr,
mean_los_resistant, mean_los_susceptible, n_patients, n_resistant,
n_susceptible, model_family, link_function, and syndrome_scope.

## Details

E(LOS_i) = exp(beta0 + beta1\*Resistance_ci + beta2\*HAI_i +
beta3\*Age_i + beta4\*Sex_i + beta5\*ICU_i + beta6\*Comorbidity_i +
centre fixed effects)

RR_LOS for a resistance class is exp(beta1).
