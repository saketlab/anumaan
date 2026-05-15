# Fit mortality model and derive adjusted relative risk of death

Fits a pathogen-class specific mixed-effects logistic regression model:

## Usage

``` r
daly_fit_mortality_rr(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  unit_type_col = "unit_type",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  age_col = "Age",
  sex_col = "Gender",
  comorbidity_col = NULL,
  death_value = "Death",
  syndrome_name = NULL,
  organism_name = NULL,
  hai_threshold_hours = 48,
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  phi_threshold = 0.7,
  min_n = 20L,
  min_deaths = 10L,
  use_random_intercept = TRUE
)
```

## Arguments

- data:

  Data frame.

- patient_id_col:

  Character. Patient identifier column.

- facility_col:

  Character or NULL. Facility / centre column. When `NULL`, no hospital
  random effect is used and a single pooled RR is returned using
  standard logistic regression (hospital column will be `NA`). When
  provided, a separate logistic model is fit per hospital and the result
  contains one row per pathogen-class-hospital combination
  (hospital-specific RR).

- organism_col:

  Character. Pathogen column.

- syndrome_col:

  Character. Syndrome column.

- infection_type_col:

  Character. Raw infection-type column.

- antibiotic_class_col:

  Character. Antibiotic class column.

- antibiotic_name_col:

  Character. Antibiotic name column.

- antibiotic_value_col:

  Character. Antibiotic susceptibility column.

- unit_type_col:

  Character or NULL. ICU / ward location column.

- date_admission_col:

  Character. Admission date column.

- date_culture_col:

  Character. First positive culture date column.

- final_outcome_col:

  Character. Final outcome column.

- final_outcome_date_col:

  Character. Final outcome date column.

- age_col:

  Character. Age column.

- sex_col:

  Character. Sex column.

- comorbidity_col:

  Character or NULL. Comorbidity column.

- death_value:

  Character. Value indicating death.

- syndrome_name:

  Character or NULL. Restrict to one syndrome.

- organism_name:

  Character vector or NULL. Restrict to these pathogens.

- hai_threshold_hours:

  Numeric. HAI derivation threshold.

- icu_values:

  Character vector. Values treated as ICU.

- phi_threshold:

  Numeric. HAI / ICU collinearity warning threshold.

- min_n:

  Integer. Minimum patients per fitted model.

- min_deaths:

  Integer. Minimum deaths per fitted model.

- use_random_intercept:

  Logical. Retained for backward compatibility; ignored when
  `facility_col` is provided (per-hospital fitting is used) or `NULL`
  (no facility information available).

## Value

Data frame with hospital (NA when `facility_col` is `NULL`), adjusted OR
and adjusted RR of death.

## Details

Death_i ~ Bernoulli(pi_i) logit(pi_i) = beta0 + beta_kd \*
Resistant_id + beta_age \* Age_i + beta_sex \* Sex_i + beta_hai \*
HAI_i + beta_icu \* ICU_i + beta_comorb \* Comorbidity_i + beta_syn \*
Syndrome_i + u_facility

where u_facility is a random intercept when more than one facility is
present.

The function returns both: - OR_death = exp(beta_kd) - RR_death =
mean(p_resistant_cf) / mean(p_susceptible_cf)

RR_death is derived from model-based predicted probabilities under
resistant and susceptible counterfactual scenarios, as required by the
burden methodology.
