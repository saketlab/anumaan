# Estimate Per-Class Mortality Odds Ratio via Mixed-Effects Logistic Regression

For each pathogen k and antibiotic class c, fits:

      logit(P(death_i)) = b0 + b1*Resistance_ci + b2*Age_i + b3*Sex_i
                        + b4*HAI_i + b5*ICU_i + b6*Comorbidity_i
                        + u_facility   [u_facility ~ N(0, sigma^2)]

and returns `OR_death_kc = exp(b1)` with 95% Wald CI.

For each pathogen x antibiotic class combination, fits a mixed-effects
logistic regression (facility random intercept) with resistance, age,
sex, and HAI as fixed effects. Returns a data frame of mortality odds
ratios (OR) and 95

## Usage

``` r
fit_mortality_rr_logistic(
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
  facility_name = NULL,
  hai_threshold_hours = 48,
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  phi_threshold = 0.7,
  min_n = 10L,
  min_deaths = 5L
)

fit_mortality_rr_logistic(
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
  facility_name = NULL,
  hai_threshold_hours = 48,
  icu_values = c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU"),
  phi_threshold = 0.7,
  min_n = 10L,
  min_deaths = 5L
)
```

## Arguments

- data:

  Data frame. Patient-level microbiology records.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- facility_col:

  Character. Default `"center_name"`.

- organism_col:

  Character. Default `"organism_name"`.

- syndrome_col:

  Character. Default `"syndrome"`.

- infection_type_col:

  Character. Default `"type_of_infection"`.

- antibiotic_class_col:

  Character. Default `"antibiotic_class"`.

- antibiotic_name_col:

  Character. Default `"antibiotic_name"`.

- antibiotic_value_col:

  Character. Default `"antibiotic_value"`.

- unit_type_col:

  Character. Default `"unit_type"`.

- date_admission_col:

  Character. Default `"date_of_admission"`.

- date_culture_col:

  Character. Default `"date_of_first_positive_culture"`.

- final_outcome_col:

  Character. Default `"final_outcome"`.

- final_outcome_date_col:

  Character. Default `"final_outcome_date"`.

- age_col:

  Character. Default `"Age"`.

- sex_col:

  Character. Default `"Gender"`.

- comorbidity_col:

  Character or `NULL`. Default `NULL`.

- death_value:

  Character. Default `"Death"`.

- syndrome_name:

  Character or `NULL`. Filter to one syndrome.

- organism_name:

  Character vector or `NULL`. Filter to pathogen(s).

- facility_name:

  Character or `NULL`. If provided, filters data to the specified
  facility before fitting.

- hai_threshold_hours:

  Numeric. Default `48`.

- icu_values:

  Character vector. ICU location labels.

- phi_threshold:

  Numeric. HAI/ICU collinearity threshold. Default `0.7`.

- min_n:

  Integer. Minimum patients per model. Default `10L`.

- min_deaths:

  Integer. Minimum deaths per model. Default `5L`.

## Value

Data frame with one row per (pathogen, antibiotic_class): `pathogen`,
`antibiotic_class`, `OR_death`, `CI_lower`, `CI_upper`, `n_patients`,
`n_deaths`, `convergence_warning`, `syndrome_scope`,
`comorbidity_encoding`.

Data frame with columns: `pathogen`, `antibiotic_class`, `OR_death`,
`CI_lower`, `CI_upper`, `n_patients`, `n_deaths`, `convergence_warning`,
`syndrome_scope`, `comorbidity_encoding`.

## Details

**Resistance** is classified at antibiotic class level (class = 1 if any
drug in the class is R), matching
[`.build_class_resistance_wide()`](https://saketlab.github.io/anumaan/reference/dot-build_class_resistance_wide.md).

**HAI/CAI** is derived by
[`derive_infection_type_for_mortality()`](https://saketlab.github.io/anumaan/reference/derive_infection_type_for_mortality.md),
which prioritises explicit labels over date-gap inference. Patients
labelled "Not Known" are excluded from each model (missing HAI
covariate).

**ICU flag** is derived using the "ever in ICU" rule by
[`.derive_icu_binary()`](https://saketlab.github.io/anumaan/reference/dot-derive_icu_binary.md):
a patient is flagged ICU = 1 if any row for that admission records an
ICU unit type.

**Comorbidity** is standardised by
[`.encode_comorbidity_mortality()`](https://saketlab.github.io/anumaan/reference/dot-encode_comorbidity_mortality.md)
(numeric Charlson, binary present/none, or ordinal
none/mild/moderate/severe).

**Pre-fitting checks** per pathogen x class:

- `min_n` patients with complete required covariates.

- `min_deaths` deaths (ensures outcome variation).

- Variation in resistance (both R and S patients present).

- HAI/ICU collinearity (phi coefficient; `phi_threshold`).
