# Estimate Per-Class LOS Relative Risk via Quasi-Poisson Regression

For each pathogen k and antibiotic class c, fits: log(E\[LOS_i\]) = b0 +
b_c \* class_ci + b_HAI \* HAI_i + sum_f gf \* I(centre=f) and returns
RR_kc_LOS = exp(b_c) with 95

## Usage

``` r
fit_los_rr_poisson(
  data,
  patient_id_col = "PatientInformation_id",
  isolate_id_col = "isolate_id",
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
  hai_threshold_hours = 48,
  centre_model = c("pooled_fe", "per_centre"),
  pool_method = c("weighted_mean", "median"),
  max_los = 365,
  min_n = 10L
)
```

## Arguments

- data:

  Data frame. Merged AMR dataset at isolate x antibiotic level.

- patient_id_col:

  Character. Default `"PatientInformation_id"`.

- isolate_id_col:

  Character. Default `"isolate_id"`.

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

  Character or NULL. Filter to syndrome (e.g. "BSI"). NULL = all
  syndromes (universal RR).

- organism_name:

  Character or NULL. Pathogen(s) to process. NULL = all.

- hai_threshold_hours:

  Numeric. HAI/CAI gap threshold. Default `48`.

- centre_model:

  Character. `"pooled_fe"` (default) or `"per_centre"`.

- pool_method:

  Character. Used when centre_model = "per_centre". `"weighted_mean"`
  (default) or `"median"`.

- max_los:

  Numeric. Default `200`.

- min_n:

  Integer. Minimum patients per model. Default `10`.

## Value

Data frame: pathogen, antibiotic_class, RR_LOS, CI_lower, CI_upper,
n_patients, centre_model, syndrome_scope. When centre_model =
"per_centre", attribute "per_centre_rr" contains per-centre estimates.

## Details

LOS is computed with HAI/CAI-specific clock (see compute_patient_los).
Resistance is classified at antibiotic class level (class = 1 if any
drug in class = R). HAI enters as a binary covariate to absorb residual
baseline severity differences beyond the LOS clock adjustment. No
PreDays covariate is used – the LOS clock change from admission to
culture date for HAI patients already removes the pre-infection period
from the outcome.

Stated limitation: when syndrome_name is supplied, RR_kc_LOS is
syndrome-specific. Downstream PAF applies this RR to all profiles of
pathogen k, assuming syndrome-invariant LOS prolongation across
infection sources. Set syndrome_name = NULL for a universal RR.
