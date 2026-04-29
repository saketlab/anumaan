# globals.R
# Package-level imports and NSE variable declarations

#' @importFrom magrittr %>%
#' @importFrom rlang .data .env := %||%
#' @importFrom stats density median na.omit pnorm profile quantile reorder setNames
#' @importFrom utils adist head object.size
NULL

# Column names used in tidyverse non-standard evaluation (NSE) pipelines
# across the package. These are not actual global variables -- they are
# column references inside dplyr::mutate(), dplyr::filter(), etc.
utils::globalVariables(c(
  # magrittr dot placeholder (not covered by @importFrom)
  ".",

  # Internal temporary columns created in dplyr pipelines
  ".abg_key", ".abx", ".abx_val", ".adm_date", ".chain_seq",
  ".cult_ep_min", ".disc", ".dt", ".dup_rank", ".ep_id", ".ep_idx",
  ".ep_start", ".inf_raw", ".is_dup", ".is_icu", ".is_missing",
  ".new_ep", ".org", ".out_date", ".prov_key", ".pt", ".raw_cult",
  ".spc", ".temp_organism", ".unit_norm", ".val",

  # Column names referenced via NSE in dplyr/tidyr pipelines
  "Age", "Antibiotic", "CFR_L", "Category", "Class", "Common",
  "commensals", "DOB", "DS_J", "D_J", "D_L_contribution", "D_x_L",
  "HAI", "ICU", "LOS_days", "M_LJ", "N_F_L", "N_F_LK", "N_NF_L",
  "N_NF_LK", "N_resistant", "N_tested", "RR", "S_J", "Syndrome",
  "Type of culture/specimen", "YLL_Kdelta",
  "age_bin_gbd", "age_confidence", "age_diff", "age_method",
  "antibiotic_class", "antibiotic_value", "calculated_age",
  "calculated_los", "cases_h", "class_resistance", "class_result",
  "class_result_event", "comorbidity_encoded", "completeness_pct",
  "contaminant_confidence", "contaminant_method", "count",
  "Age_model", "Sex_model",
  "date_of_admission", "date_of_culture", "date_of_final_outcome",
  "days_to_culture", "death", "death_weight", "deaths", "deaths_h",
  "department_confidence", "device_insertion_date", "discharged_h",
  "drug", "drug_class", "effective_DW", "event_id", "facility",
  "fatal", "final_outcome", "final_result", "gap_days",
  "has_ip", "has_op", "has_resistant", "hierarchy_rank",
  "infection_deaths", "infection_deaths_J", "infection_deaths_LJ",
  "infection_type_derived", "infection_type_method", "inferred_dept",
  "inferred_type", "latitude", "le_both", "life_expectancy",
  "life_expectancy_years", "location_name", "longitude",
  "m_r", "manual_weight", "marginal_resistance", "mdr",
  "mdr_confidence", "mdr_method", "mdr_threshold",
  "n", "n_classes", "n_deaths", "n_deaths_total", "n_drugs_in_class",
  "n_events", "n_inpatient", "n_mono", "n_op_patients", "n_op_to_ip",
  "n_patients",
  "PAF_k_mort", "PAF_kd", "PAF_mortality",
  "YLL_associated_k", "YLL_attributable_contribution", "YLL_attributable_k",
  "n_organisms", "n_records", "n_resistant", "n_resistant_categories",
  "n_resistant_classes", "n_susceptible_categories", "n_tested",
  "n_total", "organism_group", "organism_normalized", "pathogen",
  "patient_id", "patient_profile", "percentage",
  "polymicrobial_weight", "proportion",
  "ref_lower", "resistance_profile", "resistance_rank", "resistant",
  "resistant_categories", "resistant_flag", "result_clean",
  "result_count", "rr_drug", "rr_drug_norm", "rr_pathogen",
  "rr_pathogen_norm", "selection_method", "sex",
  "temp_drug_ref_lower", "temp_path_lower", "tested", "total",
  "total_mono", "total_tested", "total_weight",
  "value", "weight", "weighted_deaths", "weighted_total",
  "xdr", "xdr_confidence", "xdr_method", "yll_contribution",

  # Column names with spaces/trailing whitespace from CSV headers
  "Common commensals", "Category ",

  # Temporary dplyr/mutate columns not yet in the list above
  ".gap_days", ".prev_adm", ".present", ".syndrome_rank", ".norm_syn",
  "syndrome_selected", "poly_weight",
  "n_centres_present", "column",

  # Columns created in prep_map_diagnosis_to_icd / alethia integration
  "given_entity", "alethia_prediction", "alethia_score",
  "diagnosis_text", "icd_score", "icd_rank", "icd_prediction",
  "icd_code", "icd_method",

  "final_abx", "grand_total", "value_label", "percent", "label",
  "med_age", "med_los", "org_ordered", "organism_resist",
  "total_patients", "unique_patients", "unique_organisms",
  "infection_type", "location_type", "tests", "total_tests",
  "patients", "first",


  "year", "total_year",


  "center", "metric_type", "organism",
  "associated", "attributable",
  "yll_per_1000", "yld_per_1000",
  "yll_sum", "tot", "cell_lbl"
))
