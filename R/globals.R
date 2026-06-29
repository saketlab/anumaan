# globals.R
# Package-level imports and NSE variable declarations

#' @importFrom grDevices dev.off pdf
#' @importFrom magrittr %>%
#' @importFrom rlang .data .env := %||%
#' @importFrom stats cor density median na.omit pnorm profile quantile rbinom reorder setNames
#' @importFrom utils adist combn head object.size
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
  "final_outcome_date", "diagnosis_time_admission",
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
  "previous_hospitalisation",
  "polymicrobial_weight", "proportion",
  "ref_lower", "resistance_profile", "resistance_rank", "resistant",
  "resistant_categories", "resistant_flag", "result_clean",
  "result_count", "rr_drug", "rr_drug_norm", "rr_pathogen",
  "rr_pathogen_norm", "selection_method", "sex",
  "temp_drug_ref_lower", "temp_path_lower", "tested", "total",
  "total_mono", "total_tested", "total_weight",
  "value", "weight", "weighted_deaths", "weighted_total",
  "xdr", "xdr_confidence", "xdr_method", "yll_contribution",
  "Age_bin", "Age_years", "center_name",
  "DALY_associated", "DALY_assoc_per_1000",
  "DALY_attributable", "DALY_attr_per_1000",
  "org_group", "profile",
  "YLD_associated", "YLD_assoc_per_1000",
  "YLD_attributable", "YLD_attr_per_1000",
  "YLD_base", "YLD_base_per_1000",
  "YLL_associated", "YLL_assoc_per_1000",
  "YLL_attributable", "YLL_attr_per_1000",
  "YLL_base", "YLL_base_per_1000",

  # Column names with spaces/trailing whitespace from CSV headers
  "Common commensals", "Category ",

  # Temporary dplyr/mutate columns not yet in the list above
  ".gap_days", ".prev_adm", ".present", ".syndrome_rank", ".norm_syn",
  "syndrome_selected", "poly_weight",
  "n_centres_present", "column", "n_ast_values",

  # Columns created in prep_map_diagnosis_to_icd / alethia integration
  "given_entity", "alethia_prediction", "alethia_score",
  "diagnosis_text", "icd_score", "icd_rank", "icd_prediction",
  "icd_code", "icd_method",
  # prep_eda_plots.R -- dplyr NSE column references
  "final_abx", "grand_total", "value_label", "percent", "label",
  "med_age", "med_los", "org_ordered", "organism_resist",
  "total_patients", "unique_patients", "unique_organisms",
  "infection_type", "location_type", "tests", "total_tests",
  "patients", "first",
  "year", "total_year",
  "center", "metric_type", "organism",
  "associated", "attributable",
  "yll_per_1000", "yld_per_1000",
  "yll_sum", "tot", "cell_lbl",

  # prep_standardize_specimens.R -- temporary/derived column names
  "temp_spec_input", "temp_spec_clean",
  "specimen_normalized", "sample_category", "sterile_classification",
  "sample_type", "culture_date", "admission_date", "outcome_date",
  "unit_admission_date", "unit_duration_days", "location",

  # plot_resistance_by_agebin / plot_resistance_by_organism -- dplyr NSE columns
  "abx_call", "pct", "pt_count", "resistance",

  # profile_convex.R -- Pathway 1 pipeline NSE column names
  # Temporary column created then removed in the isolate dedup step
  ".sir_rank",
  # Marginal computation: source flag and temporary external-override column
  "marginal_source", ".ext_rate",
  # check_profile_constraints(): bare column refs inside dplyr::summarise
  "pass", "abs_residual",
  # estimate_profiles_convex() / output schema columns created in dplyr::mutate
  "profile_probability", "convergence_flag", "identifiability_flag",
  "max_abs_residual", "profile_set_type", "profile_class_set", "estimator",
  # compute_pairwise_from_data() output columns
  "antibiotic_class_1", "antibiotic_class_2", "pairwise_prevalence",
  "n_co_tested", "rho",
  # bootstrap_profiles_convex() output columns
  "probability_mean", "probability_median",
  "n_replicates_converged", "convergence_rate",
  # validate_profile_inputs() / preprocess_for_profiles() report columns
  "check_name", "n_affected", "n_rows_in", "n_rows_out",
  # compute_resistance_profiles() new stored output fields
  "constraint_targets", "constraint_names",

  # Pathway 2 — fit_bayesian_multivariate_probit()
  "resistance_binary", "antibiotic_class", ".eid", ".pt_adm_key",
  "d_idx", "h_idx", "p_idx",

  # Pathway 2 — compute_event_profile_probabilities()
  "profile_class_set", "profile_delta", "profile_probability",
  "draw_s", "R_ALL_s", "R_NF_s",

  # Pathway 2 — aggregate_profiles_for_daly()
  "R_ALL_mean", "R_ALL_lower", "R_ALL_upper",
  "R_NF_mean",  "R_NF_lower",  "R_NF_upper",
  "used_for_YLL", "used_for_YLD", "profile_set_type",
  "n_draws_used", "profile_label",

  # Pathway 2 — summarize_reserve_drugs()
  "n_total", "n_tested", "n_resistant",
  "tested_subset_resistance", "adjusted_testing_support",
  "subgroup_notes", "reserve_drug_flag", "sensitivity_flag",

  # estimate_resistance_profiles() diagnostics
  "pct_converged", "any_underdetermined", "n_profiles_estimated",
  "n_eligible_pathogens", "n_profiles",

  # aggregate_profiles_for_daly() -- created via dplyr::mutate then selected
  "low_draws_flag",

  # plot_probit_diagnostics() -- ggplot2 aes() column references
  "rhat", "ess", "metric"
))
