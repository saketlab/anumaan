# daly_rr_and_los.R
# LOS relative risk, mortality OR, and resistance-profile fitting functions

# Theoretical median of a fitted LOS distribution.
#   Weibull (shape k, scale lambda):   lambda * log(2)^(1/k)
#   Log-normal (meanlog mu, sdlog s):  exp(mu)
#   Gamma (shape a, rate b):           no closed form -> qgamma(0.5, a, b)
# Returns NA for unsupported distributions, invalid parameters, or errors.
.compute_median_from_fit <- function(fit, dist) {
  tryCatch(
    {
      med <- if (dist == "weibull") {
        k <- unname(fit$estimate["shape"])
        l <- unname(fit$estimate["scale"])
        if (!is.finite(k) || !is.finite(l) || k <= 0 || l <= 0) NA_real_
        else l * log(2)^(1 / k)
      } else if (dist == "lnorm") {
        mu <- unname(fit$estimate["meanlog"])
        sg <- unname(fit$estimate["sdlog"])
        if (!is.finite(mu) || !is.finite(sg) || sg <= 0) NA_real_
        else exp(mu)
      } else if (dist == "gamma") {
        a <- unname(fit$estimate["shape"])
        r <- unname(fit$estimate["rate"])
        if (!is.finite(a) || !is.finite(r) || a <= 0 || r <= 0) NA_real_
        else stats::qgamma(0.5, shape = a, rate = r)
      } else {
        NA_real_
      }
      if (length(med) != 1L || !is.finite(med)) NA_real_ else med
    },
    error = function(e) NA_real_
  )
}

#' @keywords internal
.safe_fitdist <- function(x, dist) safe_fit(x, dist)

# Fit candidate distributions to one LOS vector, select the lowest AIC, and
# return the theoretical median of the selected distribution.
# Non-finite and non-positive LOS values are dropped first. Returns NA
# median/dist/aic when nothing can be fitted.
.fit_best_dist_median <- function(los_vec, distributions) {
  out_na <- list(median = NA_real_, dist = NA_character_, aic = NA_real_)
  los_vec <- los_vec[is.finite(los_vec) & los_vec > 0]
  if (length(los_vec) == 0L) {
    return(out_na)
  }
  fits <- Filter(
    Negate(is.null),
    stats::setNames(lapply(distributions, .safe_fitdist, x = los_vec), distributions)
  )
  if (length(fits) == 0L) {
    return(out_na)
  }
  aics <- vapply(fits, function(f) as.numeric(f$aic), numeric(1))
  best <- names(which.min(aics))
  list(
    median = .compute_median_from_fit(fits[[best]], best),
    dist = best,
    aic = unname(aics[[best]])
  )
}


# Step 1b
#' Compute Patient-Level Length of Stay
#'
#' Computes hospital length of stay with a single definition for every
#' patient, regardless of infection type (HAI/CAI) or culture date:
#' \deqn{LOS = final\ outcome\ (discharge)\ date - admission\ date}
#' Only patients whose final outcome equals \code{final_outcome_value}
#' (default \code{"Discharged"}) are retained. LOS is \code{NA} when either
#' date is missing; such patients are excluded. A checkpoint reports how many
#' discharged patients are missing a discharge date. Rows with LOS <= 0 or
#' LOS > \code{max_los} are dropped.
#'
#' @param data Data frame.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_discharge_col Character. Final outcome / discharge date column.
#'   Default \code{"final_outcome_date"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_value Character. Default \code{"Discharged"}.
#' @param syndrome_col Character. Syndrome column name. Only used when
#'   \code{syndrome_name} is not \code{NULL}. Default \code{"syndrome"}.
#' @param syndrome_name Character or \code{NULL}. If provided, only patients
#'   with this syndrome are retained before LOS computation. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. If provided, filters data
#'   to the specified facility before computing LOS. Default \code{NULL}.
#' @param max_los Numeric. Upper cap on LOS in days; patients exceeding this
#'   are excluded. Default \code{200}.
#'
#' @return Data frame with one row per distinct patient x facility x LOS
#'   (\code{patient_id_col}, \code{facility_col}, \code{LOS_days}).
#' @export
daly_compute_patient_los <- function(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
  syndrome_col = "syndrome",
  syndrome_name = NULL,
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  max_los = 200
) {
  # Column validation
  required <- c(
    patient_id_col, facility_col,
    date_admission_col, date_discharge_col,
    final_outcome_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  # Facility filtering
  if (!is.null(facility_name)) {
    avail <- unique(data[[facility_col]])
    if (!facility_name %in% avail) {
      stop(sprintf(
        "facility_name '%s' not found in column '%s'. Available: %s",
        facility_name, facility_col,
        paste(sort(as.character(avail)), collapse = ", ")
      ))
    }
    data <- data[data[[facility_col]] == facility_name, ]
    message(sprintf(
      "Filtering to facility '%s': %d row(s) retained.",
      facility_name, nrow(data)
    ))
  } else {
    message(sprintf(
      "Computing LOS across all facilities; '%s' column included in output.",
      facility_col
    ))
  }

  # Syndrome filtering
  if (!is.null(syndrome_name)) {
    avail_syn <- unique(data[[syndrome_col]])
    if (!.name_key(syndrome_name) %in% .name_key(avail_syn)) {
      stop(sprintf(
        "syndrome_name '%s' not found in column '%s'. Available: %s",
        syndrome_name, syndrome_col,
        paste(sort(as.character(avail_syn)), collapse = ", ")
      ))
    }
    data <- data[.name_key(data[[syndrome_col]]) == .name_key(syndrome_name), ]
    message(sprintf(
      "Filtering to syndrome '%s': %d row(s) retained.",
      syndrome_name, nrow(data)
    ))
  }

  # Checkpoint: discharge date coverage
  discharged_data <- data[data[[final_outcome_col]] == final_outcome_value, ]
  n_total <- nrow(discharged_data)
  n_no_disc <- sum(is.na(as.Date(discharged_data[[date_discharge_col]])))

  if (n_no_disc > 0L) {
    message(sprintf(
      "Checkpoint: %d of %d discharged patient(s) are missing a discharge date.",
      n_no_disc, n_total
    ))
    message(
      "  -> Patients with a missing discharge date are excluded (LOS = NA)."
    )
  }

  # LOS computation: one definition for every patient, independent of
  # infection type (HAI/CAI) and culture date.
  # LOS = final outcome (discharge) date - admission date; NA if either is missing.
  los_data <- data %>%
    dplyr::filter(.data[[final_outcome_col]] == final_outcome_value) %>%
    dplyr::mutate(
      .adm = as.Date(.data[[date_admission_col]]),
      .disc = as.Date(.data[[date_discharge_col]]),
      LOS_days = as.numeric(.disc - .adm)
    )

  # Final filtering and deduplication
  los_data <- los_data %>%
    dplyr::select(-".adm", -".disc") %>%
    dplyr::filter(!is.na(LOS_days), LOS_days > 0L, LOS_days <= max_los) %>%
    dplyr::distinct(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(facility_col),
      LOS_days
    )

  message(sprintf(
    "Patient LOS computed: %d patients retained (LOS > 0 and <= %d days).",
    nrow(los_data), max_los
  ))
  return(los_data)
}


# Step 1c (internal)
#' Collapse drug-level data to antibiotic-class binary wide matrix
#'
#' Standardises antibiotic values (Intermediate -> R/S per GBD rules),
#' collapses to class level (class = 1 if ANY drug in class = R), and pivots
#' wide: one row per patient, one column per antibiotic class (0/1).
#' Column names are sanitised with make.names(); the original-to-safe name
#' mapping is stored as the "class_name_map" attribute.
#'
#' @noRd
.build_class_resistance_wide <- function(
  data,
  patient_id_col = "PatientInformation_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  untested_fill = 0L
) {
  abx_std <- data %>%
    dplyr::filter(
      !is.na(.data[[antibiotic_value_col]]),
      !is.na(.data[[antibiotic_class_col]])
    ) %>%
    dplyr::mutate(
      .abx_name = stringr::str_to_lower(stringr::str_trim(
        as.character(.data[[antibiotic_name_col]])
      )),
      .abx_val = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[antibiotic_value_col]])
      )),
      .abx_val = dplyr::case_when(
        .abx_name == "colistin" & .abx_val == "I" ~ "S",
        .abx_val == "I" ~ "R",
        .abx_val %in% c("R", "RESISTANT") ~ "R",
        .abx_val %in% c("S", "SUSCEPTIBLE", "SENSITIVE") ~ "S",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.abx_val))

  class_long <- abx_std %>%
    dplyr::group_by(
      .data[[patient_id_col]],
      .data[[antibiotic_class_col]]
    ) %>%
    dplyr::summarise(
      class_R = as.integer(any(.abx_val == "R")),
      .groups = "drop"
    )

  orig_classes <- sort(unique(class_long[[antibiotic_class_col]]))
  safe_classes <- make.names(orig_classes)
  class_name_map <- setNames(orig_classes, safe_classes)

  class_wide <- class_long %>%
    dplyr::mutate(
      .class_safe = make.names(.data[[antibiotic_class_col]])
    ) %>%
    tidyr::pivot_wider(
      id_cols     = !!rlang::sym(patient_id_col),
      names_from  = ".class_safe",
      values_from = "class_R",
      values_fill = untested_fill
    )

  attr(class_wide, "class_name_map") <- class_name_map
  return(class_wide)
}


# daly_fit_los_rr
#' Fit Relative Length of Stay with a Gamma GLM (Log Link)
#'
#' For each pathogen and antibiotic class, fits a Gamma generalised linear
#' model with log link to the length of stay (LOS) of discharged admissions:
#'
#' \deqn{\log E(LOS_i) = \beta_0 + \beta_R R_i + \beta_1 HAI_i + \beta_2 Age_i
#'   + \beta_3 Sex_i + \beta_4 ICU_i + \beta_5 Comorbidity_i + \gamma_{h(i)}}
#'
#' and returns the relative LOS, \eqn{RR_{LOS} = \exp(\beta_R)}, together with
#' the adjusted LOS of resistant admissions (\code{LOSR}) and of the same
#' admissions had they been susceptible (\code{LOSS}), per hospital and for
#' all hospitals combined.
#'
#' @section Unit of analysis:
#' One row per admission x pathogen x antibiotic class. LOS is computed once
#' per admission by \code{daly_compute_patient_los()} (final outcome date
#' minus admission date, discharged admissions only, 0 < LOS <= \code{max_los}).
#' An admission with several pathogens enters each pathogen's models.
#'
#' @section Resistance (worst phenotype):
#' For an admission and pathogen, a class is resistant (\code{R = 1}) when any
#' result of that pathogen for any drug in the class is R, and susceptible
#' (\code{R = 0}) when it was tested and no result is R. Results of other
#' pathogens in the same admission are not used. With
#' \code{intermediate_as_resistant = TRUE}, I counts as R except for colistin,
#' where I counts as S (the convention of
#' \code{prep_recode_intermediate_ast()}); with \code{FALSE}, I results are
#' dropped. Admissions not tested for a class are excluded from that class's
#' model.
#'
#' @section Covariates:
#' Taken from each admission's earliest culture of the pathogen
#' (\code{date_culture_col}): HAI (1 if \code{infection_type_col} equals
#' \code{hai_value}, 0 if it equals \code{cai_value}, otherwise missing), age
#' (linear), sex, ICU (when \code{unit_type_col} is given) and comorbidity
#' (when \code{comorbidity_col} is given). Admissions with a missing value in
#' a requested covariate are dropped from that model and counted in
#' \code{n_dropped_covariates}. A covariate with no variation in a model is
#' left out of that model.
#'
#' @section Hospital effect:
#' \describe{
#'   \item{\code{"fixed"} (default)}{One model per pathogen x class across all
#'     hospitals, with hospital as a fixed effect. Each hospital has its own
#'     LOS level; the resistance effect (\code{RR_LOS}) is shared. Small
#'     hospitals are kept and borrow the resistance effect.}
#'   \item{\code{"none"}}{One model with no hospital term.}
#'   \item{\code{"stratified"}}{A separate model in each hospital, so each
#'     hospital has its own \code{RR_LOS}. Hospitals below the minimum counts
#'     are skipped.}
#' }
#' For \code{"fixed"}, \code{RxHospital_p} is the F-test p-value for adding a
#' resistance x hospital interaction, i.e. evidence that the resistance
#' effect differs between hospitals.
#'
#' @section LOSR and LOSS:
#' Model predictions are averaged over a hospital's own resistant admissions:
#' \code{LOSR} with \code{R = 1} and \code{LOSS} with \code{R = 0} for the
#' same admissions. \code{LOSS} is therefore the counterfactual LOS of the
#' resistant admissions had they been susceptible, and
#' \code{LOSR / LOSS = RR_LOS}. A hospital with no resistant admissions in a
#' class is averaged over all its tested admissions instead
#' (\code{fallback_all_admissions = TRUE}). The all-hospital values in
#' \code{rr} are averaged over every resistant admission.
#'
#' @section Confidence intervals:
#' \describe{
#'   \item{\code{"bootstrap"} (default)}{Clusters (\code{cluster_id_col}, e.g.
#'     a patient ID; the admission when \code{NULL}) are resampled with
#'     replacement within each hospital, the model is refitted, and
#'     \code{RR_LOS}, \code{LOSR}, \code{LOSS} and the excess
#'     (\code{LOSR - LOSS}, computed within each replicate) are recomputed.
#'     Intervals are percentiles of the \code{n_boot} replicates.}
#'   \item{\code{"delta"}}{Delta method with a cluster-robust (sandwich)
#'     covariance matrix. Treats each hospital's admissions as fixed.}
#'   \item{\code{"none"}}{No intervals.}
#' }
#' \code{se_log_rr_cluster} (the cluster-robust standard error of
#' \eqn{\beta_R}) is always returned. With \code{return_draws = TRUE} the
#' bootstrap replicates are returned so they can be combined with profile
#' draws downstream.
#'
#' @param data Data frame with one row per antibiotic test result.
#' @param admission_id_col Character. Admission identifier; the unit of
#'   analysis. Default \code{"PatientInformation_id"}.
#' @param cluster_id_col Character or \code{NULL}. Identifier used to cluster
#'   repeated admissions of the same patient in the bootstrap and the
#'   sandwich covariance. \code{NULL} clusters by admission.
#' @param facility_col Character or \code{NULL}. Hospital column. \code{NULL}
#'   pools everything and forces \code{hospital_effect = "none"}.
#' @param organism_col Character. Pathogen column.
#' @param syndrome_col Character. Syndrome column (used with
#'   \code{syndrome_name}).
#' @param infection_type_col Character. Column holding the final HAI/CAI
#'   classification. Default \code{"infection_type"}.
#' @param antibiotic_class_col,antibiotic_name_col,antibiotic_value_col
#'   Character. Antibiotic class, drug name and S/I/R result columns.
#' @param date_admission_col,date_discharge_col Character. Admission and final
#'   outcome date columns.
#' @param date_culture_col Character. Culture date column, used to take
#'   covariates from each admission's earliest culture of the pathogen. When
#'   absent, the first row is used.
#' @param final_outcome_col,final_outcome_value Character. Final outcome column
#'   and the value marking discharged admissions.
#' @param age_col,sex_col Character or \code{NULL}. Age and sex columns;
#'   \code{NULL} leaves the covariate out.
#' @param unit_type_col Character or \code{NULL}. Ward / ICU column.
#' @param comorbidity_col Character or \code{NULL}. Comorbidity column.
#' @param syndrome_name Character or \code{NULL}. Restrict to one syndrome.
#' @param organism_name Character vector or \code{NULL}. Pathogens to model;
#'   \code{NULL} models all.
#' @param classes \code{NULL} (all classes), a character vector (the same
#'   classes for every pathogen), or a named list of character vectors (a
#'   panel per pathogen, names matching \code{organism_name} values).
#' @param hospital_effect One of \code{"fixed"}, \code{"none"},
#'   \code{"stratified"}. See "Hospital effect".
#' @param reference_hospital Character or \code{NULL}. Reference level of the
#'   hospital fixed effect; \code{NULL} uses the hospital with most
#'   admissions in each model. Does not change any estimate.
#' @param intermediate_as_resistant Logical. See "Resistance".
#' @param hai_value,cai_value Character. Values of \code{infection_type_col}
#'   for HAI and CAI; other values are treated as missing.
#' @param icu_values Character vector. Values of \code{unit_type_col} counted
#'   as ICU (matched case-insensitively as substrings).
#' @param max_los Numeric. Maximum LOS in days. Default 365.
#' @param min_n,min_resistant,min_susceptible Integer. Minimum admissions,
#'   resistant admissions and susceptible admissions for a model to be fitted.
#' @param flag_min Integer. Hospitals with fewer resistant or susceptible
#'   admissions than this are flagged \code{low_support} in the hospital
#'   table (they are still estimated).
#' @param ci_method One of \code{"bootstrap"}, \code{"delta"}, \code{"none"}.
#' @param n_boot Integer. Bootstrap replicates. Default 1000.
#' @param conf_level Numeric. Confidence level. Default 0.95.
#' @param seed Integer or \code{NULL}. Random seed for the bootstrap.
#' @param return_draws Logical. Return the bootstrap replicates.
#' @param patient_id_col Deprecated. Use \code{admission_id_col}.
#' @param hai_threshold_hours Deprecated and ignored. HAI/CAI is read from
#'   \code{infection_type_col}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{rr}}{One row per pathogen x class (x hospital for
#'     \code{"stratified"}): \code{RR_LOS} with \code{CI_lower},
#'     \code{CI_upper}, \code{beta_log_rr}, \code{se_log_rr_cluster}; counts;
#'     observed mean LOS of R and S admissions; all-hospital \code{LOSR},
#'     \code{LOSS} and \code{excess_LOS} with intervals; \code{dispersion};
#'     \code{RxHospital_p}; covariates used and settings.}
#'   \item{\code{hospital}}{One row per pathogen x class x hospital:
#'     counts, observed and predicted mean LOS of R and S admissions (a
#'     calibration check), \code{LOSR}, \code{LOSS}, \code{excess_LOS} with
#'     intervals, \code{RR_LOS}, \code{low_support} and
#'     \code{fallback_all_admissions}.}
#'   \item{\code{draws}}{When \code{return_draws = TRUE} and
#'     \code{ci_method = "bootstrap"}: a list named
#'     \code{"pathogen | class"} (\code{"pathogen | class | hospital"} for
#'     \code{"stratified"}), each with \code{RR_LOS} (vector),
#'     \code{LOSR}, \code{LOSS} (replicate x hospital matrices, days) and
#'     \code{LOSR_all}, \code{LOSS_all} (all-hospital values, vectors). Pass
#'     it to \code{daly_calc_yld_associated()} /
#'     \code{daly_calc_yld_attributable()} as \code{los_draws}.}
#'   \item{\code{skipped}}{Pathogen x class (x hospital) combinations not
#'     fitted, with the reason.}
#' }
#' @export
daly_fit_los_rr <- function(
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
) {
  hospital_effect <- match.arg(hospital_effect)
  ci_method <- match.arg(ci_method)
  if (!is.null(patient_id_col)) {
    warning("`patient_id_col` is deprecated; use `admission_id_col`. Using it as the admission ID.", call. = FALSE)
    admission_id_col <- patient_id_col
  }
  if (!is.null(hai_threshold_hours)) {
    warning("`hai_threshold_hours` is deprecated and ignored; HAI/CAI is read from `infection_type_col`.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  required <- c(
    admission_id_col, organism_col, infection_type_col,
    antibiotic_class_col, antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col, final_outcome_col,
    facility_col, cluster_id_col, age_col, sex_col, unit_type_col, comorbidity_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")))
  }

  data <- as.data.frame(data)
  if (is.null(facility_col)) {
    facility_col <- ".facility_pooled"
    data[[facility_col]] <- "pooled"
    if (hospital_effect != "none") {
      message("No facility_col provided: hospital_effect set to 'none'.")
      hospital_effect <- "none"
    }
  }
  if (!is.null(syndrome_name)) {
    data <- data[.name_key(data[[syndrome_col]]) == .name_key(syndrome_name), , drop = FALSE]
  }
  data <- data[!is.na(data[[organism_col]]) & trimws(as.character(data[[organism_col]])) != "", , drop = FALSE]

  # LOS: one value per discharged admission
  los_df <- suppressMessages(daly_compute_patient_los(
    data = data, patient_id_col = admission_id_col, facility_col = facility_col,
    date_admission_col = date_admission_col, date_discharge_col = date_discharge_col,
    final_outcome_col = final_outcome_col, final_outcome_value = final_outcome_value,
    max_los = max_los
  ))
  dup <- duplicated(los_df[[admission_id_col]])
  if (any(dup)) {
    warning(sprintf(
      "%d admission(s) have more than one LOS (inconsistent dates or hospital); the first is kept.",
      sum(dup)
    ), call. = FALSE)
    los_df <- los_df[!dup, , drop = FALSE]
  }
  if (nrow(los_df) == 0L) {
    warning("No valid LOS observations found after filtering.", call. = FALSE)
    return(.los_rr_empty())
  }
  message(sprintf("LOS available for %d discharged admission(s).", nrow(los_df)))
  los_tbl <- data.frame(
    .id = as.character(los_df[[admission_id_col]]),
    .hosp = as.character(los_df[[facility_col]]),
    LOS_days = los_df$LOS_days,
    stringsAsFactors = FALSE
  )

  # Class status per admission x pathogen x class (this pathogen's results only)
  status <- .los_class_status(
    data, admission_id_col, organism_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col, intermediate_as_resistant
  )
  # Covariates per admission x pathogen, from the earliest culture of the pathogen
  covs <- .los_admission_covariates(
    data, admission_id_col, organism_col, date_culture_col, infection_type_col,
    hai_value, cai_value, age_col, sex_col, unit_type_col, icu_values,
    comorbidity_col, cluster_id_col
  )
  cov_terms <- c(
    if (!is.null(infection_type_col)) "HAI",
    if (!is.null(age_col)) "Age_model",
    if (!is.null(sex_col)) "Sex_model",
    if (!is.null(unit_type_col)) "ICU",
    if (!is.null(comorbidity_col)) "Comorbidity_model"
  )

  pathogens <- if (!is.null(organism_name)) {
    unique(as.character(status$.org)[.name_key(status$.org) %in% .name_key(organism_name)])
  } else {
    sort(unique(as.character(status$.org)))
  }
  if (length(pathogens) == 0L) {
    warning("None of the requested pathogens were found.", call. = FALSE)
    return(.los_rr_empty())
  }

  alpha <- (1 - conf_level) / 2
  rr_rows <- list(); hosp_rows <- list(); draws <- list(); skipped <- list()
  skip <- function(path, cls, hosp, reason) {
    skipped[[length(skipped) + 1L]] <<- data.frame(
      pathogen = path, antibiotic_class = cls, hospital = hosp, reason = reason,
      stringsAsFactors = FALSE
    )
  }

  for (path in pathogens) {
    st_p <- status[status$.org == path, , drop = FALSE]
    path_classes <- .los_classes_for(classes, path, unique(st_p$.class))
    for (cls in path_classes) {
      st <- st_p[.name_key(st_p$.class) == .name_key(cls), c(".id", ".org", "R"), drop = FALSE]
      if (nrow(st) == 0L) {
        skip(path, cls, NA_character_, "class not tested for this pathogen")
        next
      }
      md <- merge(st, los_tbl, by = ".id")
      md <- merge(md, covs, by = c(".id", ".org"))
      n_before <- nrow(md)
      if (length(cov_terms)) md <- md[stats::complete.cases(md[, cov_terms, drop = FALSE]), , drop = FALSE]
      n_dropped <- n_before - nrow(md)

      groups <- if (hospital_effect == "stratified") sort(unique(md$.hosp)) else NA_character_
      for (grp in groups) {
        sub <- if (is.na(grp)) md else md[md$.hosp == grp, , drop = FALSE]
        n_r <- sum(sub$R == 1L); n_s <- sum(sub$R == 0L)
        if (nrow(sub) < min_n || n_r < min_resistant || n_s < min_susceptible) {
          skip(path, cls, grp, sprintf(
            "too few admissions (n = %d, R = %d, S = %d; minimum %d / %d / %d)",
            nrow(sub), n_r, n_s, min_n, min_resistant, min_susceptible
          ))
          next
        }
        res <- .los_fit_one(
          sub, cov_terms, use_hospital = hospital_effect == "fixed",
          reference_hospital = reference_hospital, ci_method = ci_method,
          n_boot = n_boot, alpha = alpha, return_draws = return_draws,
          flag_min = flag_min
        )
        if (!is.null(res$error)) {
          skip(path, cls, grp, res$error)
          next
        }
        rr_rows[[length(rr_rows) + 1L]] <- data.frame(
          pathogen = path, antibiotic_class = cls, hospital = grp,
          res$rr, n_dropped_covariates = n_dropped,
          hospital_effect = hospital_effect, ci_method = ci_method,
          model_family = "Gamma", link_function = "log",
          syndrome_scope = if (!is.null(syndrome_name)) syndrome_name else "all",
          stringsAsFactors = FALSE
        )
        hosp_rows[[length(hosp_rows) + 1L]] <- data.frame(
          pathogen = path, antibiotic_class = cls, res$hospital,
          stringsAsFactors = FALSE
        )
        if (return_draws && !is.null(res$draws)) {
          key <- paste(c(path, cls, if (!is.na(grp)) grp), collapse = " | ")
          draws[[key]] <- res$draws
        }
      }
    }
  }

  out <- list(
    rr = if (length(rr_rows)) do.call(rbind, rr_rows) else .los_rr_empty()$rr,
    hospital = if (length(hosp_rows)) do.call(rbind, hosp_rows) else .los_rr_empty()$hospital,
    draws = if (return_draws) draws else NULL,
    skipped = if (length(skipped)) do.call(rbind, skipped) else .los_rr_empty()$skipped
  )
  rownames(out$rr) <- NULL; rownames(out$hospital) <- NULL; rownames(out$skipped) <- NULL
  if (nrow(out$rr) == 0L) {
    warning("No Gamma GLM LOS models could be fitted.", call. = FALSE)
  } else {
    message(sprintf(
      "Gamma(log) LOS models fitted: %d pathogen-class estimate(s) (hospital effect: %s; CI: %s).",
      nrow(out$rr), hospital_effect, ci_method
    ))
  }
  out
}

# Empty return value of daly_fit_los_rr()
.los_rr_empty <- function() {
  list(
    rr = data.frame(pathogen = character(), antibiotic_class = character(),
                    hospital = character(), RR_LOS = numeric(), stringsAsFactors = FALSE),
    hospital = data.frame(pathogen = character(), antibiotic_class = character(),
                          hospital = character(), LOSR = numeric(), LOSS = numeric(),
                          stringsAsFactors = FALSE),
    draws = NULL,
    skipped = data.frame(pathogen = character(), antibiotic_class = character(),
                         hospital = character(), reason = character(), stringsAsFactors = FALSE)
  )
}

# Classes to model for one pathogen: NULL -> all tested; vector -> same for
# every pathogen; named list -> the pathogen's own panel.
.los_classes_for <- function(classes, path, available) {
  if (is.null(classes)) return(sort(available))
  if (is.list(classes)) {
    hit <- which(.name_key(names(classes)) == .name_key(path))
    if (!length(hit)) return(character())
    classes <- classes[[hit[1L]]]
  }
  as.character(classes)
}

# Worst phenotype per admission x pathogen x class: R = 1 if any result of
# that pathogen for any drug in the class is R. Intermediate: R (colistin: S)
# when intermediate_as_resistant, otherwise dropped.
.los_class_status <- function(data, id_col, org_col, class_col, name_col, value_col,
                              intermediate_as_resistant) {
  v <- toupper(trimws(as.character(data[[value_col]])))
  v[v %in% c("RESISTANT")] <- "R"
  v[v %in% c("SUSCEPTIBLE", "SENSITIVE")] <- "S"
  v[v %in% c("INTERMEDIATE")] <- "I"
  is_col <- grepl("colistin", tolower(as.character(data[[name_col]])), fixed = TRUE)
  if (intermediate_as_resistant) {
    v[v == "I"] <- ifelse(is_col[v == "I"], "S", "R")
  }
  keep <- v %in% c("R", "S") & !is.na(data[[class_col]]) &
    trimws(as.character(data[[class_col]])) != "" & !is.na(data[[id_col]])
  d <- data.frame(
    .id = as.character(data[[id_col]][keep]),
    .org = as.character(data[[org_col]][keep]),
    .class = as.character(data[[class_col]][keep]),
    .r = v[keep] == "R",
    stringsAsFactors = FALSE
  )
  out <- stats::aggregate(.r ~ .id + .org + .class, data = d, FUN = any)
  out$R <- as.integer(out$.r)
  out$.r <- NULL
  out
}

# Covariates per admission x pathogen, taken from the earliest culture of the
# pathogen in the admission.
.los_admission_covariates <- function(data, id_col, org_col, culture_col, inf_col,
                                      hai_value, cai_value, age_col, sex_col,
                                      unit_col, icu_values, comorb_col, cluster_col) {
  cul <- if (!is.null(culture_col) && culture_col %in% names(data)) {
    suppressWarnings(as.Date(data[[culture_col]]))
  } else {
    message(sprintf("Culture date column '%s' not found: covariates taken from the first row of each admission.", culture_col))
    rep(as.Date(NA), nrow(data))
  }
  ord <- order(as.character(data[[id_col]]), as.character(data[[org_col]]), cul, na.last = TRUE)
  d <- data[ord, , drop = FALSE]
  first <- !duplicated(paste(d[[id_col]], d[[org_col]], sep = "\r"))
  d <- d[first, , drop = FALSE]
  out <- data.frame(
    .id = as.character(d[[id_col]]),
    .org = as.character(d[[org_col]]),
    .cluster = if (!is.null(cluster_col)) as.character(d[[cluster_col]]) else as.character(d[[id_col]]),
    stringsAsFactors = FALSE
  )
  no_cl <- is.na(out$.cluster) | out$.cluster == ""
  out$.cluster[no_cl] <- out$.id[no_cl]
  if (!is.null(inf_col)) {
    it <- toupper(trimws(as.character(d[[inf_col]])))
    out$HAI <- ifelse(it == toupper(hai_value), 1L, ifelse(it == toupper(cai_value), 0L, NA_integer_))
  }
  if (!is.null(age_col)) out$Age_model <- suppressWarnings(as.numeric(d[[age_col]]))
  if (!is.null(sex_col)) {
    s <- trimws(as.character(d[[sex_col]]))
    s[s == "" | toupper(s) %in% c("NA", "UNKNOWN")] <- NA
    out$Sex_model <- s
  }
  if (!is.null(unit_col)) {
    u <- tolower(trimws(as.character(d[[unit_col]])))
    miss <- is.na(u) | u %in% c("", "na", "null")
    out$ICU <- ifelse(miss, NA_integer_,
                      as.integer(grepl(paste(tolower(icu_values), collapse = "|"), u)))
  }
  if (!is.null(comorb_col)) {
    enc <- suppressMessages(.encode_comorbidity_mortality(
      data.frame(.id = out$.id, d[comorb_col], stringsAsFactors = FALSE),
      comorbidity_col = comorb_col, patient_id_col = ".id"
    ))
    out$Comorbidity_model <- suppressWarnings(as.numeric(enc$comorbidity_encoded))
  }
  out
}

# Cluster-robust (sandwich) covariance for a Gamma log-link GLM. With the log
# link and Gamma variance the working weights are 1, so the bread is (X'X)^-1
# and each observation's score is x_i (y_i - mu_i) / mu_i.
.los_vcov_cluster <- function(X, y, mu, cluster) {
  r <- (y - mu) / mu
  bread <- solve(crossprod(X))
  meat <- rowsum(X * r, cluster)
  G <- nrow(meat)
  adj <- if (G > 1L) G / (G - 1L) else 1
  bread %*% crossprod(meat) %*% bread * adj
}

# LOSR / LOSS over `rows`: mean prediction with R = 1 and R = 0, and their
# gradients with respect to the coefficients (for the delta method).
.los_gcomp <- function(X, b, rows) {
  X1 <- X[rows, , drop = FALSE]; X1[, "R"] <- 1
  X0 <- X1; X0[, "R"] <- 0
  m1 <- exp(drop(X1 %*% b)); m0 <- exp(drop(X0 %*% b))
  list(LOSR = mean(m1), LOSS = mean(m0),
       g1 = colMeans(X1 * m1), g0 = colMeans(X0 * m0))
}

# Fit one Gamma(log) model and compute RR_LOS, LOSR/LOSS (all hospitals and
# per hospital) and their intervals.
.los_fit_one <- function(md, cov_terms, use_hospital, reference_hospital,
                         ci_method, n_boot, alpha, return_draws, flag_min = 10L) {
  terms <- "R"
  for (v in cov_terms) {
    if (length(unique(md[[v]])) > 1L) terms <- c(terms, v)
  }
  if ("Sex_model" %in% terms) {
    tab <- sort(table(md$Sex_model), decreasing = TRUE)
    md$Sex_model <- factor(md$Sex_model, levels = names(tab))
  }
  hosp_levels <- sort(unique(md$.hosp))
  ref <- NA_character_
  if (use_hospital && length(hosp_levels) > 1L) {
    counts <- table(md$.hosp)
    ref <- if (!is.null(reference_hospital) && reference_hospital %in% hosp_levels) {
      reference_hospital
    } else {
      names(counts)[which.max(counts)]
    }
    md$.hosp_f <- stats::relevel(factor(md$.hosp), ref = ref)
    terms <- c(terms, ".hosp_f")
  }
  form <- stats::as.formula(paste("LOS_days ~", paste(terms, collapse = " + ")))
  fit <- tryCatch(
    stats::glm(form, data = md, family = stats::Gamma(link = "log"),
               control = stats::glm.control(maxit = 100)),
    error = function(e) NULL
  )
  if (is.null(fit) || !fit$converged) return(list(error = "model did not converge"))
  b_all <- stats::coef(fit)
  if (is.na(b_all[["R"]])) return(list(error = "resistance term not estimable"))

  X <- stats::model.matrix(form, data = md)
  ok_col <- !is.na(b_all)
  X <- X[, ok_col, drop = FALSE]; b <- b_all[ok_col]
  y <- md$LOS_days; mu <- drop(exp(X %*% b))
  V <- .los_vcov_cluster(X, y, mu, md$.cluster)
  se_r <- sqrt(V["R", "R"])
  z <- stats::qnorm(1 - alpha)

  rx_p <- NA_real_
  if (use_hospital && length(hosp_levels) > 1L) {
    fit_i <- tryCatch(
      stats::glm(stats::update(form, . ~ . + R:.hosp_f), data = md,
                 family = stats::Gamma(link = "log"), control = stats::glm.control(maxit = 100)),
      error = function(e) NULL
    )
    if (!is.null(fit_i)) {
      rx_p <- tryCatch(stats::anova(fit, fit_i, test = "F")[["Pr(>F)"]][2L], error = function(e) NA_real_)
    }
  }

  is_r <- md$R == 1L
  hosp_rows <- function(h, rows_r, rows_h) {
    i <- rows_h & rows_r
    if (!any(i)) list(rows = rows_h, fallback = TRUE) else list(rows = i, fallback = FALSE)
  }
  net <- .los_gcomp(X, b, is_r)
  per_h <- lapply(hosp_levels, function(h) {
    sel <- hosp_rows(h, is_r, md$.hosp == h)
    c(.los_gcomp(X, b, sel$rows), fallback = sel$fallback)
  })
  names(per_h) <- hosp_levels

  # intervals
  ci <- function(est, lo, hi) c(est = est, lower = lo, upper = hi)
  na_ci <- function(est) ci(est, NA_real_, NA_real_)
  delta_ci <- function(g) {
    se1 <- sqrt(drop(g$g1 %*% V %*% g$g1)); se0 <- sqrt(drop(g$g0 %*% V %*% g$g0))
    seE <- sqrt(drop((g$g1 - g$g0) %*% V %*% (g$g1 - g$g0)))
    list(
      LOSR = ci(g$LOSR, exp(log(g$LOSR) - z * se1 / g$LOSR), exp(log(g$LOSR) + z * se1 / g$LOSR)),
      LOSS = ci(g$LOSS, exp(log(g$LOSS) - z * se0 / g$LOSS), exp(log(g$LOSS) + z * se0 / g$LOSS)),
      excess = ci(g$LOSR - g$LOSS, g$LOSR - g$LOSS - z * seE, g$LOSR - g$LOSS + z * seE)
    )
  }
  rr_est <- exp(b[["R"]])
  rr_ci <- na_ci(rr_est)
  net_ci <- list(LOSR = na_ci(net$LOSR), LOSS = na_ci(net$LOSS), excess = na_ci(net$LOSR - net$LOSS))
  h_ci <- lapply(per_h, function(g) list(LOSR = na_ci(g$LOSR), LOSS = na_ci(g$LOSS), excess = na_ci(g$LOSR - g$LOSS)))
  draws <- NULL
  n_boot_ok <- NA_integer_

  if (ci_method == "delta") {
    rr_ci <- ci(rr_est, exp(b[["R"]] - z * se_r), exp(b[["R"]] + z * se_r))
    net_ci <- delta_ci(net)
    h_ci <- lapply(per_h, delta_ci)
  } else if (ci_method == "bootstrap") {
    clusters <- md$.cluster
    cl_rows <- split(seq_len(nrow(md)), clusters)
    cl_hosp <- vapply(cl_rows, function(i) md$.hosp[i[1L]], character(1))
    by_h <- split(names(cl_rows), cl_hosp)
    b_rr <- rep(NA_real_, n_boot)
    b_net <- matrix(NA_real_, n_boot, 2L)
    b_r <- b_s <- matrix(NA_real_, n_boot, length(hosp_levels), dimnames = list(NULL, hosp_levels))
    start <- b
    for (k in seq_len(n_boot)) {
      ids <- unlist(lapply(by_h, function(p) p[sample.int(length(p), length(p), replace = TRUE)]), use.names = FALSE)
      rows <- unlist(cl_rows[ids], use.names = FALSE)
      bk <- tryCatch(
        stats::glm.fit(X[rows, , drop = FALSE], y[rows], family = stats::Gamma(link = "log"),
                       start = start, control = stats::glm.control(maxit = 100)),
        error = function(e) NULL
      )
      if (is.null(bk) || !bk$converged || is.na(bk$coefficients[["R"]])) next
      bb <- bk$coefficients; bb[is.na(bb)] <- 0
      Xk <- X[rows, , drop = FALSE]; rk <- md$R[rows] == 1L; hk <- md$.hosp[rows]
      g <- .los_gcomp(Xk, bb, rk)
      b_rr[k] <- exp(bb[["R"]]); b_net[k, ] <- c(g$LOSR, g$LOSS)
      for (h in hosp_levels) {
        sel <- hosp_rows(h, rk, hk == h)
        if (!any(sel$rows)) next
        gh <- .los_gcomp(Xk, bb, sel$rows)
        b_r[k, h] <- gh$LOSR; b_s[k, h] <- gh$LOSS
      }
    }
    ok <- !is.na(b_rr)
    n_boot_ok <- sum(ok)
    if (n_boot_ok < 2L) return(list(error = "bootstrap failed"))
    q <- function(x) unname(stats::quantile(x, c(alpha, 1 - alpha), na.rm = TRUE))
    rr_ci <- ci(rr_est, q(b_rr)[1L], q(b_rr)[2L])
    qci <- function(est, x) { qq <- q(x); ci(est, qq[1L], qq[2L]) }
    net_ci <- list(LOSR = qci(net$LOSR, b_net[, 1L]), LOSS = qci(net$LOSS, b_net[, 2L]),
                   excess = qci(net$LOSR - net$LOSS, b_net[, 1L] - b_net[, 2L]))
    h_ci <- lapply(hosp_levels, function(h) {
      g <- per_h[[h]]
      list(LOSR = qci(g$LOSR, b_r[, h]), LOSS = qci(g$LOSS, b_s[, h]),
           excess = qci(g$LOSR - g$LOSS, b_r[, h] - b_s[, h]))
    })
    names(h_ci) <- hosp_levels
    if (return_draws) {
      draws <- list(
        RR_LOS = b_rr[ok], LOSR = b_r[ok, , drop = FALSE], LOSS = b_s[ok, , drop = FALSE],
        LOSR_all = b_net[ok, 1L], LOSS_all = b_net[ok, 2L]
      )
    }
  }

  obs_mean <- function(x) if (length(x)) mean(x) else NA_real_
  rr <- data.frame(
    n_admissions = nrow(md), n_resistant = sum(is_r), n_susceptible = sum(!is_r),
    n_hospitals = length(hosp_levels),
    RR_LOS = rr_est, CI_lower = rr_ci[["lower"]], CI_upper = rr_ci[["upper"]],
    beta_log_rr = b[["R"]], se_log_rr_cluster = se_r,
    mean_los_resistant = obs_mean(y[is_r]), mean_los_susceptible = obs_mean(y[!is_r]),
    LOSR = net$LOSR, LOSR_lower = net_ci$LOSR[["lower"]], LOSR_upper = net_ci$LOSR[["upper"]],
    LOSS = net$LOSS, LOSS_lower = net_ci$LOSS[["lower"]], LOSS_upper = net_ci$LOSS[["upper"]],
    excess_LOS = net$LOSR - net$LOSS,
    excess_lower = net_ci$excess[["lower"]], excess_upper = net_ci$excess[["upper"]],
    dispersion = summary(fit)$dispersion, RxHospital_p = rx_p,
    covariates = paste(setdiff(terms, c("R", ".hosp_f")), collapse = ", "),
    reference_hospital = ref, n_boot_ok = n_boot_ok,
    stringsAsFactors = FALSE
  )
  hospital <- do.call(rbind, lapply(hosp_levels, function(h) {
    i <- md$.hosp == h; g <- per_h[[h]]; cc <- h_ci[[h]]
    n_r <- sum(i & is_r); n_s <- sum(i & !is_r)
    data.frame(
      hospital = h, n_resistant = n_r, n_susceptible = n_s,
      obs_mean_R = obs_mean(y[i & is_r]), obs_mean_S = obs_mean(y[i & !is_r]),
      pred_mean_R = obs_mean(mu[i & is_r]), pred_mean_S = obs_mean(mu[i & !is_r]),
      LOSR = g$LOSR, LOSR_lower = cc$LOSR[["lower"]], LOSR_upper = cc$LOSR[["upper"]],
      LOSS = g$LOSS, LOSS_lower = cc$LOSS[["lower"]], LOSS_upper = cc$LOSS[["upper"]],
      excess_LOS = g$LOSR - g$LOSS,
      excess_lower = cc$excess[["lower"]], excess_upper = cc$excess[["upper"]],
      RR_LOS = g$LOSR / g$LOSS,
      low_support = n_r < flag_min || n_s < flag_min,
      fallback_all_admissions = isTRUE(g$fallback),
      stringsAsFactors = FALSE
    )
  }))
  list(rr = rr, hospital = hospital, draws = draws)
}


# Step 2b
#' Estimate Per-Class LOS Relative Risk via Parametric Distribution Fitting
#'
#' For each pathogen x antibiotic class, splits patients into resistant (R) and
#' susceptible (S) groups, fits the candidate parametric distribution(s)
#' (default: gamma) to each group's LOS, selects the best fit by AIC, and takes
#' the \strong{theoretical median} of the selected distribution:
#'   RR_LOS = median_LOS(R) / median_LOS(S)
#'
#' Theoretical medians: Weibull (shape k, scale lambda)
#' \eqn{\lambda (\ln 2)^{1/k}}; log-normal (meanlog mu) \eqn{e^{\mu}};
#' gamma (shape alpha, rate beta) \code{qgamma(0.5, alpha, beta)} (no closed
#' form). A class is skipped when either group has fewer than \code{min_n}
#' patients or its distribution fit fails.
#'
#' This is an alternative to \code{daly_fit_los_rr()} that avoids regression
#' and models the full LOS distribution per group. The output flat data frame
#' is fully compatible with \code{assign_rr_to_profiles()} and
#' \code{filter_profiles_to_rr_classes()}.
#'
#' @param data Data frame of isolate/patient rows (long format).
#' @param patient_id_col Character. Patient identifier column.
#' @param facility_col Character. Facility/centre column.
#' @param organism_col Character. Organism name column.
#' @param syndrome_col Character. Syndrome column (used only when
#'   \code{syndrome_name} is not \code{NULL}).
#' @param infection_type_col Character. Column already holding the final
#'   HAI/CAI classification (e.g. \code{prep_derive_hai_cai()}'s output).
#'   Used as-is -- \strong{not} re-derived from admission/culture dates.
#'   Default \code{"infection_type"}.
#' @param antibiotic_class_col Character. Antibiotic class column.
#' @param antibiotic_name_col Character. Antibiotic name column.
#' @param antibiotic_value_col Character. Antibiotic value column (S/I/R).
#' @param date_admission_col Character. Date of admission column.
#' @param date_discharge_col Character. Date of discharge column.
#' @param date_culture_col Character. Unused -- retained only for backward API
#'   compatibility. LOS is computed as final outcome date - admission date and
#'   does not depend on the culture date.
#' @param final_outcome_col Character. Final outcome column.
#' @param final_outcome_value Character. Value indicating discharge.
#' @param syndrome_name Character or \code{NULL}. Restrict to this syndrome.
#' @param organism_name Character vector or \code{NULL}. Restrict to these
#'   pathogen(s); otherwise all pathogens in the filtered data.
#' @param facility_name Character or \code{NULL}. If provided, filters data
#'   to the specified facility before fitting. Default \code{NULL}.
#' @param hai_threshold_hours Numeric. Unused -- retained only for backward
#'   API compatibility. HAI/CAI is read directly from \code{infection_type_col}
#'   rather than re-derived from a date gap.
#' @param distributions Character vector. Candidate distributions to fit.
#'   Default \code{"gamma"}.
#' @param max_los Numeric. Maximum plausible LOS in days. Default \code{200}.
#' @param min_n Integer. Minimum patients required in both R and S groups.
#'   Default \code{10L}.
#'
#' @return A data frame with one row per pathogen x class:
#'   \describe{
#'     \item{pathogen}{Pathogen name.}
#'     \item{antibiotic_class}{Antibiotic class name.}
#'     \item{RR_LOS}{Ratio of fitted medians: median_LOS(R) / median_LOS(S).}
#'     \item{median_LOS_R}{Fitted median LOS (days) for resistant patients.}
#'     \item{median_LOS_S}{Fitted median LOS (days) for susceptible patients.}
#'     \item{median_LOS_R_years, median_LOS_S_years}{The same medians in years
#'       (days / 365).}
#'     \item{best_dist_R}{Distribution fitted to the R group.}
#'     \item{best_dist_S}{Distribution fitted to the S group.}
#'     \item{n_R}{Number of resistant patients.}
#'     \item{n_S}{Number of susceptible patients.}
#'     \item{syndrome_scope}{Syndrome filter applied, or "all".}
#'   }
#'
#' @section Two options for profile-level LOS:
#' This function is \strong{option 2} (per-class with max rule): its output is
#' passed to \code{\link{daly_assign_rr_to_profiles}()} (with
#' \code{los_r_col = "median_LOS_R"}, \code{los_s_col = "median_LOS_S"}), which
#' gives each resistance profile the estimate of its dominant resistant class.
#' The S group is "susceptible to that class", not fully susceptible.
#' \strong{Option 1} is \code{\link{daly_fit_los_distn_by_profile}()}, which follows
#' Section 5.3 directly by fitting LOS per patient-level resistance profile
#' against the fully susceptible profile; use it when patient-level profiles
#' are available and adequately populated.
#'
#' @seealso \code{\link{daly_fit_los_rr}}, \code{\link{daly_assign_rr_to_profiles}},
#'   \code{\link{daly_fit_los_distn_by_profile}} (option 1)
#' @export
daly_fit_los_rr_distribution <- function(
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
) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required: install.packages('fitdistrplus')")
  }

  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col,
    final_outcome_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # -- HAI / CAI classification is already computed upstream (e.g. by
  # prep_derive_hai_cai()) and available in infection_type_col -- read it
  # directly rather than re-deriving it from the admission/culture date gap.
  data$infection_type_derived <- data[[infection_type_col]]

  # -- Filter to discharged patients + optional syndrome + optional facility --
  df <- dplyr::filter(data, .data[[final_outcome_col]] == final_outcome_value)
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .name_key(.data[[syndrome_col]]) == .name_key(syndrome_name))
  }
  if (!is.null(facility_name)) {
    df <- dplyr::filter(df, .data[[facility_col]] == facility_name)
  }

  if (nrow(df) == 0L) {
    warning(sprintf(
      "No discharged patients remain after filters (outcome='%s'%s%s). Returning empty data frame.",
      final_outcome_value,
      if (!is.null(syndrome_name)) sprintf(", syndrome='%s'", syndrome_name) else "",
      if (!is.null(facility_name)) sprintf(", facility='%s'", facility_name) else ""
    ))
    return(data.frame())
  }

  # Resolve pathogens
  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  # Inner helper: fit candidate distribution(s), select the best by AIC, and
  # return the theoretical median of the selected distribution
  .fit_dist_median <- function(los_vec) {
    los_vec <- los_vec[is.finite(los_vec) & los_vec > 0]
    if (length(los_vec) == 0L) {
      return(list(median = NA_real_, dist = NA_character_))
    }
    fits <- Filter(
      Negate(is.null),
      setNames(
        lapply(distributions, .safe_fitdist, x = los_vec),
        distributions
      )
    )
    if (length(fits) == 0L) {
      return(list(median = NA_real_, dist = NA_character_))
    }
    aics <- sapply(fits, `[[`, "aic")
    best <- names(which.min(aics))
    list(median = .compute_median_from_fit(fits[[best]], best), dist = best)
  }

  all_rr <- list()

  for (path in pathogens) {
    path_df <- df %>%
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
          stringr::str_to_lower(stringr::str_trim(path))
      )
    if (nrow(path_df) == 0L) {
      message(sprintf("'%s': no data after filters, skipping.", path))
      next
    }

    # Patient LOS (final outcome date - admission date)
    los_pat <- daly_compute_patient_los(
      path_df,
      patient_id_col      = patient_id_col,
      facility_col        = facility_col,
      facility_name       = NULL,
      date_admission_col  = date_admission_col,
      date_discharge_col  = date_discharge_col,
      final_outcome_col   = final_outcome_col,
      final_outcome_value = final_outcome_value,
      max_los             = max_los
    )
    if (nrow(los_pat) == 0L) {
      message(sprintf("'%s': no patients with valid LOS, skipping.", path))
      next
    }

    # Class resistance wide matrix (patient x class, binary 0/1)
    # untested_fill = NA_integer_: patients not tested for a class get NA,
    # so they are excluded from both n_R and n_S for that class (not treated
    # as susceptible). Only patients actually tested contribute to the fit.
    resist_wide <- .build_class_resistance_wide(
      path_df,
      patient_id_col       = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col  = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col,
      untested_fill        = NA_integer_
    )
    class_name_map <- attr(resist_wide, "class_name_map")
    class_safe <- setdiff(names(resist_wide), patient_id_col)

    model_data <- dplyr::inner_join(los_pat, resist_wide, by = patient_id_col)

    if (nrow(model_data) < min_n) {
      message(sprintf(
        "'%s': only %d patients after merging, skipping.",
        path, nrow(model_data)
      ))
      next
    }

    n_fitted <- 0L

    # Per-class R vs S distribution fit
    for (cls in class_safe) {
      orig_class <- class_name_map[[cls]]

      los_r <- model_data$LOS_days[
        !is.na(model_data[[cls]]) & model_data[[cls]] == 1L
      ]
      los_s <- model_data$LOS_days[
        !is.na(model_data[[cls]]) & model_data[[cls]] == 0L
      ]

      n_r <- length(los_r)
      n_s <- length(los_s)

      if (n_r < min_n || n_s < min_n) {
        message(sprintf(
          "'%s' | class '%s': n_R=%d / n_S=%d -- one group below min_n=%d, skipping.",
          path, orig_class, n_r, n_s, min_n
        ))
        next
      }

      est_r <- .fit_dist_median(los_r)
      est_s <- .fit_dist_median(los_s)

      if (is.na(est_r$median) || is.na(est_s$median) || est_s$median <= 0) {
        message(sprintf(
          "'%s' | class '%s': distribution fitting failed (median_R=%s, median_S=%s), skipping.",
          path, orig_class,
          if (is.na(est_r$median)) "NA" else round(est_r$median, 3L),
          if (is.na(est_s$median)) "NA" else round(est_s$median, 3L)
        ))
        next
      }

      # Relative LOS = ratio of fitted medians (resistant / susceptible)
      rr_los <- est_r$median / est_s$median

      if (rr_los > 10) {
        message(sprintf(
          "'%s' | class '%s': extreme RR_LOS = %.2f (median_R=%.2f / median_S=%.2f) -- check data.",
          path, orig_class, rr_los, est_r$median, est_s$median
        ))
      } else if (rr_los < 1) {
        message(sprintf(
          "'%s' | class '%s': RR_LOS = %.4f < 1 (resistant patients have shorter LOS than susceptible) -- verify.",
          path, orig_class, rr_los
        ))
      }

      all_rr[[length(all_rr) + 1L]] <- data.frame(
        pathogen         = path,
        antibiotic_class = orig_class,
        facility_name    = if (!is.null(facility_name)) facility_name else NA_character_,
        RR_LOS           = round(rr_los, 4L),
        median_LOS_R       = round(est_r$median, 3L),
        median_LOS_S       = round(est_s$median, 3L),
        median_LOS_R_years = round(est_r$median / 365, 6L),
        median_LOS_S_years = round(est_s$median / 365, 6L),
        best_dist_R      = est_r$dist,
        best_dist_S      = est_s$dist,
        n_R              = n_r,
        n_S              = n_s,
        syndrome_scope   = if (!is.null(syndrome_name)) syndrome_name else "all",
        stringsAsFactors = FALSE
      )
      n_fitted <- n_fitted + 1L
    }

    message(sprintf(
      "'%s': %d class relative LOS(s) estimated (fitted median LOS, syndrome=%s).",
      path, n_fitted,
      if (!is.null(syndrome_name)) syndrome_name else "all"
    ))
  }

  if (length(all_rr) == 0L) {
    warning("No relative LOS values could be estimated. Check data, filters, and min_n.")
    return(data.frame(
      pathogen = character(), antibiotic_class = character(),
      facility_name = character(), RR_LOS = numeric(),
      median_LOS_R = numeric(), median_LOS_S = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  return(dplyr::bind_rows(all_rr))
}


# Step 2c
#' Estimate Profile-Specific Median LOS (Section 5.3)
#'
#' For each resistance profile \eqn{\delta} of each pathogen \eqn{k}, fits
#' parametric distributions (default: Weibull, log-normal, gamma) to the LOS
#' of the patients with that profile, \eqn{\{LOS_i : \delta_i = \delta\}},
#' selects the distribution with the lowest AIC, and takes the
#' \strong{theoretical median} of the selected distribution:
#' \deqn{LOS_{k\delta} = Median_{f_\delta}(LOS)}
#' The same procedure applied to the fully susceptible reference profile
#' \eqn{\delta_S} gives \eqn{LOS_{k\delta_S}}, and
#' \code{RR_LOS = LOS_{k delta} / LOS_{k delta_S}}.
#'
#' \strong{Profiles are not derived here.} Each patient's resistance profile
#' must be supplied in \code{profile_data} (e.g. the \code{profile_delta}
#' labels from the \code{profile_*} pipeline). Patients without a profile, or
#' with conflicting profiles, are excluded and reported in a message.
#'
#' LOS is computed by \code{\link{daly_compute_patient_los}()}
#' (final outcome date - admission date, discharged patients only,
#' 0 < LOS <= \code{max_los}). Theoretical medians: Weibull
#' \eqn{\lambda (\ln 2)^{1/k}}; log-normal \eqn{e^{\mu}}; gamma
#' \code{qgamma(0.5, shape, rate)}.
#'
#' @param data Data frame used for LOS (patient, facility, organism, dates,
#'   final outcome).
#' @param profile_data Data frame with one profile per patient (or per
#'   patient x pathogen): columns \code{profile_id_col}, \code{profile_col},
#'   and optionally \code{profile_pathogen_col}.
#' @param patient_id_col Character. Patient identifier in \code{data}.
#' @param facility_col Character. Facility / centre column in \code{data}.
#' @param organism_col Character. Organism column in \code{data}; used to
#'   assign pathogens when \code{profile_pathogen_col} is \code{NULL}.
#' @param syndrome_col Character. Syndrome column (used only with
#'   \code{syndrome_name}).
#' @param date_admission_col Character. Admission date column.
#' @param date_discharge_col Character. Final outcome / discharge date column.
#' @param final_outcome_col Character. Final outcome column.
#' @param final_outcome_value Character. Value indicating discharge.
#' @param profile_id_col Character. ID column in \code{profile_data} matching
#'   \code{patient_id_col}. Default \code{patient_id_col}.
#' @param profile_col Character. Profile label column in \code{profile_data}
#'   (e.g. \code{"RRSS"}). Default \code{"profile_delta"}.
#' @param profile_pathogen_col Character or \code{NULL}. Pathogen column in
#'   \code{profile_data}, when profiles are pathogen-specific. If \code{NULL},
#'   pathogens are taken from \code{organism_col} in \code{data}.
#' @param ref_profile Character or \code{NULL}. Label of the fully susceptible
#'   reference profile. If \code{NULL}, the label consisting only of
#'   \code{"S"} characters (e.g. \code{"SSSS"}) is used per pathogen.
#' @param syndrome_name Character or \code{NULL}. Restrict to this syndrome.
#' @param organism_name Character vector or \code{NULL}. Restrict to these
#'   pathogens.
#' @param facility_name Character or \code{NULL}. Restrict to this facility.
#' @param distributions Character vector of candidate distributions.
#'   Default \code{c("weibull", "lnorm", "gamma")}.
#' @param max_los Numeric. Maximum LOS in days. Default \code{365}.
#' @param min_n Integer. Minimum LOS observations required to fit a profile.
#'   Default \code{10L}.
#'
#' @return Data frame with one row per pathogen x profile: \code{pathogen},
#'   \code{profile_delta}, \code{is_reference}, \code{n_patients}, \code{n_los}
#'   (LOS observations), \code{median_LOS} (days), \code{median_LOS_years},
#'   \code{best_dist}, \code{aic}, \code{median_LOS_ref}, \code{RR_LOS},
#'   \code{status} (\code{"fitted"}, \code{"insufficient_n"},
#'   \code{"fit_failed"}) and \code{syndrome_scope}. Profiles that cannot be
#'   fitted are kept with \code{NA} estimates.
#'
#' @section Two options for profile-level LOS:
#' The package offers two ways to obtain LOS for each resistance profile:
#' \enumerate{
#'   \item \strong{Profile-specific (this function; Section 5.3).} Fits the
#'     LOS of the patients who actually have each profile \eqn{\delta} and
#'     compares it with the fully susceptible profile \eqn{\delta_S}.
#'     Requires patient-level profiles (e.g. from the \code{profile_*}
#'     pipeline) and enough patients per profile (\code{min_n}).
#'   \item \strong{Per-class with max rule.}
#'     \code{\link{daly_fit_los_rr_distribution}()} fits resistant vs
#'     susceptible LOS separately for each antibiotic class, then
#'     \code{\link{daly_assign_rr_to_profiles}()} gives each profile the
#'     estimate of its dominant (highest-RR) resistant class. Needs no
#'     patient-level profiles and pools data across profiles, so it is a
#'     fallback or sensitivity analysis when profiles are sparse; its
#'     reference group is "susceptible to that class", not fully susceptible.
#' }
#' Note: the YLD functions (\code{\link{daly_calc_yld_associated}()},
#' \code{\link{daly_calc_yld_attributable}()}) currently read the
#' \code{LOSR_years}/\code{LOSS_years} columns produced by option 2; the output
#' of this function is not yet joined into that format.
#'
#' @seealso \code{\link{daly_fit_los_rr_distribution}} and
#'   \code{\link{daly_assign_rr_to_profiles}} (option 2)
#' @export
daly_fit_los_distn_by_profile <- function(
  data,
  profile_data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  profile_id_col = patient_id_col,
  profile_col = "profile_delta",
  profile_pathogen_col = NULL,
  ref_profile = NULL,
  syndrome_name = NULL,
  organism_name = NULL,
  facility_name = NULL,
  distributions = c("weibull", "lnorm", "gamma"),
  max_los = 365,
  min_n = 10L
) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required: install.packages('fitdistrplus')")
  }

  # Column validation
  prof_required <- c(profile_id_col, profile_col, profile_pathogen_col)
  missing_prof <- setdiff(prof_required, names(profile_data))
  if (length(missing_prof) > 0L) {
    stop(sprintf(
      "Column(s) not found in profile_data: %s",
      paste(missing_prof, collapse = ", ")
    ))
  }
  if (is.null(profile_pathogen_col) && !organism_col %in% names(data)) {
    stop(sprintf(
      "organism_col '%s' not found in data (needed when profile_pathogen_col is NULL).",
      organism_col
    ))
  }

  # Step 1: LOS (final outcome date - admission date)
  los <- daly_compute_patient_los(
    data,
    patient_id_col = patient_id_col,
    facility_col = facility_col,
    facility_name = facility_name,
    syndrome_col = syndrome_col,
    syndrome_name = syndrome_name,
    date_admission_col = date_admission_col,
    date_discharge_col = date_discharge_col,
    final_outcome_col = final_outcome_col,
    final_outcome_value = final_outcome_value,
    max_los = max_los
  )
  los <- data.frame(
    .pid = as.character(los[[patient_id_col]]),
    LOS_days = los$LOS_days,
    stringsAsFactors = FALSE
  )

  # Step 2: supplied profiles, one per patient (x pathogen)
  prof <- data.frame(
    .pid = as.character(profile_data[[profile_id_col]]),
    .path = if (is.null(profile_pathogen_col)) NA_character_ else
      as.character(profile_data[[profile_pathogen_col]]),
    .profile = trimws(as.character(profile_data[[profile_col]])),
    stringsAsFactors = FALSE
  )
  prof <- unique(prof[!is.na(prof$.pid) & !is.na(prof$.profile) & prof$.profile != "", ])
  prof_key <- paste(prof$.pid, prof$.path, sep = "||")
  conflict <- prof_key %in% prof_key[duplicated(prof_key)]
  if (any(conflict)) {
    message(sprintf(
      "Excluding %d patient(s) with more than one profile in profile_data.",
      length(unique(prof_key[conflict]))
    ))
    prof <- prof[!conflict, ]
  }

  # Pathogen assignment
  if (is.null(profile_pathogen_col)) {
    orgs <- unique(data.frame(
      .pid = as.character(data[[patient_id_col]]),
      .org = as.character(data[[organism_col]]),
      stringsAsFactors = FALSE
    ))
    orgs <- orgs[!is.na(orgs$.org) & trimws(orgs$.org) != "", ]
    prof$.path <- NULL
    prof <- merge(prof, orgs, by = ".pid")
    names(prof)[names(prof) == ".org"] <- ".path"
  }

  n_los_pat <- length(unique(los$.pid))
  model_data <- merge(los, prof, by = ".pid")
  n_unmatched <- n_los_pat - length(unique(model_data$.pid))
  if (n_unmatched > 0L) {
    message(sprintf(
      "%d patient(s) with valid LOS have no profile in profile_data and are excluded.",
      n_unmatched
    ))
  }
  if (nrow(model_data) == 0L) {
    warning("No patients with both a valid LOS and a supplied profile.")
    return(data.frame())
  }

  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(model_data$.path))
  }

  results <- list()
  for (path in pathogens) {
    pdat <- model_data[.name_key(model_data$.path) == .name_key(path), ]
    if (nrow(pdat) == 0L) {
      message(sprintf("'%s': no patients with LOS and profile, skipping.", path))
      next
    }

    profiles <- sort(unique(pdat$.profile))
    ref <- if (!is.null(ref_profile)) {
      ref_profile
    } else {
      cand <- profiles[grepl("^S+$", profiles)]
      if (length(cand) == 1L) cand else NA_character_
    }
    if (is.na(ref) || !ref %in% profiles) {
      message(sprintf(
        "'%s': fully susceptible reference profile not found; RR_LOS will be NA.",
        path
      ))
    }

    rows <- lapply(profiles, function(pr) {
      los_vec <- pdat$LOS_days[pdat$.profile == pr]
      n_los <- length(los_vec)
      est <- if (n_los >= min_n) {
        .fit_best_dist_median(los_vec, distributions)
      } else {
        list(median = NA_real_, dist = NA_character_, aic = NA_real_)
      }
      status <- if (n_los < min_n) {
        "insufficient_n"
      } else if (is.na(est$median)) {
        "fit_failed"
      } else {
        "fitted"
      }
      data.frame(
        pathogen = path,
        profile_delta = pr,
        is_reference = identical(pr, ref),
        n_patients = length(unique(pdat$.pid[pdat$.profile == pr])),
        n_los = n_los,
        median_LOS = round(est$median, 3L),
        median_LOS_years = round(est$median / 365, 6L),
        best_dist = est$dist,
        aic = est$aic,
        status = status,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)

    ref_median <- if (!is.na(ref) && ref %in% out$profile_delta) {
      out$median_LOS[out$profile_delta == ref]
    } else {
      NA_real_
    }
    out$median_LOS_ref <- ref_median
    out$RR_LOS <- if (is.na(ref_median) || ref_median <= 0) {
      NA_real_
    } else {
      round(out$median_LOS / ref_median, 4L)
    }
    out$syndrome_scope <- if (!is.null(syndrome_name)) syndrome_name else "all"
    results[[length(results) + 1L]] <- out

    message(sprintf(
      "'%s': %d of %d profile(s) fitted (reference '%s').",
      path, sum(out$status == "fitted"), nrow(out), ref
    ))
  }

  if (length(results) == 0L) {
    warning("No profile-specific LOS could be estimated.")
    return(data.frame())
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out[, c(
    "pathogen", "profile_delta", "is_reference", "n_patients", "n_los",
    "median_LOS", "median_LOS_years", "best_dist", "aic",
    "median_LOS_ref", "RR_LOS", "status", "syndrome_scope"
  )]
}


# HAI/CAI classification for mortality cohort
#' Classify HAI/CAI Infection Type for the Mortality Cohort
#'
#' Derives the healthcare-associated infection (HAI) vs community-acquired
#' infection (CAI) classification for each patient record in the mortality
#' cohort, using admission-to-culture date differences. Also performs a
#' data-quality check for deceased patients missing an outcome date.
#'
#' @param data                 Data frame of patient-level records.
#' @param infection_type_col   Character. Column containing the raw infection
#'   type label. Default \code{"type_of_infection"}.
#' @param date_admission_col   Character. Admission date column.
#'   Default \code{"date_of_admission"}.
#' @param date_culture_col     Character. First positive culture date column.
#'   Default \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col    Character. Final outcome column.
#'   Default \code{"final_outcome"}.
#' @param final_outcome_date_col Character. Final outcome date column.
#'   Default \code{"final_outcome_date"}.
#' @param death_value          Character. Value in \code{final_outcome_col}
#'   indicating a fatal outcome. Default \code{"Death"}.
#' @param hai_threshold_hours  Numeric. Hours after admission at which a
#'   positive culture is classified as HAI. Default \code{48}.
#' @param patient_id_col       Character. Patient identifier column.
#'   Default \code{"PatientInformation_id"}.
#'
#' @return The input data frame with a derived \code{infection_type} column
#'   (values \code{"HAI"} or \code{"CAI"}).
#' @export

daly_derive_hai_cai_for_mortality <- function(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  death_value = "Death",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
) {
  # Column validation
  required <- c(
    infection_type_col, date_admission_col,
    date_culture_col, final_outcome_col
  )
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  has_outcome_date <- final_outcome_date_col %in% names(data)

  # Death-date data-quality check
  # Flag deceased patients who are missing their outcome date but do have
  # admission or culture date recorded.  HAI/CAI derivation proceeds normally
  # (it only needs admission vs culture date); this is purely a data-quality
  # report so the analyst can trace incomplete records.
  if (has_outcome_date) {
    is_dead <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[final_outcome_col]])
    )) == stringr::str_to_upper(death_value)

    missing_outcome_dt <- is.na(suppressWarnings(
      as.Date(as.character(data[[final_outcome_date_col]]))
    ))
    has_admit_or_culture <- !is.na(data[[date_admission_col]]) |
      !is.na(data[[date_culture_col]])

    dead_no_date <- is_dead & missing_outcome_dt & has_admit_or_culture

    if (any(dead_no_date)) {
      n_flag <- sum(dead_no_date)
      message(sprintf(
        paste0(
          "[Mortality] Data-quality flag: %d patient row(s) have ",
          "final_outcome='%s' but no outcome date recorded ",
          "(admission/culture dates are present). ",
          "HAI/CAI derivation is unaffected. ",
          "Review these records for date completeness."
        ),
        n_flag, death_value
      ))
      if (patient_id_col %in% names(data)) {
        flagged_ids <- unique(
          data[dead_no_date, patient_id_col, drop = TRUE]
        )
        message(sprintf(
          "  Flagged patient IDs (first 20 of %d): %s",
          length(flagged_ids),
          paste(head(flagged_ids, 20L), collapse = ", ")
        ))
      }
      flagged_df <- data[dead_no_date, , drop = FALSE]
    } else {
      flagged_df <- data[0L, , drop = FALSE]
    }
  } else {
    message(sprintf(
      "[Mortality] Column '%s' not found; skipping death-date quality check.",
      final_outcome_date_col
    ))
    flagged_df <- data[0L, , drop = FALSE]
  }

  # Identify ambiguous labels
  ambiguous_mask <- {
    inf_up <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[infection_type_col]])
    ))
    inf_up %in% c("NOT KNOWN", "NOT_KNOWN", "UNKNOWN", "NULL", "NA", "") |
      is.na(data[[infection_type_col]])
  }
  missing_admit <- is.na(data[[date_admission_col]])
  missing_culture <- is.na(data[[date_culture_col]])
  cannot_infer <- ambiguous_mask & (missing_admit | missing_culture)

  if (any(cannot_infer)) {
    n_cannot <- sum(cannot_infer)
    message(sprintf(
      paste0(
        "[Mortality] Cannot infer HAI/CAI for %d patient row(s): ",
        "label is ambiguous AND >=1 date is missing. ",
        "Assigned 'Not Known' -- these rows are excluded from the ",
        "mortality regression (missing HAI covariate)."
      ),
      n_cannot
    ))
    if (patient_id_col %in% names(data)) {
      ids_flag <- unique(data[cannot_infer, patient_id_col, drop = TRUE])
      message(sprintf(
        "  Patient IDs (first 20 of %d): %s",
        length(ids_flag),
        paste(head(ids_flag, 20L), collapse = ", ")
      ))
    }
    message(sprintf(
      "  Breakdown -- missing %s: %d | missing %s: %d",
      date_admission_col, sum(ambiguous_mask & missing_admit),
      date_culture_col,   sum(ambiguous_mask & missing_culture)
    ))
  }

  # Derive infection type
  # Explicit labels take priority; date-gap only used when label is ambiguous.
  data <- data %>%
    dplyr::mutate(
      .inf_raw = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[infection_type_col]])
      )),
      .gap_h = as.numeric(difftime(
        as.Date(suppressWarnings(
          as.character(.data[[date_culture_col]])
        )),
        as.Date(suppressWarnings(
          as.character(.data[[date_admission_col]])
        )),
        units = "hours"
      )),
      .cannot_infer = (.inf_raw %in% c(
        "NOT KNOWN", "NOT_KNOWN",
        "UNKNOWN", "NULL", "NA", ""
      ) |
        is.na(.data[[infection_type_col]])),
      infection_type_derived = dplyr::case_when(
        # Explicit HAI labels
        .inf_raw %in% c(
          "HAI",
          "HOSPITAL ACQUIRED",
          "HOSPITAL-ACQUIRED",
          "HOSPITAL_ACQUIRED",
          "NOSOCOMIAL",
          "HOSPITAL ACQUIRED INFECTION",
          "HOSPITAL-ACQUIRED INFECTION"
        ) ~ "HAI",
        # Explicit CAI labels
        .inf_raw %in% c(
          "CAI",
          "COMMUNITY ACQUIRED",
          "COMMUNITY-ACQUIRED",
          "COMMUNITY_ACQUIRED",
          "COMMUNITY ACQUIRED INFECTION",
          "COMMUNITY-ACQUIRED INFECTION"
        ) ~ "CAI",
        # Date-gap derivation (both dates present)
        .cannot_infer &
          !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours,
            "CAI", "HAI"
          ),
        # Cannot infer
        .cannot_infer ~ "Not Known",
        # Unrecognised non-empty label
        TRUE ~ "Not Known"
      )
    ) %>%
    dplyr::select(-".inf_raw", -".gap_h", -".cannot_infer")

  n_hai <- sum(data$infection_type_derived == "HAI", na.rm = TRUE)
  n_cai <- sum(data$infection_type_derived == "CAI", na.rm = TRUE)
  n_not_known <- sum(data$infection_type_derived == "Not Known", na.rm = TRUE)
  message(sprintf(
    "[Mortality] Infection type: HAI = %d | CAI = %d | Not Known = %d (excluded from model).",
    n_hai, n_cai, n_not_known
  ))

  attr(data, "missing_death_date_patients") <- flagged_df
  return(data)
}


# Step M2
#' Derive ICU Binary Flag per Patient
#'
#' Collapses unit-type data (one row per drug test per patient) to a
#' patient-level binary ICU indicator using the \strong{"ever in ICU"} rule:
#' if a patient has \emph{any} row with an ICU-type location during that
#' admission, \code{ICU = 1}; otherwise \code{ICU = 0}.
#'
#' Missing unit-type values are coded as \code{NA}. When the proportion of
#' patients with \emph{all} unit-type rows missing exceeds
#' \code{missing_threshold}, the function returns a three-level ordered factor
#' (\code{"ICU"} / \code{"Ward"} / \code{"Unknown"}) and emits a message
#' recommending that the \code{"Unknown"} level be included in the model.
#'
#' @param data Data frame at drug-test level (before patient collapse).
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param unit_type_col Character. Column containing ward/ICU labels.
#'   Default \code{"unit_type"}.
#' @param icu_values Character vector. Values (matched case-insensitively) that
#'   indicate ICU admission.
#'   Default \code{c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU")}.
#' @param missing_threshold Numeric in [0, 1]. Proportion of patients with
#'   entirely missing unit-type above which a 3-level factor is returned.
#'   Default \code{0.10}.
#'
#' @return Patient-level data frame: \code{patient_id_col} and \code{ICU}
#'   (integer 0/1, or factor \code{"ICU"/"Ward"/"Unknown"} when missing is
#'   high).
#' @keywords internal
.derive_icu_binary <- function(
  data,
  patient_id_col = "PatientInformation_id",
  unit_type_col = "unit_type",
  icu_values = c(
    "ICU", "Intensive Care", "Critical Care",
    "PICU", "NICU"
  ),
  missing_threshold = 0.10
) {
  if (!unit_type_col %in% names(data)) {
    message(sprintf(
      "[ICU] Column '%s' not found -- ICU covariate will be omitted.",
      unit_type_col
    ))
    return(
      dplyr::distinct(data, !!rlang::sym(patient_id_col)) %>%
        dplyr::mutate(ICU = NA_integer_)
    )
  }

  icu_pattern <- paste(stringr::str_to_lower(icu_values), collapse = "|")

  per_patient <- data %>%
    dplyr::mutate(
      .unit_norm = stringr::str_to_lower(stringr::str_trim(
        as.character(.data[[unit_type_col]])
      )),
      .is_missing = is.na(.data[[unit_type_col]]) |
        .unit_norm %in% c("na", "", "null"),
      .is_icu = !.is_missing &
        grepl(icu_pattern, .unit_norm, fixed = FALSE)
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::summarise(
      .any_icu     = any(.is_icu, na.rm = TRUE),
      .all_missing = all(.is_missing),
      .groups      = "drop"
    )

  n_patients <- nrow(per_patient)
  prop_missing <- sum(per_patient$.all_missing) / n_patients

  if (prop_missing > missing_threshold) {
    message(sprintf(
      paste0(
        "[ICU] %.1f%% of patients have entirely missing unit_type. ",
        "Returning 3-level factor (ICU / Ward / Unknown). ",
        "Consider including 'Unknown' as a separate level in the model."
      ),
      prop_missing * 100
    ))
    per_patient <- per_patient %>%
      dplyr::mutate(
        ICU = dplyr::case_when(
          .all_missing ~ "Unknown",
          .any_icu ~ "ICU",
          TRUE ~ "Ward"
        ),
        ICU = factor(ICU, levels = c("Ward", "ICU", "Unknown"))
      )
    n_icu <- sum(per_patient$ICU == "ICU", na.rm = TRUE)
    n_ward <- sum(per_patient$ICU == "Ward", na.rm = TRUE)
    n_unk <- sum(per_patient$ICU == "Unknown", na.rm = TRUE)
  } else {
    per_patient <- per_patient %>%
      dplyr::mutate(
        ICU = dplyr::case_when(
          .all_missing ~ NA_integer_,
          .any_icu ~ 1L,
          TRUE ~ 0L
        )
      )
    n_icu <- sum(per_patient$ICU == 1L, na.rm = TRUE)
    n_ward <- sum(per_patient$ICU == 0L, na.rm = TRUE)
    n_unk <- sum(is.na(per_patient$ICU))
  }

  message(sprintf(
    "[ICU] Per-patient ICU flag: ICU = %d | Ward = %d | Unknown/NA = %d.",
    n_icu, n_ward, n_unk
  ))

  per_patient %>%
    dplyr::select(!!rlang::sym(patient_id_col), ICU)
}


# Step M3
#' Encode Comorbidity Column for Mortality Model
#'
#' Standardises a free-text or numeric comorbidity column to a consistent
#' coding for use as a covariate in \code{daly_fit_mortality_or()}.
#'
#' Three encoding strategies are applied in order:
#' \enumerate{
#'   \item \strong{Numeric} (Charlson / Elixhauser index already present):
#'         used as-is.
#'   \item \strong{Binary text} (\code{"present"} / \code{"none"}, etc.):
#'         recoded to integer 0 / 1.
#'   \item \strong{Ordinal text} (\code{"none"} / \code{"mild"} /
#'         \code{"moderate"} / \code{"severe"}): recoded to an ordered factor.
#' }
#' Missing / unknown values are set to \code{NA} in all cases.
#'
#' @param data Patient-level data frame.
#' @param comorbidity_col Character. Default \code{"comorbidities"}.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with column \code{comorbidity_encoded} added.
#'   Attribute \code{"comorbidity_encoding"} records the strategy used:
#'   \code{"numeric"}, \code{"binary"}, \code{"ordinal"}, or \code{"absent"}.
#' @keywords internal
.encode_comorbidity_mortality <- function(
  data,
  comorbidity_col = "comorbidities",
  patient_id_col = "PatientInformation_id"
) {
  if (!comorbidity_col %in% names(data)) {
    message(sprintf(
      "[Comorbidity] Column '%s' not found -- covariate will be omitted.",
      comorbidity_col
    ))
    data$comorbidity_encoded <- NA_real_
    attr(data, "comorbidity_encoding") <- "absent"
    return(data)
  }

  raw <- data[[comorbidity_col]]

  # Strategy 1: already a numeric index (Charlson / Elixhauser)
  if (is.numeric(raw)) {
    data$comorbidity_encoded <- raw
    attr(data, "comorbidity_encoding") <- "numeric"
    message(sprintf(
      "[Comorbidity] Numeric index detected (range %.1f-%.1f). Used as-is.",
      min(raw, na.rm = TRUE), max(raw, na.rm = TRUE)
    ))
    return(data)
  }

  # Normalise to uppercase for pattern matching
  raw_up <- stringr::str_to_upper(stringr::str_trim(as.character(raw)))
  # Sentinel missing values
  na_vals <- c(
    "NA", "NULL", "", "UNKNOWN", "NOT KNOWN", "NOT_KNOWN",
    "MISSING", "N/A", "NONE RECORDED"
  )
  raw_up[raw_up %in% na_vals | is.na(raw)] <- NA_character_

  # Strategy 2: binary present / none
  binary_present <- c(
    "PRESENT", "YES", "1", "TRUE", "POSITIVE",
    "COMORBID", "COMORBIDITY PRESENT"
  )
  binary_none <- c(
    "NONE", "NO", "0", "FALSE", "ABSENT",
    "NO COMORBIDITY", "COMORBIDITY ABSENT"
  )
  non_na_vals <- raw_up[!is.na(raw_up)]

  if (length(non_na_vals) > 0L &&
    all(non_na_vals %in% c(binary_present, binary_none))) {
    data$comorbidity_encoded <- dplyr::case_when(
      raw_up %in% binary_present ~ 1L,
      raw_up %in% binary_none ~ 0L,
      TRUE ~ NA_integer_
    )
    n_present <- sum(data$comorbidity_encoded == 1L, na.rm = TRUE)
    n_none <- sum(data$comorbidity_encoded == 0L, na.rm = TRUE)
    n_miss <- sum(is.na(data$comorbidity_encoded))
    attr(data, "comorbidity_encoding") <- "binary"
    message(sprintf(
      "[Comorbidity] Binary encoding: present = %d | none = %d | NA = %d.",
      n_present, n_none, n_miss
    ))
    return(data)
  }

  # Strategy 3: ordinal none / mild / moderate / severe
  ord_none <- c("NONE", "NO COMORBIDITY", "ABSENT", "0")
  ord_mild <- c("MILD", "LOW", "MINOR", "MINIMAL", "1")
  ord_moderate <- c("MODERATE", "MEDIUM", "MODERATE COMORBIDITY", "2")
  ord_severe <- c(
    "SEVERE", "HIGH", "MAJOR", "HEAVY",
    "SIGNIFICANT", "MULTIPLE", "3"
  )

  data$comorbidity_encoded <- dplyr::case_when(
    raw_up %in% ord_none ~ "none",
    raw_up %in% ord_mild ~ "mild",
    raw_up %in% ord_moderate ~ "moderate",
    raw_up %in% ord_severe ~ "severe",
    TRUE ~ NA_character_
  )

  n_unrecog <- sum(
    !is.na(raw_up) &
      !(raw_up %in% c(ord_none, ord_mild, ord_moderate, ord_severe))
  )
  if (n_unrecog > 0L) {
    message(sprintf(
      "[Comorbidity] %d unrecognised value(s) set to NA. Examples: %s",
      n_unrecog,
      paste(head(unique(
        raw_up[!is.na(raw_up) &
          !(raw_up %in% c(
            ord_none, ord_mild,
            ord_moderate, ord_severe
          ))]
      ), 5L), collapse = ", ")
    ))
  }

  data$comorbidity_encoded <- factor(
    data$comorbidity_encoded,
    levels  = c("none", "mild", "moderate", "severe"),
    ordered = TRUE
  )
  attr(data, "comorbidity_encoding") <- "ordinal"
  message(sprintf(
    "[Comorbidity] Ordinal encoding: %s.",
    paste(
      paste0(
        levels(data$comorbidity_encoded), " = ",
        as.integer(table(data$comorbidity_encoded))
      ),
      collapse = " | "
    )
  ))
  return(data)
}


# Step M4
#' Check Collinearity Between HAI and ICU Covariates
#'
#' Computes the phi (Pearson) correlation coefficient for the 2x2 contingency
#' table of HAI x ICU. When ICU is hospital-acquired, the two covariates can
#' be highly correlated, making regression coefficients unstable.
#'
#' Perfect separation (any 2x2 marginal = 0) is detected and reported
#' separately, as it guarantees model non-convergence.
#'
#' @param df Patient-level data frame with integer 0/1 columns.
#' @param hai_col Character. Default \code{"HAI"}.
#' @param icu_col Character. Default \code{"ICU"}.
#' @param phi_threshold Numeric. Absolute phi above which a warning is issued.
#'   Default \code{0.7}.
#'
#' @return Named list: \code{phi}, \code{tbl} (2x2 table),
#'   \code{warning_issued} (logical).
#' @keywords internal
.check_hai_icu_collinearity <- function(df,
                                        hai_col = "HAI",
                                        icu_col = "ICU",
                                        phi_threshold = 0.7) {
  if (!hai_col %in% names(df) || !icu_col %in% names(df)) {
    return(list(phi = NA_real_, tbl = NULL, warning_issued = FALSE))
  }

  h <- as.integer(df[[hai_col]])
  ic <- as.integer(df[[icu_col]])
  ok <- !is.na(h) & !is.na(ic)
  h <- h[ok]
  ic <- ic[ok]

  if (length(h) < 10L) {
    message("[Collinearity] Too few complete observations to assess HAI/ICU correlation.")
    return(list(phi = NA_real_, tbl = NULL, warning_issued = FALSE))
  }

  tbl <- table(HAI = h, ICU = ic)
  message("[Collinearity] HAI \u00d7 ICU cross-table:")
  .msg_print(tbl)

  # Perfect separation: any marginal is zero
  if (any(rowSums(tbl) == 0L) || any(colSums(tbl) == 0L)) {
    warning(
      "[Collinearity] Perfect separation detected between HAI and ICU: ",
      "at least one cell combination has zero patients. ",
      "Model will likely fail to converge. ",
      "Consider dropping one of these two covariates."
    )
    return(list(phi = NA_real_, tbl = tbl, warning_issued = TRUE))
  }

  phi <- tryCatch(
    stats::cor(h, ic, method = "pearson"),
    error = function(e) NA_real_
  )

  warn_issued <- FALSE
  if (!is.na(phi) && abs(phi) >= phi_threshold) {
    warning(sprintf(
      paste0(
        "[Collinearity] HAI and ICU are highly correlated ",
        "(phi = %.3f \u2265 threshold %.2f). ",
        "Regression coefficients may be unstable. ",
        "Consider removing one covariate or checking your data."
      ),
      phi, phi_threshold
    ))
    warn_issued <- TRUE
  } else if (!is.na(phi)) {
    message(sprintf(
      "[Collinearity] HAI/ICU phi = %.3f (below threshold %.2f \u2014 OK).",
      phi, phi_threshold
    ))
  }

  list(phi = phi, tbl = tbl, warning_issued = warn_issued)
}


#' Fit mortality model and derive adjusted relative risk of death
#'
#' Fits a pathogen-class specific mixed-effects logistic regression model:
#'
#'   Death_i ~ Bernoulli(pi_i)
#'   logit(pi_i) =
#'     beta0 +
#'     beta_kd * Resistant_id +
#'     beta_age * Age_i +
#'     beta_sex * Sex_i +
#'     beta_hai * HAI_i +
#'     beta_icu * ICU_i +
#'     beta_comorb * Comorbidity_i +
#'     beta_syn * Syndrome_i +
#'     u_facility
#'
#' where u_facility is a random intercept when more than one facility is present.
#'
#' The function returns both:
#'   - OR_death = exp(beta_kd)
#'   - RR_death = mean(p_resistant_cf) / mean(p_susceptible_cf)
#'
#' RR_death is derived from model-based predicted probabilities under
#' resistant and susceptible counterfactual scenarios, as required by
#' the burden methodology.
#'
#' @param data Data frame.
#' @param patient_id_col Character. Patient identifier column.
#' @param facility_col Character or NULL. Facility / centre column. When
#'   \code{NULL}, no hospital random effect is used and a single pooled RR is
#'   returned using standard logistic regression (hospital column will be
#'   \code{NA}). When provided, a separate logistic model is fit per hospital
#'   and the result contains one row per pathogen-class-hospital combination
#'   (hospital-specific RR).
#' @param organism_col Character. Pathogen column.
#' @param syndrome_col Character. Syndrome column.
#' @param infection_type_col Character. Raw infection-type column.
#' @param antibiotic_class_col Character. Antibiotic class column.
#' @param antibiotic_name_col Character. Antibiotic name column.
#' @param antibiotic_value_col Character. Antibiotic susceptibility column.
#' @param unit_type_col Character or NULL. ICU / ward location column.
#' @param date_admission_col Character. Admission date column.
#' @param date_culture_col Character. First positive culture date column.
#' @param final_outcome_col Character. Final outcome column.
#' @param final_outcome_date_col Character. Final outcome date column.
#' @param age_col Character. Age column.
#' @param sex_col Character. Sex column.
#' @param comorbidity_col Character or NULL. Comorbidity column.
#' @param death_value Character. Value indicating death.
#' @param syndrome_name Character or NULL. Restrict to one syndrome.
#' @param organism_name Character vector or NULL. Restrict to these pathogens.
#' @param hai_threshold_hours Numeric. HAI derivation threshold.
#' @param icu_values Character vector. Values treated as ICU.
#' @param phi_threshold Numeric. HAI / ICU collinearity warning threshold.
#' @param min_n Integer. Minimum patients per fitted model.
#' @param min_deaths Integer. Minimum deaths per fitted model.
#'
#' @return Data frame with hospital (NA when \code{facility_col} is
#'   \code{NULL}), adjusted OR and adjusted RR of death.
#' @export
daly_fit_mortality_rr <- function(
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
  min_deaths = 10L
) {
  required_cols <- c(
    patient_id_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_culture_col,
    final_outcome_col, age_col, sex_col
  )

  if (!is.null(facility_col)) {
    required_cols <- c(required_cols, facility_col)
  }
  if (!is.null(syndrome_name)) {
    required_cols <- c(required_cols, syndrome_col)
  }

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Missing required column(s): %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # When facility_col is NULL, no hospital random intercept is used and a
  # pooled (standard logistic) RR is returned. A dummy column is added
  # internally so helper functions remain compatible.
  use_facility <- !is.null(facility_col)
  if (!use_facility) {
    facility_col <- ".facility_pooled"
    data[[facility_col]] <- "pooled"
    message("No facility_col provided: hospital random effect will not be included. Returning pooled RR.")
  }

  # Step 1: derive HAI / CAI for mortality model
  data <- daly_derive_hai_cai_for_mortality(
    data = data,
    infection_type_col = infection_type_col,
    date_admission_col = date_admission_col,
    date_culture_col = date_culture_col,
    final_outcome_col = final_outcome_col,
    final_outcome_date_col = final_outcome_date_col,
    death_value = death_value,
    hai_threshold_hours = hai_threshold_hours,
    patient_id_col = patient_id_col
  )

  # Step 2: optional syndrome restriction
  df <- data
  if (!is.null(syndrome_name)) {
    df <- df[.name_key(df[[syndrome_col]]) == .name_key(syndrome_name), , drop = FALSE]
  }

  if (nrow(df) == 0L) {
    warning("No rows remain after syndrome filtering.")
    return(data.frame())
  }

  # Step 3: ICU indicator
  icu_tbl <- .derive_icu_binary(
    data = df,
    patient_id_col = patient_id_col,
    unit_type_col = unit_type_col,
    icu_values = icu_values
  )

  # Step 4: patient-level covariates
  # one row per patient
  patient_covars <- df |>
    dplyr::group_by(.data[[patient_id_col]]) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      !!patient_id_col := .data[[patient_id_col]],
      !!facility_col := .data[[facility_col]],
      death = as.integer(
        stringr::str_to_upper(stringr::str_trim(as.character(.data[[final_outcome_col]]))) %in%
          stringr::str_to_upper(death_value)
      ),
      HAI = as.integer(.data[["infection_type_derived"]] == "HAI"),
      Age_model = suppressWarnings(as.numeric(.data[[age_col]])),
      Sex_model = as.factor(.data[[sex_col]])
    )

  if (syndrome_col %in% names(df)) {
    syndrome_tbl <- df |>
      dplyr::group_by(.data[[patient_id_col]]) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        !!patient_id_col := .data[[patient_id_col]],
        Syndrome_model = as.factor(.data[[syndrome_col]])
      )
    patient_covars <- dplyr::left_join(patient_covars, syndrome_tbl, by = patient_id_col)
  } else {
    patient_covars$Syndrome_model <- factor(NA_character_)
  }

  patient_covars <- dplyr::left_join(patient_covars, icu_tbl, by = patient_id_col)

  # Step 5: comorbidity encoding
  if (!is.null(comorbidity_col) && comorbidity_col %in% names(df)) {
    comorb_df <- df |>
      dplyr::group_by(.data[[patient_id_col]]) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(c(patient_id_col, comorbidity_col)))

    comorb_df <- .encode_comorbidity_mortality(
      data = comorb_df,
      comorbidity_col = comorbidity_col,
      patient_id_col = patient_id_col
    )

    patient_covars <- dplyr::left_join(
      patient_covars,
      comorb_df |>
        dplyr::select(.data[[patient_id_col]], comorbidity_encoded),
      by = patient_id_col
    )
  } else {
    patient_covars$comorbidity_encoded <- NA_real_
  }

  # Step 6: complete patient-level covariates. Class resistance is built per
  # pathogen inside the loop below from that pathogen's own results, so a
  # patient with several pathogens enters every pathogen's models and one
  # pathogen's results never set another pathogen's resistance.
  patient_covars <- patient_covars |>
    dplyr::filter(!is.na(death), !is.na(HAI), !is.na(Age_model), !is.na(Sex_model))

  if (nrow(patient_covars) == 0L) {
    warning("No complete patient-level rows available for mortality modelling.")
    return(data.frame())
  }

  org_present <- unique(as.character(df[[organism_col]]))
  org_present <- org_present[!is.na(org_present) & trimws(org_present) != ""]
  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(org_present)
  }

  out <- list()

  # helper: derive standardized RR from counterfactual predictions
  .compute_counterfactual_rr <- function(model, data_cf, resistance_var) {
    d1 <- data_cf
    d0 <- data_cf

    d1[[resistance_var]] <- 1L
    d0[[resistance_var]] <- 0L

    p1 <- stats::predict(model, newdata = d1, type = "response", allow.new.levels = TRUE)
    p0 <- stats::predict(model, newdata = d0, type = "response", allow.new.levels = TRUE)

    mean_p1 <- mean(p1, na.rm = TRUE)
    mean_p0 <- mean(p0, na.rm = TRUE)

    rr <- if (is.finite(mean_p0) && mean_p0 > 0) mean_p1 / mean_p0 else NA_real_

    list(
      rr = rr,
      mean_p_resistant = mean_p1,
      mean_p_susceptible = mean_p0
    )
  }

  for (path in pathogens) {
    # Step 7: class-level resistance from this pathogen's results only
    df_path <- df[.name_key(df[[organism_col]]) == .name_key(path), , drop = FALSE]
    if (nrow(df_path) == 0L) {
      next
    }
    resistance_wide <- .build_class_resistance_wide(
      data = df_path,
      patient_id_col = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col,
      untested_fill = NA_integer_
    )
    class_name_map <- attr(resistance_wide, "class_name_map")
    class_safe_cols <- setdiff(names(resistance_wide), patient_id_col)

    path_df <- patient_covars |>
      dplyr::inner_join(resistance_wide, by = patient_id_col)

    if (nrow(path_df) == 0L) {
      next
    }

    # collinearity check on pathogen-specific cohort
    if ("ICU" %in% names(path_df)) {
      invisible(.check_hai_icu_collinearity(
        df = path_df,
        hai_col = "HAI",
        icu_col = "ICU",
        phi_threshold = phi_threshold
      ))
    }

    for (cls in class_safe_cols) {
      orig_class <- class_name_map[[cls]]

      sub <- path_df |>
        dplyr::filter(!is.na(.data[[cls]]))

      if (nrow(sub) == 0L) {
        next
      }

      # When facility_col is provided, fit a hospital-specific logistic model
      # per hospital (no random intercept needed within a single hospital).
      # When facility_col is NULL (pooled), fit one standard logistic model.
      hospitals <- if (use_facility) {
        sort(unique(sub[[facility_col]]))
      } else {
        NA_character_
      }

      for (hosp in hospitals) {
        hosp_sub <- if (use_facility) {
          sub[sub[[facility_col]] == hosp, ]
        } else {
          sub
        }

        if (nrow(hosp_sub) < min_n) {
          next
        }

        n_deaths <- sum(hosp_sub$death == 1L, na.rm = TRUE)
        n_survivors <- sum(hosp_sub$death == 0L, na.rm = TRUE)
        n_resistant <- sum(hosp_sub[[cls]] == 1L, na.rm = TRUE)
        n_susceptible <- sum(hosp_sub[[cls]] == 0L, na.rm = TRUE)

        if (n_deaths < min_deaths || n_survivors < min_deaths) {
          next
        }

        if (n_resistant == 0L || n_susceptible == 0L) {
          next
        }

        fixed_terms <- c(sprintf("`%s`", cls), "Age_model", "Sex_model", "HAI")

        if ("ICU" %in% names(hosp_sub) &&
          length(unique(stats::na.omit(as.character(hosp_sub$ICU)))) > 1L) {
          fixed_terms <- c(fixed_terms, "ICU")
        }

        if ("comorbidity_encoded" %in% names(hosp_sub) &&
          length(unique(stats::na.omit(as.character(hosp_sub$comorbidity_encoded)))) > 1L) {
          fixed_terms <- c(fixed_terms, "comorbidity_encoded")
        }

        if ("Syndrome_model" %in% names(hosp_sub) &&
          length(unique(stats::na.omit(as.character(hosp_sub$Syndrome_model)))) > 1L &&
          is.null(syndrome_name)) {
          fixed_terms <- c(fixed_terms, "Syndrome_model")
        }

        # No random intercept: when use_facility is TRUE the model is already
        # restricted to one hospital; when use_facility is FALSE there is no
        # facility information to include as a random effect.
        fmla <- stats::as.formula(
          paste("death ~", paste(fixed_terms, collapse = " + "))
        )

        fit <- tryCatch(
          suppressWarnings(
            stats::glm(
              formula = fmla,
              data = hosp_sub,
              family = stats::binomial(link = "logit")
            )
          ),
          error = function(e) NULL
        )

        if (is.null(fit)) {
          next
        }

        coefs <- stats::coef(fit)
        vcov_mat <- as.matrix(stats::vcov(fit))

        cls_pattern <- paste0("^`?", gsub(".", "\\\\.", cls, fixed = TRUE), "`?$")
        coef_name <- grep(cls_pattern, names(coefs), value = TRUE)

        if (length(coef_name) == 0L) {
          next
        }

        beta <- coefs[[coef_name[1L]]]
        se <- sqrt(vcov_mat[coef_name[1L], coef_name[1L]])

        or_death <- exp(beta)
        or_lower <- exp(beta - 1.96 * se)
        or_upper <- exp(beta + 1.96 * se)

        rr_obj <- .compute_counterfactual_rr(
          model = fit,
          data_cf = hosp_sub,
          resistance_var = cls
        )

        out[[length(out) + 1L]] <- data.frame(
          pathogen = path,
          antibiotic_class = orig_class,
          hospital = if (use_facility) hosp else NA_character_,
          beta_logit = unname(beta),
          OR_death = unname(or_death),
          OR_CI_lower = unname(or_lower),
          OR_CI_upper = unname(or_upper),
          RR_death = unname(rr_obj$rr),
          mean_p_resistant_cf = unname(rr_obj$mean_p_resistant),
          mean_p_susceptible_cf = unname(rr_obj$mean_p_susceptible),
          n_patients = nrow(hosp_sub),
          n_deaths = n_deaths,
          n_resistant = n_resistant,
          n_susceptible = n_susceptible,
          model_type = "logistic",
          syndrome_scope = if (!is.null(syndrome_name)) syndrome_name else "all",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(out) == 0L) {
    warning("No mortality models could be fitted.")
    return(data.frame())
  }

  result <- dplyr::bind_rows(out)

  message(sprintf(
    "Mortality models fitted: %d pathogen-class estimates returned.",
    nrow(result)
  ))

  result
}

# Length-of-stay distribution fitting and comparison

#' Safely Fit a Distribution
#'
#' Wrapper around \code{fitdistrplus::fitdist} that returns \code{NULL}
#' instead of throwing an error when fitting fails.
#'
#' @param x Numeric vector of positive values (e.g., LOS in days).
#' @param dist Character. Distribution name: \code{"weibull"}, \code{"lnorm"},
#'   or \code{"gamma"}.
#'
#' @return A \code{fitdist} object, or \code{NULL} on failure.
#' @export
#'
#' @examples
#' safe_fit(rlnorm(100), "lnorm")
safe_fit <- function(x, dist) {
  x <- x[is.finite(x)]
  if (length(x) < 5 || length(unique(x)) < 2) {
    return(NULL)
  }
  tryCatch(fitdistrplus::fitdist(x, dist), error = function(e) NULL)
}


#' Fit Multiple Distributions
#'
#' Fits Weibull, Lognormal, and Gamma distributions to a numeric vector,
#' returning all fit objects. Failed fits are \code{NULL}.
#'
#' @param x Numeric vector of positive values (e.g., LOS in days).
#'
#' @return Named list with elements \code{weibull}, \code{lnorm}, \code{gamma},
#'   each a \code{fitdist} object or \code{NULL}.
#' @export
#'
#' @examples
#' fits <- fit_distributions(rlnorm(200, 2, 0.5))
#' fits$lnorm$estimate
fit_distributions <- function(x) {
  list(
    weibull = safe_fit(x, "weibull"),
    lnorm   = safe_fit(x, "lnorm"),
    gamma   = safe_fit(x, "gamma")
  )
}


#' Compare Distribution Fits by AIC
#'
#' Fits Weibull, Lognormal, and Gamma distributions to a numeric vector and
#' returns their AIC values for comparison.
#'
#' @param x Numeric vector of positive values (e.g., LOS in days).
#' @param fits Optional. Pre-computed fits from \code{fit_distributions()}.
#'
#' @return Data frame with columns \code{Weibull_AIC}, \code{Lognormal_AIC},
#'   \code{Gamma_AIC}.
#' @export
#'
#' @examples
#' compare_distribution_aic(rlnorm(200, meanlog = 2, sdlog = 0.5))
compare_distribution_aic <- function(x, fits = NULL) {
  if (is.null(fits)) fits <- fit_distributions(x)

  data.frame(
    Weibull_AIC   = if (is.null(fits$weibull)) NA_real_ else fits$weibull$aic,
    Lognormal_AIC = if (is.null(fits$lnorm)) NA_real_ else fits$lnorm$aic,
    Gamma_AIC     = if (is.null(fits$gamma)) NA_real_ else fits$gamma$aic
  )
}


#' Summarise a Fitted Distribution
#'
#' Extracts mean, median, SD, and parameter values from a \code{fitdist}
#' object for Weibull, Lognormal, or Gamma distributions.
#'
#' @param fit A \code{fitdist} object from \code{fitdistrplus::fitdist()}.
#' @param dist Character. One of \code{"weibull"}, \code{"lnorm"}, or
#'   \code{"gamma"}.
#'
#' @return Data frame with columns \code{Mean_LOS}, \code{Median_LOS},
#'   \code{SD_LOS}, \code{Parameters}.
#' @export
#'
#' @examples
#' fit <- fitdistrplus::fitdist(rlnorm(200), "lnorm")
#' summarise_distribution(fit, "lnorm")
summarise_distribution <- function(fit, dist) {
  if (dist == "weibull") {
    k <- fit$estimate["shape"]
    l <- fit$estimate["scale"]
    mean_val <- l * gamma(1 + 1 / k)
    median_val <- l * (log(2))^(1 / k)
    sd_val <- sqrt(l^2 * (gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2))
    params <- paste0("shape=", round(k, 3), ", scale=", round(l, 3))
  } else if (dist == "lnorm") {
    mu <- fit$estimate["meanlog"]
    sg <- fit$estimate["sdlog"]
    mean_val <- exp(mu + sg^2 / 2)
    median_val <- exp(mu)
    sd_val <- sqrt((exp(sg^2) - 1) * exp(2 * mu + sg^2))
    params <- paste0("meanlog=", round(mu, 3), ", sdlog=", round(sg, 3))
  } else if (dist == "gamma") {
    a <- fit$estimate["shape"]
    r <- fit$estimate["rate"]
    mean_val <- a / r
    median_val <- stats::qgamma(0.5, shape = a, rate = r)
    sd_val <- sqrt(a) / r
    params <- paste0("shape=", round(a, 3), ", rate=", round(r, 3))
  } else {
    stop("dist must be one of 'weibull', 'lnorm', 'gamma'")
  }

  data.frame(
    Mean_LOS   = round(mean_val, 2),
    Median_LOS = round(median_val, 2),
    SD_LOS     = round(sd_val, 2),
    Parameters = params
  )
}


#' Plot LOS Distribution with Fitted Overlays
#'
#' Creates a histogram of LOS values with Weibull, Lognormal, and Gamma
#' density curves overlaid.
#'
#' @param los_vec Numeric vector of LOS values (days).
#' @param title Character. Plot title.
#' @param bins Integer. Number of histogram bins. Default 35.
#' @param fits Optional. Pre-computed fits from \code{fit_distributions()}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' plot_los_distributions(rlnorm(200, 2, 0.5), "Example LOS Distribution")
plot_los_distributions <- function(los_vec, title, bins = 35, fits = NULL) {
  if (is.null(fits)) fits <- fit_distributions(los_vec)

  overlay_cfg <- list(
    weibull = list(
      dfun = stats::dweibull, color = "red",
      linetype = "solid", label = "Weibull"
    ),
    lnorm = list(
      dfun = stats::dlnorm, color = "blue",
      linetype = "dashed", label = "Lognormal"
    ),
    gamma = list(
      dfun = stats::dgamma, color = "darkgreen",
      linetype = "dotdash", label = "Gamma"
    )
  )

  p <- ggplot2::ggplot(data.frame(los = los_vec), ggplot2::aes(x = .data$los)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = "#d9d9d9", color = "black"
    )

  fitted_labels <- character(0)
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (is.null(fit)) next
    cfg <- overlay_cfg[[nm]]
    p <- p + ggplot2::stat_function(
      fun = cfg$dfun,
      args = as.list(fit$estimate),
      color = cfg$color, linetype = cfg$linetype, linewidth = 1.2
    )
    fitted_labels <- c(
      fitted_labels,
      paste0(cfg$label, " (", cfg$color, ")")
    )
  }

  subtitle <- if (length(fitted_labels) > 0) {
    paste(fitted_labels, collapse = " | ")
  } else {
    "No distributions could be fitted"
  }

  p + ggplot2::labs(
    title = title, subtitle = subtitle,
    x = "LOS (days)", y = "Density"
  ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 18),
      plot.subtitle = ggplot2::element_text(size = 14, color = "grey40")
    )
}


#' Extract LOS Vectors by Resistance Status
#'
#' For a given organism (and optionally center), computes per-patient resistance
#' status and returns separate LOS vectors for resistant and susceptible
#' exposures.
#'
#' @param abx_data Data frame of antibiotic susceptibility results with columns
#'   for patient ID, organism, antibiotic name, and susceptibility value
#'   (already cleaned to "R"/"S").
#' @param los_data Data frame from \code{\link{daly_compute_patient_los}()}
#'   (columns \code{patient_id_col}, \code{LOS_days}, and \code{center_col}
#'   when \code{center} is supplied).
#' @param organism Character. Organism name (lowercase, trimmed) to filter on.
#' @param center Optional character. Center name to filter on. Default
#'   \code{NULL} (all centers).
#' @param patient_id_col Character. Patient ID column. Default
#'   \code{"PatientInformation_id"}.
#' @param organism_col Character. Organism column. Default
#'   \code{"organism_clean"}.
#' @param center_col Character. Center column. Default \code{"center_name"}.
#' @param antibiotic_col Character. Antibiotic name column. Default
#'   \code{"antibiotic_name"}.
#' @param sus_col Character. Susceptibility column (values "R"/"S"). Default
#'   \code{"sus"}.
#'
#' @return List with elements \code{R} (LOS vector for resistant exposures) and
#'   \code{S} (LOS vector for susceptible exposures).
#' @export
get_los_by_resistance <- function(abx_data,
                                  los_data,
                                  organism,
                                  center = NULL,
                                  patient_id_col = "PatientInformation_id",
                                  organism_col = "organism_clean",
                                  center_col = "center_name",
                                  antibiotic_col = "antibiotic_name",
                                  sus_col = "sus") {
  org_abx <- abx_data %>%
    dplyr::filter(.name_key(.data[[organism_col]]) == .name_key(organism))

  if (!is.null(center)) {
    org_abx <- dplyr::filter(org_abx, .data[[center_col]] == center)
    los_data <- dplyr::filter(los_data, .data[[center_col]] == center)
  }

  patient_resistance <- org_abx %>%
    dplyr::group_by(.data[[patient_id_col]], .data[[antibiotic_col]]) %>%
    dplyr::summarise(
      final_sus = ifelse(any(.data[[sus_col]] == "R"), "R", "S"),
      .groups = "drop"
    ) %>%
    dplyr::count(.data[[patient_id_col]], .data$final_sus) %>%
    tidyr::pivot_wider(
      names_from = "final_sus", values_from = "n",
      values_fill = 0
    )

  # Handle case where R or S column may not exist
  if (!"R" %in% names(patient_resistance)) {
    patient_resistance$R <- 0L
  }
  if (!"S" %in% names(patient_resistance)) {
    patient_resistance$S <- 0L
  }

  patient_resistance <- dplyr::rename(patient_resistance,
    R_count = "R", S_count = "S"
  )

  org_los <- dplyr::inner_join(los_data, patient_resistance,
    by = patient_id_col
  )

  # rep() instead of tidyr::uncount(): same expanded LOS vectors, no tibble round-trip.
  los_R <- rep(org_los$LOS_days[org_los$R_count > 0], org_los$R_count[org_los$R_count > 0])
  los_S <- rep(org_los$LOS_days[org_los$S_count > 0], org_los$S_count[org_los$S_count > 0])

  list(R = los_R, S = los_S)
}


# Relative risk reference data, pathogen-drug mappings, and profile-level RR
# assignment. Profile RR uses the GBD max rule: the highest class-level RR
# among a profile's resistant classes, with profiles re-normalised after
# dropping any whose resistant classes carry no RR estimate.

#' Load RR (Relative Risk) Reference Data
#'
#' Loads the bundled rr_list_gbd.csv file containing pathogen-drug combinations
#' with relative risk values for burden estimation.
#'
#' @return Data frame with Pathogen, Drug, and RR columns
#' @export
#'
#' @examples
#' rr_data <- daly_load_rr_reference()
#' head(rr_data)
daly_load_rr_reference <- function() {
  csv_path <- find_extdata_file("rr_list_gbd.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    stop("rr_list_gbd.csv not found in inst/extdata/")
  }

  rr_data <- utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Clean column names
  names(rr_data) <- c(
    "pathogen",
    "drug",
    "sample_size",
    "mean_rr",
    "lower_bound",
    "upper_bound"
  )

  # Remove rows with NA pathogen (like footer notes)
  rr_data <- rr_data %>%
    dplyr::filter(!is.na(pathogen), !is.na(drug))
  rr_data$pathogen <- .name_key(rr_data$pathogen) # matches rr_pathogen keys

  message(sprintf(
    "[v] Loaded %d pathogen-drug RR combinations",
    nrow(rr_data)
  ))

  return(rr_data)
}


#' Add RR Pathogen and Drug Class Mappings
#'
#' Maps organism names to GBD RR pathogen categories and antibiotic class names
#' to RR drug categories in a single pass. Either mapping is skipped silently
#' when the corresponding column is absent from \code{data}.
#'
#' Pathogen mapping uses the organism taxonomy (\code{get_organism_taxonomy()});
#' drug class mapping uses \code{inst/extdata/WHO_aware_class.csv}.
#'
#' @param data Data frame with organism and/or antibiotic class columns.
#' @param organism_col Character. Column with organism names to map to
#'   \code{rr_pathogen}. Set to \code{NULL} to skip. Default
#'   \code{"organism_normalized"}.
#' @param class_col Character. Column with antibiotic class names to map to
#'   \code{rr_drug}. Set to \code{NULL} to skip. Default
#'   \code{"antibiotic_class"}.
#'
#' @return Data frame with \code{rr_pathogen} and/or \code{rr_drug} columns
#'   appended.
#' @export
daly_add_rr_mappings <- function(data,
                                 organism_col = "organism_normalized",
                                 class_col = "antibiotic_class") {
  # Pathogen mapping
  if (!is.null(organism_col) && organism_col %in% names(data)) {
    taxonomy <- get_organism_taxonomy()
    if (nrow(taxonomy) > 0) {
      key <- .name_key(data[[organism_col]])
      data$rr_pathogen <- ifelse(key %in% taxonomy$organism_name, key, NA_character_)
      n_mapped <- sum(!is.na(data$rr_pathogen))
      message(sprintf(
        "Mapped to RR pathogen: %d/%d (%.1f%%)",
        n_mapped, nrow(data), 100 * n_mapped / nrow(data)
      ))
      unmapped <- unique(data[[organism_col]][is.na(data$rr_pathogen)])
      if (length(unmapped) > 0) {
        message("Unmapped organisms: ", paste(head(unmapped, 10), collapse = ", "))
      }
    } else {
      warning("Organism taxonomy is empty. Skipping pathogen mapping.")
    }
  }

  # Drug class mapping
  if (!is.null(class_col) && class_col %in% names(data)) {
    csv_path <- find_extdata_file("WHO_aware_class.csv")
    if (csv_path == "") {
      warning("WHO_aware_class.csv not found. Skipping drug class mapping.")
    } else {
      who <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
      if (!"rr_drug" %in% names(who)) who$rr_drug <- who$Class
      class_map <- unique(who[, c("Class", "rr_drug")])
      data <- dplyr::left_join(
        data, class_map,
        by = stats::setNames("Class", class_col)
      )
      n_mapped <- sum(!is.na(data$rr_drug))
      message(sprintf(
        "Mapped to RR drug: %d/%d (%.1f%%)",
        n_mapped, nrow(data), 100 * n_mapped / nrow(data)
      ))
    }
  }

  return(data)
}


#' Assign Per-Class LOS RR to Resistance Profiles (Max Rule)
#'
#' For each resistance profile delta, determines the profile-level relative
#' LOS using the GBD max rule:
#'   RR_kd_LOS = max over c in C_R(d) of RR_kc_LOS   [if C_R(d) non-empty]
#'             = 1                                      [if d = all-susceptible]
#' where C_R(d) is the set of classes the profile is resistant to. The class
#' attaining the maximum is the profile's dominant class d*. The CI reported
#' for each profile is that of its dominant class.
#'
#' \strong{Ties}: when several resistant classes share the maximum RR, the
#' first of them in the class order wins (\code{profile_class_order} for
#' string profiles, the \code{classes} element for legacy profiles). A
#' message reports how many profiles had a tie.
#'
#' If \code{rr_table} was fit with a \code{syndrome_name} filter, its RR
#' values are syndrome-specific but get applied here to all profiles of a
#' pathogen -- i.e. it assumes syndrome-invariant LOS prolongation.
#'
#' @section Profile input formats:
#' \describe{
#'   \item{Legacy}{The named list from \code{compute_resistance_profiles()}
#'     (\code{$profiles} with 0/1 class columns and \code{$classes}).}
#'   \item{Profile strings}{A data frame with one row per (hospital x)
#'     pathogen x profile, where the profile is a string of R/S letters such
#'     as \code{"RRSR"} -- e.g. the output of
#'     \code{profile_aggregate_for_daly()} (\code{group}, \code{pathogen},
#'     \code{profile_delta}, \code{R_mean}). \code{profile_class_order} gives
#'     the class of each letter.}
#' }
#'
#' @section LOS for each profile:
#' When \code{rr_table} is the list returned by \code{daly_fit_los_rr()},
#' the profile takes the adjusted LOS of its dominant class: from
#' \code{rr_table$hospital} for the profile's own hospital when the profiles
#' carry a hospital column (\code{los_level = "hospital"}), otherwise the
#' all-hospital values in \code{rr_table$rr} (\code{los_level = "network"}).
#' Days are converted to years. With a data-frame \code{rr_table}, absolute
#' LOS is taken from \code{los_r_col} / \code{los_s_col} as before.
#' \code{LOSR_years} / \code{LOSS_years} are \code{NA} for the
#' all-susceptible profile, which has no dominant class.
#'
#' @section Two options for profile-level LOS:
#' This function is \strong{option 2}: per-class LOS estimates mapped to
#' profiles by the max rule. \strong{Option 1} is
#' \code{\link{daly_fit_los_distn_by_profile}()} (Section 5.3), which estimates LOS
#' directly for each patient-level resistance profile against the fully
#' susceptible profile.
#'
#' @param profiles_output Legacy list from \code{compute_resistance_profiles()}
#'   or a data frame of profile strings (see "Profile input formats").
#' @param rr_table The list returned by \code{daly_fit_los_rr()}, or a data
#'   frame with columns \code{pathogen_col}, \code{class_col}, \code{rr_col}
#'   and optionally \code{CI_lower} / \code{CI_upper} (e.g. from
#'   \code{daly_fit_los_rr_distribution()}).
#' @param pathogen_col Character. Pathogen column in \code{rr_table} (and in
#'   string profiles). Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param rr_col Character. Default \code{"RR_LOS"}.
#' @param fallback_rr Numeric. RR for resistant classes with no match.
#'   Default \code{1} (no attributable effect).
#' @param los_r_col,los_s_col Character or \code{NULL}. Columns in a
#'   data-frame \code{rr_table} with absolute LOS for resistant / susceptible
#'   patients (e.g. \code{"median_LOS_R"} / \code{"median_LOS_S"}). Ignored
#'   when \code{rr_table} is a \code{daly_fit_los_rr()} list.
#' @param los_unit Character. Unit of \code{los_r_col} / \code{los_s_col}:
#'   \code{"years"} or \code{"days"}. Default \code{"years"}.
#' @param profile_col Character. Profile-string column. Default
#'   \code{"profile_delta"}.
#' @param profile_probability_col Character. Profile probability column in
#'   string profiles. Default \code{"R_mean"}. Copied to \code{probability}.
#' @param profile_hospital_col Character or \code{NULL}. Hospital column in
#'   string profiles. Default \code{"group"}; \code{NULL} or absent means
#'   pooled profiles. Copied to \code{hospital}.
#' @param profile_class_order Character vector (the class of each letter, in
#'   order) or a named list of such vectors, one per pathogen. Required for
#'   string profiles.
#'
#' @return Named list (one entry per pathogen): the profiles data frame
#'   augmented with \code{RR_LOS_profile}, \code{dominant_class}, (if
#'   available) \code{CI_lower_profile} / \code{CI_upper_profile},
#'   \code{LOSR_years} / \code{LOSS_years}, and for string profiles
#'   \code{probability}, \code{hospital} and \code{los_level}.
#' @export
daly_assign_rr_to_profiles <- function(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1,
  los_r_col = NULL,
  los_s_col = NULL,
  los_unit = c("years", "days"),
  profile_col = "profile_delta",
  profile_probability_col = "R_mean",
  profile_hospital_col = "group",
  profile_class_order = NULL
) {
  los_unit <- match.arg(los_unit)

  # daly_fit_los_rr() list: RR from $rr, adjusted LOS from $hospital / $rr (days)
  los_hosp <- NULL
  fit_list <- is.list(rr_table) && !is.data.frame(rr_table) &&
    all(c("rr", "hospital") %in% names(rr_table))
  if (fit_list) {
    los_hosp <- as.data.frame(rr_table$hospital)
    rr_table <- as.data.frame(rr_table$rr)
    pathogen_col <- "pathogen"
    class_col <- "antibiotic_class"
    rr_col <- "RR_LOS"
    los_r_col <- "LOSR"
    los_s_col <- "LOSS"
    los_unit <- "days"
  }

  missing_rr <- setdiff(c(pathogen_col, class_col, rr_col), names(rr_table))
  if (length(missing_rr) > 0L) {
    stop(sprintf(
      "Column(s) not found in rr_table: %s",
      paste(missing_rr, collapse = ", ")
    ))
  }
  do_los <- !is.null(los_r_col) || !is.null(los_s_col)
  if (do_los) {
    if (is.null(los_r_col) || is.null(los_s_col)) {
      stop("los_r_col and los_s_col must both be supplied together, or both left NULL.")
    }
    missing_los <- setdiff(c(los_r_col, los_s_col), names(rr_table))
    if (length(missing_los) > 0L) {
      stop(sprintf(
        "Column(s) not found in rr_table: %s",
        paste(missing_los, collapse = ", ")
      ))
    }
  }
  to_years <- function(x) if (los_unit == "days") x / 365 else x
  has_ci <- all(c("CI_lower", "CI_upper") %in% names(rr_table))

  if (is.data.frame(profiles_output)) {
    return(.assign_rr_string_profiles(
      profiles_output, rr_table, los_hosp, pathogen_col, class_col, rr_col,
      fallback_rr, do_los, los_r_col, los_s_col, to_years, has_ci,
      profile_col, profile_probability_col, profile_hospital_col,
      profile_class_order
    ))
  }

  if (!is.list(profiles_output) ||
    !all(sapply(
      profiles_output,
      function(x) all(c("profiles", "classes") %in% names(x))
    ))) {
    stop(
      "profiles_output must be the list returned by compute_resistance_profiles() ",
      "or a data frame of profile strings."
    )
  }

  out <- list()
  n_ties <- 0L

  for (path in names(profiles_output)) {
    prof_list <- profiles_output[[path]]
    profiles <- prof_list$profiles
    classes <- prof_list$classes # alphabetical; matches profile cols

    rr_k <- rr_table[.name_key(rr_table[[pathogen_col]]) == .name_key(path), , drop = FALSE]
    rr_lookup <- setNames(rr_k[[rr_col]], rr_k[[class_col]])
    ci_lo_lookup <- if (has_ci) setNames(rr_k$CI_lower, rr_k[[class_col]]) else NULL
    ci_hi_lookup <- if (has_ci) setNames(rr_k$CI_upper, rr_k[[class_col]]) else NULL
    los_r_lookup <- if (do_los) setNames(rr_k[[los_r_col]], rr_k[[class_col]]) else NULL
    los_s_lookup <- if (do_los) setNames(rr_k[[los_s_col]], rr_k[[class_col]]) else NULL

    n_prof <- nrow(profiles)
    rr_profile <- numeric(n_prof)
    dom_class <- character(n_prof)
    ci_lo_prof <- if (has_ci) numeric(n_prof) else NULL
    ci_hi_prof <- if (has_ci) numeric(n_prof) else NULL
    losr_prof <- if (do_los) rep(NA_real_, n_prof) else NULL
    loss_prof <- if (do_los) rep(NA_real_, n_prof) else NULL

    for (i in seq_len(n_prof)) {
      resist_cls <- classes[as.integer(profiles[i, classes]) == 1L]

      if (length(resist_cls) == 0L) {
        rr_profile[i] <- 1.0
        dom_class[i] <- "all_susceptible"
        if (has_ci) {
          ci_lo_prof[i] <- 1.0
          ci_hi_prof[i] <- 1.0
        }
        # LOSR_years/LOSS_years left NA: no dominant class for this profile.
        next
      }

      rrs_c <- sapply(resist_cls, function(cc) {
        ifelse(cc %in% names(rr_lookup), rr_lookup[[cc]], fallback_rr)
      })
      max_idx <- which.max(rrs_c) # ties: first in class order
      if (sum(rrs_c == rrs_c[max_idx]) > 1L) n_ties <- n_ties + 1L
      rr_profile[i] <- rrs_c[max_idx]
      dom_class[i] <- resist_cls[max_idx]

      if (has_ci) {
        dc <- dom_class[i]
        ci_lo_prof[i] <- ifelse(dc %in% names(ci_lo_lookup),
          ci_lo_lookup[[dc]], fallback_rr
        )
        ci_hi_prof[i] <- ifelse(dc %in% names(ci_hi_lookup),
          ci_hi_lookup[[dc]], fallback_rr
        )
      }

      if (do_los) {
        dc <- dom_class[i]
        losr_prof[i] <- if (dc %in% names(los_r_lookup)) los_r_lookup[[dc]] else NA_real_
        loss_prof[i] <- if (dc %in% names(los_s_lookup)) los_s_lookup[[dc]] else NA_real_
      }
    }

    profiles$RR_LOS_profile <- round(rr_profile, 4L)
    profiles$dominant_class <- dom_class
    if (has_ci) {
      profiles$CI_lower_profile <- round(ci_lo_prof, 4L)
      profiles$CI_upper_profile <- round(ci_hi_prof, 4L)
    }
    if (do_los) {
      profiles$LOSR_years <- round(to_years(losr_prof), 6L)
      profiles$LOSS_years <- round(to_years(loss_prof), 6L)
    }

    out[[path]] <- profiles
    message(sprintf(
      "'%s': relative LOS assigned to %d profiles. Max relative LOS = %.4f (dominant class: %s).%s",
      path, n_prof, max(rr_profile), dom_class[which.max(rr_profile)],
      if (do_los) {
        sprintf(
          " Absolute LOS also assigned (%d/%d resistant profiles with a matched dominant-class LOS).",
          sum(!is.na(losr_prof) & dom_class != "all_susceptible"),
          sum(dom_class != "all_susceptible")
        )
      } else {
        ""
      }
    ))
  }
  if (n_ties > 0L) {
    message(sprintf("%d profile(s) had tied maximum RR; the first class in class order was used.", n_ties))
  }

  return(out)
}

# Max-rule assignment for profiles given as R/S strings (e.g. "RRSR"), one row
# per (hospital x) pathogen x profile. Returns a named list per pathogen.
.assign_rr_string_profiles <- function(prof, rr_table, los_hosp, pathogen_col,
                                       class_col, rr_col, fallback_rr, do_los,
                                       los_r_col, los_s_col, to_years, has_ci,
                                       profile_col, probability_col, hospital_col,
                                       class_order) {
  prof <- as.data.frame(prof, stringsAsFactors = FALSE)
  if (is.null(class_order)) {
    stop("profile_class_order is required for profile strings: the class of each letter, in order.")
  }
  need <- c(profile_col, probability_col, pathogen_col)
  miss <- setdiff(need, names(prof))
  if (length(miss)) {
    stop(sprintf("Column(s) not found in profiles_output: %s", paste(miss, collapse = ", ")))
  }
  has_hosp <- !is.null(hospital_col) && hospital_col %in% names(prof)
  prof$probability <- as.numeric(prof[[probability_col]])
  prof$hospital <- if (has_hosp) as.character(prof[[hospital_col]]) else NA_character_

  # RR is shared per pathogen x class, or per hospital for stratified fits
  rr_hosp <- if ("hospital" %in% names(rr_table)) as.character(rr_table$hospital) else rep(NA_character_, nrow(rr_table))
  rr_key <- function(p, c, h) paste(.name_key(p), .name_key(c), ifelse(is.na(h), "", h), sep = "\r")
  rr_k <- rr_key(rr_table[[pathogen_col]], rr_table[[class_col]], rr_hosp)
  lk <- function(col, p, c, h) {
    i <- match(rr_key(p, c, h), rr_k)
    if (is.na(i)) i <- match(rr_key(p, c, NA_character_), rr_k)
    if (is.na(i)) NA_real_ else rr_table[[col]][i]
  }
  hl_key <- if (!is.null(los_hosp)) {
    paste(.name_key(los_hosp$pathogen), .name_key(los_hosp$antibiotic_class), los_hosp$hospital, sep = "\r")
  } else {
    character()
  }

  out <- list()
  n_ties <- 0L; n_net <- 0L
  for (path in unique(as.character(prof[[pathogen_col]]))) {
    order_p <- class_order
    if (is.list(class_order)) {
      hit <- which(.name_key(names(class_order)) == .name_key(path))
      if (!length(hit)) stop(sprintf("profile_class_order has no entry for pathogen '%s'.", path))
      order_p <- class_order[[hit[1L]]]
    }
    df <- prof[as.character(prof[[pathogen_col]]) == path, , drop = FALSE]
    n <- nrow(df)
    rr_p <- rep(NA_real_, n); dom <- character(n)
    ci_lo <- ci_hi <- rep(NA_real_, n)
    losr <- loss <- rep(NA_real_, n); lvl <- rep(NA_character_, n)
    for (i in seq_len(n)) {
      letters_i <- strsplit(toupper(trimws(df[[profile_col]][i])), "")[[1]]
      if (length(letters_i) != length(order_p) || !all(letters_i %in% c("R", "S"))) {
        stop(sprintf(
          "Profile '%s' for '%s' does not match profile_class_order (%d classes).",
          df[[profile_col]][i], path, length(order_p)
        ))
      }
      res_cls <- order_p[letters_i == "R"]
      h <- df$hospital[i]
      if (!length(res_cls)) {
        rr_p[i] <- 1; dom[i] <- "all_susceptible"; ci_lo[i] <- 1; ci_hi[i] <- 1
        next
      }
      rrs <- vapply(res_cls, function(cc) {
        v <- lk(rr_col, path, cc, h)
        if (is.na(v)) fallback_rr else v
      }, numeric(1))
      j <- which.max(rrs) # ties: first in profile_class_order
      if (sum(rrs == rrs[j]) > 1L) n_ties <- n_ties + 1L
      rr_p[i] <- rrs[j]; dom[i] <- res_cls[j]
      if (has_ci) {
        ci_lo[i] <- lk("CI_lower", path, dom[i], h)
        ci_hi[i] <- lk("CI_upper", path, dom[i], h)
      }
      if (do_los) {
        hi <- if (!is.null(los_hosp) && !is.na(h)) match(paste(.name_key(path), .name_key(dom[i]), h, sep = "\r"), hl_key) else NA_integer_
        if (!is.na(hi)) {
          losr[i] <- los_hosp$LOSR[hi]; loss[i] <- los_hosp$LOSS[hi]; lvl[i] <- "hospital"
        } else {
          losr[i] <- lk(los_r_col, path, dom[i], h); loss[i] <- lk(los_s_col, path, dom[i], h)
          lvl[i] <- "network"
          if (!is.null(los_hosp) && !is.na(h)) n_net <- n_net + 1L
        }
      }
    }
    df$RR_LOS_profile <- round(rr_p, 4L)
    df$dominant_class <- dom
    if (has_ci) {
      df$CI_lower_profile <- round(ci_lo, 4L)
      df$CI_upper_profile <- round(ci_hi, 4L)
    }
    if (do_los) {
      df$LOSR_years <- round(to_years(losr), 6L)
      df$LOSS_years <- round(to_years(loss), 6L)
      df$los_level <- lvl
    }
    out[[path]] <- df
    message(sprintf(
      "'%s': relative LOS assigned to %d profile row(s) across %d hospital(s).",
      path, n, length(unique(df$hospital))
    ))
  }
  if (n_ties > 0L) {
    message(sprintf("%d profile row(s) had tied maximum RR; the first class in profile_class_order was used.", n_ties))
  }
  if (n_net > 0L) {
    warning(sprintf(
      "%d profile row(s) had no hospital-level LOS for their dominant class; the all-hospital LOS was used (los_level = 'network').",
      n_net
    ), call. = FALSE)
  }
  out
}

#' Filter Resistance Profiles to Classes with RR Estimates
#'
#' Drops profiles whose resistant classes have no RR estimate in
#' \code{rr_table} and re-normalises the remaining probabilities. The
#' all-susceptible reference profile is always kept.
#'
#' @param profiles_with_rr Named list returned by
#'   \code{daly_assign_rr_to_profiles()}.
#' @param rr_table Data frame of RR estimates.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param probability_col Character. Default \code{"probability"}.
#' @param fallback_rr Numeric. Default \code{1}.
#'
#' @return Named list (one entry per pathogen) with filtered and
#'   re-normalised profiles.
#' @export
daly_filter_profiles_to_rr_classes <- function(
  profiles_with_rr,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  probability_col = "probability",
  fallback_rr = 1
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by daly_assign_rr_to_profiles().")
  }
  if (!all(c(pathogen_col, class_col) %in% names(rr_table))) {
    stop(sprintf(
      "Columns '%s' and '%s' must be present in rr_table.",
      pathogen_col, class_col
    ))
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    # Classes that have a real RR estimate for this pathogen
    rr_classes <- rr_table[rr_table[[pathogen_col]] == path, class_col]
    rr_classes <- unique(rr_classes[!is.na(rr_classes)])

    # Identify class indicator columns present in this profile data frame
    non_class_cols <- c(
      "profile", probability_col, "RR_LOS_profile",
      "dominant_class", "CI_lower_profile",
      "CI_upper_profile", "numerator", "PAF_LOS",
      "denominator"
    )
    class_cols <- setdiff(names(df), non_class_cols)

    # Matchable classes: those present as indicator columns AND in rr_table
    matchable <- intersect(class_cols, rr_classes)

    if (length(matchable) == 0L) {
      warning(sprintf(
        "'%s': no profile classes overlap with relative LOS table -- all profiles dropped.",
        path
      ))
      next
    }

    # Keep profiles where:
    #   (a) the profile is all-susceptible (reference, always kept), OR
    #   (b) at least one matchable class is resistant (== 1)
    # This ensures the all-susceptible reference profile is never dropped.
    all_class_cols <- setdiff(names(df), non_class_cols)
    is_all_s <- if (length(all_class_cols) > 0L) {
      rowSums(df[, all_class_cols, drop = FALSE] == 1L, na.rm = TRUE) == 0L
    } else {
      rep(TRUE, nrow(df))
    }

    if (length(matchable) == 1L) {
      has_matchable_r <- df[[matchable]] == 1L
    } else {
      has_matchable_r <- rowSums(df[, matchable, drop = FALSE] == 1L) > 0L
    }

    keep <- is_all_s | has_matchable_r

    n_before <- nrow(df)
    df_kept <- df[keep, , drop = FALSE]
    n_after <- nrow(df_kept)
    n_dropped <- n_before - n_after

    if (n_after == 0L) {
      warning(sprintf("'%s': all profiles dropped after matching -- skipping.", path))
      next
    }

    # Re-normalise probabilities
    prob_sum <- sum(df_kept[[probability_col]], na.rm = TRUE)
    if (prob_sum <= 0) {
      warning(sprintf("'%s': probability sum is zero after filtering -- skipping.", path))
      next
    }
    df_kept[[probability_col]] <- df_kept[[probability_col]] / prob_sum

    out[[path]] <- df_kept

    message(sprintf(
      "'%s': %d -> %d profiles kept (%d dropped, no matchable relative LOS class). Prob re-normalised.",
      path, n_before, n_after, n_dropped
    ))
  }

  return(out)
}

#' Calculate Fatal Resistance Prevalence (R_k)
#'
#' Computes the profile-weighted expected proportion of deaths attributable to
#' resistance (R_k) and the expected odds ratio of death (E[OR_death]) for
#' each pathogen.
#'
#' @param profiles_with_rr Named list returned by
#'   \code{daly_assign_rr_to_profiles()} or
#'   \code{daly_filter_profiles_to_rr_classes()}.
#' @param probability_col Character. Profile probability column.
#'   Default \code{"probability"}.
#' @param rr_profile_col Character. Profile-level OR column.
#'   Default \code{"RR_LOS_profile"}.
#'
#' @return Named list (one entry per pathogen) each containing
#'   \code{per_profile} data frame, \code{R_k}, and \code{E_OR_k}.
#' @export
daly_calc_resistance_prevalence_fatal <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by daly_assign_rr_to_profiles().")
  }

  # Columns that are metadata / computed -- NOT class indicator columns
  non_class_cols <- c(
    "profile", probability_col, rr_profile_col,
    "dominant_class", "CI_lower_profile", "CI_upper_profile",
    "numerator", "PAF_LOS", "denominator",
    "numerator_assoc", "fraction_assoc",
    "numerator_R_kd", "R_kd"
  )

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_K_delta  (sums to 1)
    or <- df[[rr_profile_col]] # OR_death_K_delta (1.0 for all-susceptible)

    # E[OR_k] = sum_delta  R'_K_delta * OR_K_delta
    E_OR_k <- sum(p * or, na.rm = TRUE)

    if (!is.finite(E_OR_k) || E_OR_k <= 0) {
      warning(sprintf(
        "'%s': E[OR_death] = %.6g (must be > 0) -- all mortality ORs may be <= 1 or NA; skipping.",
        path, E_OR_k
      ))
      next
    }

    numerator_R_kd <- p * or
    R_kd_vec <- numerator_R_kd / E_OR_k

    df$numerator_R_kd <- round(numerator_R_kd, 6L)
    df$R_kd <- round(R_kd_vec, 6L)

    # Resistant profiles: at least one binary class indicator column == 1.
    # The all-susceptible profile has every class column = 0 and OR = 1.
    class_cols <- setdiff(names(df), non_class_cols)
    if (length(class_cols) > 0L) {
      is_resistant <- rowSums(
        df[, class_cols, drop = FALSE] == 1L,
        na.rm = TRUE
      ) > 0L
    } else {
      is_resistant <- or > 1.0
    }

    R_k <- sum(R_kd_vec[is_resistant], na.rm = TRUE)

    out[[path]] <- list(
      per_profile = df,
      R_k         = round(R_k, 6L),
      E_OR_k      = round(E_OR_k, 6L)
    )

    message(sprintf(
      "'%s': R_k (fatal resistance prevalence) = %.4f | E[OR_death] = %.4f | %d profiles (%d resistant).",
      path, R_k, E_OR_k, nrow(df), sum(is_resistant)
    ))
  }

  return(out)
}
