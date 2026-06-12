# daly_rr_and_los.R
# LOS relative risk, mortality OR, and resistance-profile fitting functions

.compute_mean_from_fit <- function(fit, dist) {
  tryCatch(
    {
      if (dist == "weibull") {
        k <- fit$estimate["shape"]
        l <- fit$estimate["scale"]
        l * gamma(1 + 1 / k)
      } else if (dist == "lnorm") {
        mu <- fit$estimate["meanlog"]
        sg <- fit$estimate["sdlog"]
        exp(mu + sg^2 / 2)
      } else if (dist == "gamma") {
        a <- fit$estimate["shape"]
        r <- fit$estimate["rate"]
        a / r
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
}

#' @keywords internal
.safe_fitdist <- function(x, dist) safe_fit(x, dist)


# -- Step 1a -------------------------------------------------------------------

#' Derive Infection Type (HAI / CAI) per Patient
#'
#' Classifies each row as HAI or CAI. Uses \code{infection_type_col} when it
#' contains a valid value. For rows where that column is \code{NA},
#' \code{"Not known"}, or \code{"NULL"}, derives the classification from the
#' gap between \code{date_culture_col} and \code{date_admission_col}:
#' \itemize{
#'   \item gap <= \code{hai_threshold_hours} -> \strong{CAI}
#'   \item gap >  \code{hai_threshold_hours} -> \strong{HAI}
#' }
#'
#' @param data Data frame.
#' @param infection_type_col Character. Raw infection type column.
#'   Default \code{"type_of_infection"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Date of first positive culture.
#'   Default \code{"date_of_first_positive_culture"}.
#' @param hai_threshold_hours Numeric. Gap threshold in hours. Default \code{48}.
#' @param patient_id_col Character. Unique patient identifier column.
#'   Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with column \code{infection_type_derived}
#'   (\code{"HAI"} / \code{"CAI"} / \code{"Unknown"}).
#' @export
daly_derive_hai_cai_for_los <- function(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
) {
  required <- c(infection_type_col, date_admission_col, date_culture_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  # Identify patients with ambiguous infection type AND at least one missing date.
  # These patients cannot have HAI/CAI inferred and are assigned "Not Known".
  ambiguous_mask <- {
    inf_raw_check <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[infection_type_col]])
    ))
    inf_raw_check %in% c("NOT KNOWN", "NOT_KNOWN", "UNKNOWN", "NULL", "NA", "") |
      is.na(data[[infection_type_col]])
  }
  missing_admit <- is.na(data[[date_admission_col]])
  missing_culture <- is.na(data[[date_culture_col]])
  cannot_infer <- ambiguous_mask & (missing_admit | missing_culture)

  if (any(cannot_infer)) {
    n_cannot <- sum(cannot_infer)
    message(sprintf(
      paste0(
        "Cannot infer infection type for %d patient(s): ",
        "infection type is unknown/null AND at least one date is missing. ",
        "Assigning 'Not Known'."
      ),
      n_cannot
    ))
    if (patient_id_col %in% names(data)) {
      flagged <- data[cannot_infer, patient_id_col, drop = TRUE]
      message("  Patient IDs with missing date(s):")
      message(paste0("    ", paste(flagged, collapse = ", ")))
    } else {
      message(
        "  Row indices with missing date(s): ",
        paste(which(cannot_infer), collapse = ", ")
      )
    }
    n_miss_admit <- sum(ambiguous_mask & missing_admit)
    n_miss_culture <- sum(ambiguous_mask & missing_culture)
    message(sprintf(
      "  Breakdown: missing %s = %d | missing %s = %d",
      date_admission_col, n_miss_admit,
      date_culture_col,   n_miss_culture
    ))
  }

  data <- data %>%
    dplyr::mutate(
      .inf_raw = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[infection_type_col]])
      )),
      .gap_h = as.numeric(difftime(
        as.Date(.data[[date_culture_col]]),
        as.Date(.data[[date_admission_col]]),
        units = "hours"
      )),
      .cannot_infer = .inf_raw %in% c(
        "NOT KNOWN", "NOT_KNOWN", "UNKNOWN",
        "NULL", "NA", ""
      ) |
        is.na(.data[[infection_type_col]]),
      infection_type_derived = dplyr::case_when(
        .inf_raw %in% c(
          "HAI", "HOSPITAL ACQUIRED",
          "HOSPITAL-ACQUIRED", "HOSPITAL_ACQUIRED"
        ) |
          grepl("HOSPITAL.ACQUIRED|HEALTH.CARE.ASSOCIATED|HEALTHCARE.ASSOCIATED|\\bHAI\\b",
            .inf_raw,
            perl = TRUE
          ) ~ "HAI",
        .inf_raw %in% c(
          "CAI", "COMMUNITY ACQUIRED",
          "COMMUNITY-ACQUIRED", "COMMUNITY_ACQUIRED"
        ) |
          grepl("COMMUNITY.ACQUIRED|\\bCAI\\b",
            .inf_raw,
            perl = TRUE
          ) ~ "CAI",
        # Ambiguous type (Not Known/NA) but BOTH dates present: infer from gap
        .cannot_infer &
          !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours, "CAI", "HAI"),
        # Ambiguous type and at least one date missing: cannot infer
        .cannot_infer ~ "Not Known",
        # Unrecognised label but BOTH dates present: infer from gap
        !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours, "CAI", "HAI"),
        # Unrecognised label and dates missing: cannot resolve
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::select(-".inf_raw", -".gap_h", -".cannot_infer")

  n_hai <- sum(data$infection_type_derived == "HAI", na.rm = TRUE)
  n_cai <- sum(data$infection_type_derived == "CAI", na.rm = TRUE)
  n_not_known <- sum(data$infection_type_derived == "Not Known", na.rm = TRUE)
  n_unknown <- sum(data$infection_type_derived == "Unknown", na.rm = TRUE)
  message(sprintf(
    "Infection type derived: HAI = %d | CAI = %d | Not Known = %d | Unknown = %d.",
    n_hai, n_cai, n_not_known, n_unknown
  ))
  return(data)
}


# -- Step 1b -------------------------------------------------------------------

#' Compute Patient-Level Post-Infection LOS
#'
#' Computes LOS with infection-type-specific clock start:
#' \itemize{
#'   \item \strong{CAI}: LOS = date_discharge - date_admission
#'   \item \strong{HAI}: LOS = date_discharge - date_culture
#'   \item \strong{Unknown / Not known}: LOS = date_discharge - date_culture
#'     (if culture date present), else date_discharge - date_admission. If
#'     neither reference date is available the patient is excluded.
#' }
#' Only discharged patients are retained. A checkpoint reports how many are
#' missing a discharge date; those patients fall back to \code{los_col} when
#' provided, otherwise they are excluded. Rows with LOS <= 0 or
#' LOS > \code{max_los} are dropped. Returns one row per patient.
#'
#' @param data Data frame (after \code{daly_derive_hai_cai_for_los()} has been run).
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param organism_col Character. Column containing organism/pathogen names.
#'   Used for episode grouping: a new episode is only created when the culture
#'   date gap exceeds 14 days AND the organism differs from the episode's first
#'   organism. Default \code{"organism_name"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_discharge_col Character. Default \code{"final_outcome_date"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_value Character. Default \code{"Discharged"}.
#' @param infection_type_derived_col Character. Column from
#'   \code{daly_derive_hai_cai_for_los()}. Default \code{"infection_type_derived"}.
#' @param syndrome_col Character. Syndrome column name. Only used when
#'   \code{syndrome_name} is not \code{NULL}. Default \code{"syndrome"}.
#' @param syndrome_name Character or \code{NULL}. If provided, only patients
#'   with this syndrome are retained before LOS computation. Default \code{NULL}.
#' @param los_col Character or \code{NULL}. Optional name of a pre-computed
#'   LOS column (e.g. recorded directly in the data). Used only for patients
#'   whose discharge date is missing. Values that are \code{NA}, zero, or
#'   negative are treated as invalid and those patients are excluded.
#'   Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. If provided, filters data
#'   to the specified facility before computing LOS. Default \code{NULL}.
#' @param max_los Numeric. Upper cap on LOS in days; patients exceeding this
#'   are excluded. Default \code{200}.
#'
#' @return Data frame: one row per patient with \code{patient_id_col},
#'   \code{facility_col}, \code{infection_type_derived_col}, \code{LOS_days}.
#' @export
daly_compute_patient_los <- function(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  syndrome_name = NULL,
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  infection_type_derived_col = "infection_type_derived",
  los_col = NULL,
  max_los = 200
) {
  # -- Column validation ------------------------------------------------------
  required <- c(
    patient_id_col, facility_col, organism_col,
    date_admission_col, date_discharge_col, date_culture_col,
    final_outcome_col, infection_type_derived_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  if (!is.null(los_col) && !los_col %in% names(data)) {
    stop(sprintf("los_col '%s' not found in data.", los_col))
  }

  # -- Facility filtering -----------------------------------------------------
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

  # -- Syndrome filtering ----------------------------------------------------
  if (!is.null(syndrome_name)) {
    avail_syn <- unique(data[[syndrome_col]])
    if (!syndrome_name %in% avail_syn) {
      stop(sprintf(
        "syndrome_name '%s' not found in column '%s'. Available: %s",
        syndrome_name, syndrome_col,
        paste(sort(as.character(avail_syn)), collapse = ", ")
      ))
    }
    data <- data[data[[syndrome_col]] == syndrome_name, ]
    message(sprintf(
      "Filtering to syndrome '%s': %d row(s) retained.",
      syndrome_name, nrow(data)
    ))
  }

  # -- Checkpoint: discharge date coverage -----------------------------------
  discharged_data <- data[data[[final_outcome_col]] == final_outcome_value, ]
  n_total <- nrow(discharged_data)
  n_no_disc <- sum(is.na(as.Date(discharged_data[[date_discharge_col]])))

  if (n_no_disc > 0L) {
    message(sprintf(
      "Checkpoint: %d of %d discharged patient(s) are missing a discharge date.",
      n_no_disc, n_total
    ))
    if (!is.null(los_col)) {
      message(sprintf(
        "  -> Will fall back to '%s' column for those patients (values <= 0 or NA will be excluded).",
        los_col
      ))
    } else {
      message(
        "  -> No los_col provided; patients with missing discharge date will be excluded."
      )
    }
  }

  # -- Pre-compute episode-minimum culture dates ------------------------------
  # A new episode starts only when BOTH conditions are met:
  #   1. Gap from the current episode's start date > 14 days
  #   2. Organism differs from the current episode's first organism
  # If the same organism returns after > 14 days it is treated as the same
  # infection event; the original (first) culture date is used for LOS.
  episode_min_cult <- data %>%
    dplyr::filter(
      .data[[final_outcome_col]] == final_outcome_value,
      !is.na(.data[[date_culture_col]])
    ) %>%
    dplyr::distinct(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(date_culture_col),
      !!rlang::sym(organism_col)
    ) %>%
    dplyr::mutate(.raw_cult = as.Date(.data[[date_culture_col]])) %>%
    dplyr::arrange(!!rlang::sym(patient_id_col), .raw_cult) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::mutate(
      .ep_id = {
        n <- dplyr::n()
        ep <- integer(n)
        org_vec <- .data[[organism_col]]
        ep_start <- .raw_cult[1L]
        ep_org <- org_vec[1L]
        current_ep <- 0L
        for (i in seq_len(n)) {
          if (i > 1L &&
            as.numeric(.raw_cult[i] - ep_start) > 14L &&
            !identical(org_vec[i], ep_org)) {
            current_ep <- current_ep + 1L
            ep_start <- .raw_cult[i]
            ep_org <- org_vec[i]
          }
          ep[i] <- current_ep
        }
        ep
      }
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col), .ep_id) %>%
    dplyr::mutate(.cult_ep_min = min(.raw_cult)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(date_culture_col),
      !!rlang::sym(organism_col),
      .cult_ep_min
    )

  # -- LOS computation -------------------------------------------------------
  unknown_types_lc <- c("unknown", "not known")

  los_data <- data %>%
    dplyr::filter(.data[[final_outcome_col]] == final_outcome_value) %>%
    dplyr::left_join(episode_min_cult,
      by = c(patient_id_col, date_culture_col, organism_col)
    ) %>%
    dplyr::mutate(
      .adm = as.Date(.data[[date_admission_col]]),
      .disc = as.Date(.data[[date_discharge_col]]),
      .cult = .cult_ep_min,
      .inf_lc = tolower(trimws(as.character(.data[[infection_type_derived_col]]))),
      LOS_days = dplyr::case_when(
        # -- Patients WITH a discharge date ----------------------------
        # HAI: clock starts at culture date
        !is.na(.disc) & .inf_lc == "hai" ~
          as.numeric(.disc - .cult),
        # CAI: clock starts at admission date
        !is.na(.disc) & .inf_lc == "cai" ~
          as.numeric(.disc - .adm),
        # Unknown / Not known: prefer culture date, fall back to admission date
        !is.na(.disc) & .inf_lc %in% unknown_types_lc & !is.na(.cult) ~
          as.numeric(.disc - .cult),
        !is.na(.disc) & .inf_lc %in% unknown_types_lc & is.na(.cult) & !is.na(.adm) ~
          as.numeric(.disc - .adm),
        # Unknown / Not known with no reference date: exclude
        !is.na(.disc) & .inf_lc %in% unknown_types_lc ~
          NA_real_,
        # All other infection types with discharge date: fallback to admission
        !is.na(.disc) ~
          as.numeric(.disc - .adm),
        # -- Patients WITHOUT a discharge date: resolved below ---------
        TRUE ~ NA_real_
      )
    )

  # For patients missing a discharge date, fall back to los_col if provided.
  # Values <= 0 or NA in los_col are treated as invalid and excluded.
  if (!is.null(los_col)) {
    los_data <- los_data %>%
      dplyr::mutate(
        LOS_days = dplyr::if_else(
          is.na(.disc),
          dplyr::if_else(
            !is.na(.data[[los_col]]) & as.numeric(.data[[los_col]]) > 0,
            as.numeric(.data[[los_col]]),
            NA_real_
          ),
          LOS_days
        )
      )
  }

  # -- Final filtering and deduplication -------------------------------------
  los_data <- los_data %>%
    dplyr::select(-".adm", -".disc", -".cult", -".cult_ep_min", -".inf_lc") %>%
    dplyr::filter(!is.na(LOS_days), LOS_days > 0L, LOS_days <= max_los) %>%
    dplyr::distinct(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(facility_col),
      !!rlang::sym(infection_type_derived_col),
      LOS_days
    )

  message(sprintf(
    "Patient LOS computed: %d patients retained (LOS > 0 and <= %d days).",
    nrow(los_data), max_los
  ))
  return(los_data)
}


# -- Step 1c (internal) --------------------------------------------------------

#' Collapse drug-level data to antibiotic-class binary wide matrix
#'
#' Standardises antibiotic values (Intermediate -> R/S per GBD rules),
#' collapses to class level (class = 1 if ANY drug in class = R), and pivots
#' wide: one row per patient, one column per antibiotic class (0/1).
#' Column names are sanitised with make.names(); the original-to-safe name
#' mapping is stored as the "class_name_map" attribute.
#'
#' @keywords internal
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


# -- daly_fit_los_rr -----------------------------------------------------------

#' Fit relative LOS using Gamma GLM with log link
#'
#' Estimates class-specific relative length of stay (RR_LOS) for each pathogen
#' using a Gamma generalized linear model with log link:
#'
#'   E(LOS_i) = exp(beta0 + beta1*Resistance_ci + beta2*HAI_i +
#'                  beta3*Age_i + beta4*Sex_i + beta5*ICU_i +
#'                  beta6*Comorbidity_i + centre fixed effects)
#'
#' RR_LOS for a resistance class is exp(beta1).
#'
#'
#' @param data Data frame.
#' @param patient_id_col Character. Patient identifier column.
#' @param facility_col Character or NULL. Facility / centre column. When
#'   \code{NULL}, no hospital effect is included in the model and a single
#'   pooled RR is returned (hospital column will be \code{NA}). When provided,
#'   the model is fit separately per hospital and the result contains one row
#'   per pathogen-class-hospital combination (hospital-specific RR).
#' @param organism_col Character. Pathogen column.
#' @param syndrome_col Character. Syndrome column.
#' @param infection_type_col Character. Raw infection type column.
#' @param antibiotic_class_col Character. Antibiotic class column.
#' @param antibiotic_name_col Character. Antibiotic name column.
#' @param antibiotic_value_col Character. Antibiotic susceptibility column.
#' @param date_admission_col Character. Admission date column.
#' @param date_discharge_col Character. Discharge / final outcome date column.
#' @param date_culture_col Character. First positive culture date column.
#' @param final_outcome_col Character. Final outcome column.
#' @param final_outcome_value Character. Value indicating discharged patients.
#' @param age_col Character. Age column.
#' @param sex_col Character. Sex column.
#' @param unit_type_col Character or NULL. ICU / ward location column.
#' @param comorbidity_col Character or NULL. Comorbidity column.
#' @param syndrome_name Character or NULL. Restrict to one syndrome.
#' @param organism_name Character vector or NULL. Restrict to these pathogens.
#' @param hai_threshold_hours Numeric. HAI derivation threshold.
#' @param max_los Numeric. Maximum retained LOS in days.
#' @param min_n Integer. Minimum patients required in a class model.
#' @param min_resistant Integer. Minimum resistant patients required.
#' @param min_susceptible Integer. Minimum susceptible patients required.
#' @param add_centre_fe Logical. Retained for backward compatibility; ignored
#'   when \code{facility_col} is provided (per-hospital fitting is used) or
#'   \code{NULL} (no facility information available).
#'
#' @return Data frame with pathogen, antibiotic_class, hospital (NA when
#'   \code{facility_col} is \code{NULL}), RR_LOS, CI_lower,
#'   CI_upper, beta_log_rr, mean_los_resistant, mean_los_susceptible,
#'   n_patients, n_resistant, n_susceptible, model_family, link_function,
#'   and syndrome_scope.
#' @export
daly_fit_los_rr <- function(
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
  min_susceptible = 5L,
  add_centre_fe = TRUE
) {
  required_cols <- c(
    patient_id_col, organism_col,
    antibiotic_class_col, antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col, date_culture_col,
    final_outcome_col
  )

  if (!is.null(facility_col)) {
    required_cols <- c(required_cols, facility_col)
  }
  if (!is.null(age_col)) {
    required_cols <- c(required_cols, age_col)
  }
  if (!is.null(sex_col)) {
    required_cols <- c(required_cols, sex_col)
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

  # When facility_col is NULL, no hospital random/fixed effect is used.
  # A dummy column is added internally so helper functions remain compatible.
  use_facility <- !is.null(facility_col)
  if (!use_facility) {
    facility_col <- ".facility_pooled"
    data[[facility_col]] <- "pooled"
    message("No facility_col provided: hospital effect will not be included. Returning pooled RR.")
  }

  # Step 1: derive HAI / CAI
  data <- daly_derive_hai_cai_for_los(
    data = data,
    infection_type_col = infection_type_col,
    date_admission_col = date_admission_col,
    date_culture_col = date_culture_col,
    hai_threshold_hours = hai_threshold_hours,
    patient_id_col = patient_id_col
  )

  # Step 2: compute patient-level LOS
  # HAI: LOS = discharge date - first positive culture date
  # CAI: LOS = discharge date - admission date
  los_df <- daly_compute_patient_los(
    data = data,
    patient_id_col = patient_id_col,
    facility_col = facility_col,
    organism_col = organism_col,
    syndrome_col = syndrome_col,
    syndrome_name = syndrome_name,
    date_admission_col = date_admission_col,
    date_discharge_col = date_discharge_col,
    date_culture_col = date_culture_col,
    final_outcome_col = final_outcome_col,
    final_outcome_value = final_outcome_value,
    infection_type_derived_col = "infection_type_derived",
    max_los = max_los
  )

  if (nrow(los_df) == 0L) {
    warning("No valid LOS observations found after filtering.")
    return(data.frame())
  }

  # Step 3: build class-level resistance wide matrix
  # class = 1 if any drug in class is resistant
  # class = 0 if tested and no drug in class is resistant
  # class = NA if class not tested
  resistance_wide <- .build_class_resistance_wide(
    data = data,
    patient_id_col = patient_id_col,
    antibiotic_class_col = antibiotic_class_col,
    antibiotic_name_col = antibiotic_name_col,
    antibiotic_value_col = antibiotic_value_col,
    untested_fill = NA_integer_
  )

  class_name_map <- attr(resistance_wide, "class_name_map")
  class_safe_cols <- setdiff(names(resistance_wide), patient_id_col)

  # Step 4: patient-level covariates
  patient_covars <- data |>
    dplyr::group_by(.data[[patient_id_col]]) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      !!patient_id_col := .data[[patient_id_col]],
      !!facility_col := .data[[facility_col]],
      !!organism_col := .data[[organism_col]],
      infection_type_derived = .data[["infection_type_derived"]],
      HAI = dplyr::if_else(.data[["infection_type_derived"]] == "HAI", 1L, 0L),
      Age_model = if (!is.null(age_col) && age_col %in% names(data)) {
        suppressWarnings(as.numeric(.data[[age_col]]))
      } else {
        NA_real_
      },
      Sex_model = if (!is.null(sex_col) && sex_col %in% names(data)) {
        as.factor(.data[[sex_col]])
      } else {
        factor(NA_character_)
      }
    )

  # ICU covariate
  if (!is.null(unit_type_col) && unit_type_col %in% names(data)) {
    icu_tbl <- .derive_icu_binary(
      data = data,
      patient_id_col = patient_id_col,
      unit_type_col = unit_type_col
    )
    patient_covars <- dplyr::left_join(patient_covars, icu_tbl, by = patient_id_col)
  } else {
    patient_covars$ICU <- NA_integer_
  }

  # Comorbidity covariate
  if (!is.null(comorbidity_col) && comorbidity_col %in% names(data)) {
    comorb_df <- data |>
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

  # Step 5: merge LOS + covariates + resistance
  model_base <- los_df |>
    dplyr::inner_join(
      patient_covars,
      by = c(patient_id_col, facility_col, "infection_type_derived")
    ) |>
    dplyr::inner_join(resistance_wide, by = patient_id_col)

  if (nrow(model_base) == 0L) {
    warning("No observations available after joining LOS, covariates, and resistance matrix.")
    return(data.frame())
  }

  # Step 6: resolve pathogen list
  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(model_base[[organism_col]]))
  }

  results <- list()

  # Step 7: fit Gamma(log) model per pathogen x class
  for (path in pathogens) {
    path_df <- model_base |>
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
          stringr::str_to_lower(stringr::str_trim(path))
      )

    if (nrow(path_df) == 0L) {
      next
    }

    for (cls in class_safe_cols) {
      orig_class <- class_name_map[[cls]]

      sub <- path_df |>
        dplyr::filter(!is.na(.data[[cls]]))

      if (nrow(sub) == 0L) {
        next
      }

      # When facility_col is provided, fit a hospital-specific model per hospital.
      # When facility_col is NULL (pooled), fit a single model across all patients.
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

        n_resistant <- sum(hosp_sub[[cls]] == 1L, na.rm = TRUE)
        n_susceptible <- sum(hosp_sub[[cls]] == 0L, na.rm = TRUE)

        if (n_resistant < min_resistant || n_susceptible < min_susceptible) {
          next
        }

        if (length(unique(hosp_sub[[cls]])) < 2L) {
          next
        }

        if (all(is.na(hosp_sub$LOS_days)) || any(hosp_sub$LOS_days <= 0, na.rm = TRUE)) {
          next
        }

        fixed_terms <- c(sprintf("`%s`", cls), "HAI")

        if ("Age_model" %in% names(hosp_sub) && sum(!is.na(hosp_sub$Age_model)) > 0L) {
          fixed_terms <- c(fixed_terms, "Age_model")
        }

        if ("Sex_model" %in% names(hosp_sub) &&
          length(unique(stats::na.omit(as.character(hosp_sub$Sex_model)))) > 1L) {
          fixed_terms <- c(fixed_terms, "Sex_model")
        }

        if ("ICU" %in% names(hosp_sub) &&
          length(unique(stats::na.omit(as.character(hosp_sub$ICU)))) > 1L) {
          fixed_terms <- c(fixed_terms, "ICU")
        }

        if ("comorbidity_encoded" %in% names(hosp_sub) &&
          length(unique(stats::na.omit(as.character(hosp_sub$comorbidity_encoded)))) > 1L) {
          fixed_terms <- c(fixed_terms, "comorbidity_encoded")
        }

        # No centre fixed effect: when use_facility is TRUE the model is already
        # restricted to one hospital; when use_facility is FALSE there is no
        # hospital information to adjust for.

        fmla <- stats::as.formula(
          paste("LOS_days ~", paste(fixed_terms, collapse = " + "))
        )

        fit <- tryCatch(
          stats::glm(
            formula = fmla,
            data = hosp_sub,
            family = stats::Gamma(link = "log")
          ),
          error = function(e) NULL
        )

        if (is.null(fit)) {
          next
        }

        coefs <- stats::coef(fit)
        vcov_mat <- stats::vcov(fit)

        cls_pattern <- paste0("^`?", gsub(".", "\\\\.", cls, fixed = TRUE), "`?$")
        coef_name <- grep(cls_pattern, names(coefs), value = TRUE)

        if (length(coef_name) == 0L) {
          next
        }

        beta <- coefs[[coef_name[1L]]]
        se <- sqrt(vcov_mat[coef_name[1L], coef_name[1L]])

        rr_los <- exp(beta)
        ci_lower <- exp(beta - 1.96 * se)
        ci_upper <- exp(beta + 1.96 * se)

        mean_los_resistant <- mean(hosp_sub$LOS_days[hosp_sub[[cls]] == 1L], na.rm = TRUE)
        mean_los_susceptible <- mean(hosp_sub$LOS_days[hosp_sub[[cls]] == 0L], na.rm = TRUE)

        results[[length(results) + 1L]] <- data.frame(
          pathogen = path,
          antibiotic_class = orig_class,
          hospital = if (use_facility) hosp else NA_character_,
          beta_log_rr = unname(beta),
          RR_LOS = unname(rr_los),
          CI_lower = unname(ci_lower),
          CI_upper = unname(ci_upper),
          mean_los_resistant = mean_los_resistant,
          mean_los_susceptible = mean_los_susceptible,
          n_patients = nrow(hosp_sub),
          n_resistant = n_resistant,
          n_susceptible = n_susceptible,
          model_family = "Gamma",
          link_function = "log",
          syndrome_scope = if (!is.null(syndrome_name)) syndrome_name else "all",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0L) {
    warning("No Gamma GLM LOS models could be fitted.")
    return(data.frame())
  }

  out <- dplyr::bind_rows(results)

  message(sprintf(
    "Gamma(log) LOS models fitted: %d pathogen-class estimates returned.",
    nrow(out)
  ))

  out
}


# -- Step 2b -------------------------------------------------------------------

#' Estimate Per-Class LOS Relative Risk via Parametric Distribution Fitting
#'
#' For each pathogen x antibiotic class, splits patients into resistant (R) and
#' susceptible (S) groups, fits a parametric distribution (default: gamma) to
#' each group's LOS, and returns:
#'   RR_LOS = mean_LOS(R) / mean_LOS(S)
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
#' @param infection_type_col Character. Column used to derive HAI/CAI.
#' @param antibiotic_class_col Character. Antibiotic class column.
#' @param antibiotic_name_col Character. Antibiotic name column.
#' @param antibiotic_value_col Character. Antibiotic value column (S/I/R).
#' @param date_admission_col Character. Date of admission column.
#' @param date_discharge_col Character. Date of discharge column.
#' @param date_culture_col Character. Date of culture column.
#' @param final_outcome_col Character. Final outcome column.
#' @param final_outcome_value Character. Value indicating discharge.
#' @param syndrome_name Character or \code{NULL}. Restrict to this syndrome.
#' @param organism_name Character vector or \code{NULL}. Restrict to these
#'   pathogen(s); otherwise all pathogens in the filtered data.
#' @param facility_name Character or \code{NULL}. If provided, filters data
#'   to the specified facility before fitting. Default \code{NULL}.
#' @param hai_threshold_hours Numeric. Hours after admission before a culture
#'   is classified as HAI. Default \code{48}.
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
#'     \item{RR_LOS}{Fitted mean ratio: mean_LOS(R) / mean_LOS(S).}
#'     \item{mean_LOS_R}{Fitted mean LOS for resistant patients.}
#'     \item{mean_LOS_S}{Fitted mean LOS for susceptible patients.}
#'     \item{best_dist_R}{Distribution fitted to the R group.}
#'     \item{best_dist_S}{Distribution fitted to the S group.}
#'     \item{n_R}{Number of resistant patients.}
#'     \item{n_S}{Number of susceptible patients.}
#'     \item{syndrome_scope}{Syndrome filter applied, or "all".}
#'   }
#'
#' @seealso \code{\link{daly_fit_los_rr}}, \code{daly_assign_rr_to_profiles}
#' @export
daly_fit_los_rr_distribution <- function(
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
) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required: install.packages('fitdistrplus')")
  }

  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col,
    date_culture_col, final_outcome_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # -- Derive infection type (HAI / CAI) --------------------------------------
  data <- daly_derive_hai_cai_for_los(
    data,
    infection_type_col  = infection_type_col,
    date_admission_col  = date_admission_col,
    date_culture_col    = date_culture_col,
    hai_threshold_hours = hai_threshold_hours
  )

  # -- Filter to discharged patients + optional syndrome + optional facility --
  df <- dplyr::filter(data, .data[[final_outcome_col]] == final_outcome_value)
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .data[[syndrome_col]] == syndrome_name)
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

  # -- Resolve pathogens ------------------------------------------------------
  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  # -- Inner helper: fit distribution(s) and return analytical mean -----------
  .fit_dist_mean <- function(los_vec) {
    fits <- Filter(
      Negate(is.null),
      setNames(
        lapply(distributions, .safe_fitdist, x = los_vec),
        distributions
      )
    )
    if (length(fits) == 0L) {
      return(list(mean = NA_real_, dist = NA_character_))
    }
    aics <- sapply(fits, `[[`, "aic")
    best <- names(which.min(aics))
    list(mean = .compute_mean_from_fit(fits[[best]], best), dist = best)
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

    # -- Patient LOS --------------------------------------------------------
    los_pat <- daly_compute_patient_los(
      path_df,
      patient_id_col             = patient_id_col,
      facility_col               = facility_col,
      facility_name              = NULL,
      date_admission_col         = date_admission_col,
      date_discharge_col         = date_discharge_col,
      date_culture_col           = date_culture_col,
      final_outcome_col          = final_outcome_col,
      final_outcome_value        = final_outcome_value,
      infection_type_derived_col = "infection_type_derived",
      max_los                    = max_los
    )
    if (nrow(los_pat) == 0L) {
      message(sprintf("'%s': no patients with valid LOS, skipping.", path))
      next
    }

    # -- Class resistance wide matrix (patient x class, binary 0/1) ----------
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

    # -- Per-class R vs S distribution fit ----------------------------------
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

      est_r <- .fit_dist_mean(los_r)
      est_s <- .fit_dist_mean(los_s)

      if (is.na(est_r$mean) || is.na(est_s$mean) || est_s$mean <= 0) {
        message(sprintf(
          "'%s' | class '%s': distribution fitting failed (mean_R=%s, mean_S=%s), skipping.",
          path, orig_class,
          if (is.na(est_r$mean)) "NA" else round(est_r$mean, 3L),
          if (is.na(est_s$mean)) "NA" else round(est_s$mean, 3L)
        ))
        next
      }

      rr_los <- est_r$mean / est_s$mean

      if (rr_los > 10) {
        message(sprintf(
          "'%s' | class '%s': extreme RR_LOS = %.2f (mean_R=%.2f / mean_S=%.2f) -- check data.",
          path, orig_class, rr_los, est_r$mean, est_s$mean
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
        mean_LOS_R       = round(est_r$mean, 3L),
        mean_LOS_S       = round(est_s$mean, 3L),
        mean_LOS_R_years = round(est_r$mean / 365, 6L),
        mean_LOS_S_years = round(est_s$mean / 365, 6L),
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
      "'%s': %d class relative LOS(s) estimated (gamma distribution, syndrome=%s).",
      path, n_fitted,
      if (!is.null(syndrome_name)) syndrome_name else "all"
    ))
  }

  if (length(all_rr) == 0L) {
    warning("No relative LOS values could be estimated. Check data, filters, and min_n.")
    return(data.frame(
      pathogen = character(), antibiotic_class = character(),
      facility_name = character(), RR_LOS = numeric(),
      mean_LOS_R = numeric(), mean_LOS_S = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  return(dplyr::bind_rows(all_rr))
}


# -- HAI/CAI classification for mortality cohort --------------------------------

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
  # -- Column validation ------------------------------------------------------
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

  # -- Death-date data-quality check -----------------------------------------
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

  # -- Identify ambiguous labels ----------------------------------------------
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

  # -- Derive infection type --------------------------------------------------
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
        # -- Explicit HAI labels ----------------------------------------
        .inf_raw %in% c(
          "HAI",
          "HOSPITAL ACQUIRED",
          "HOSPITAL-ACQUIRED",
          "HOSPITAL_ACQUIRED",
          "NOSOCOMIAL",
          "HOSPITAL ACQUIRED INFECTION",
          "HOSPITAL-ACQUIRED INFECTION"
        ) ~ "HAI",
        # -- Explicit CAI labels ----------------------------------------
        .inf_raw %in% c(
          "CAI",
          "COMMUNITY ACQUIRED",
          "COMMUNITY-ACQUIRED",
          "COMMUNITY_ACQUIRED",
          "COMMUNITY ACQUIRED INFECTION",
          "COMMUNITY-ACQUIRED INFECTION"
        ) ~ "CAI",
        # -- Date-gap derivation (both dates present) ------------------
        .cannot_infer &
          !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours,
            "CAI", "HAI"
          ),
        # -- Cannot infer ----------------------------------------------
        .cannot_infer ~ "Not Known",
        # -- Unrecognised non-empty label ------------------------------
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


# -- Step M2 -------------------------------------------------------------------

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


# -- Step M3 -------------------------------------------------------------------

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


# -- Step M4 -------------------------------------------------------------------

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
  print(tbl)

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


# -- daly_fit_mortality_rr -----------------------------------------------------

# daly_rr_and_los.R
# Mortality model for adjusted OR and adjusted RR of death

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
#' @param use_random_intercept Logical. Retained for backward compatibility;
#'   ignored when \code{facility_col} is provided (per-hospital fitting is
#'   used) or \code{NULL} (no facility information available).
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
  min_deaths = 10L,
  use_random_intercept = TRUE
) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for daly_fit_mortality_rr().")
  }

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

  # ---------------------------------------------------------------------------
  # Step 1: derive HAI / CAI for mortality model
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Step 2: optional syndrome restriction
  # ---------------------------------------------------------------------------
  df <- data
  if (!is.null(syndrome_name)) {
    df <- df[df[[syndrome_col]] == syndrome_name, , drop = FALSE]
  }

  if (nrow(df) == 0L) {
    warning("No rows remain after syndrome filtering.")
    return(data.frame())
  }

  # ---------------------------------------------------------------------------
  # Step 3: ICU indicator
  # ---------------------------------------------------------------------------
  icu_tbl <- .derive_icu_binary(
    data = df,
    patient_id_col = patient_id_col,
    unit_type_col = unit_type_col,
    icu_values = icu_values
  )

  # ---------------------------------------------------------------------------
  # Step 4: patient-level covariates
  # one row per patient
  # ---------------------------------------------------------------------------
  patient_covars <- df |>
    dplyr::group_by(.data[[patient_id_col]]) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      !!patient_id_col := .data[[patient_id_col]],
      !!facility_col := .data[[facility_col]],
      !!organism_col := .data[[organism_col]],
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

  # ---------------------------------------------------------------------------
  # Step 5: comorbidity encoding
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # Step 6: build class-level resistance matrix
  # ---------------------------------------------------------------------------
  resistance_wide <- .build_class_resistance_wide(
    data = df,
    patient_id_col = patient_id_col,
    antibiotic_class_col = antibiotic_class_col,
    antibiotic_name_col = antibiotic_name_col,
    antibiotic_value_col = antibiotic_value_col,
    untested_fill = NA_integer_
  )

  class_name_map <- attr(resistance_wide, "class_name_map")
  class_safe_cols <- setdiff(names(resistance_wide), patient_id_col)

  # ---------------------------------------------------------------------------
  # Step 7: merge patient-level model data
  # ---------------------------------------------------------------------------
  model_base <- patient_covars |>
    dplyr::inner_join(resistance_wide, by = patient_id_col) |>
    dplyr::filter(!is.na(death), !is.na(HAI), !is.na(Age_model), !is.na(Sex_model))

  if (nrow(model_base) == 0L) {
    warning("No complete patient-level rows available for mortality modelling.")
    return(data.frame())
  }

  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(model_base[[organism_col]]))
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
    path_df <- model_base |>
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
          stringr::str_to_lower(stringr::str_trim(path))
      )

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
# los.R
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


#' Prepare LOS Dataset
#'
#' Calculates length of stay from admission and outcome dates, filters to
#' discharged patients with valid LOS within a plausible range.
#'
#' @param data Data frame with patient-level records.
#' @param admission_col Character. Column name for admission date.
#'   Default \code{"date_of_admission"}.
#' @param outcome_date_col Character. Column name for outcome/discharge date.
#'   Default \code{"final_outcome_date"}.
#' @param outcome_col Character. Column name for outcome status.
#'   Default \code{"final_outcome"}.
#' @param patient_id_col Character. Column name for patient identifier.
#'   Default \code{"PatientInformation_id"}.
#' @param max_los Numeric. Maximum plausible LOS in days. Default 200.
#'
#' @return Data frame with one row per patient-admission, including
#'   \code{LOS_days}.
#' @export
prepare_los_data <- function(data,
                             admission_col = "date_of_admission",
                             outcome_date_col = "final_outcome_date",
                             outcome_col = "final_outcome",
                             patient_id_col = "PatientInformation_id",
                             max_los = 200) {
  data %>%
    dplyr::mutate(
      .adm_date = as.Date(.data[[admission_col]]),
      .out_date = as.Date(.data[[outcome_date_col]]),
      LOS_days  = as.numeric(.out_date - .adm_date)
    ) %>%
    dplyr::filter(
      .data[[outcome_col]] == "Discharged",
      !is.na(.data$LOS_days),
      .data$LOS_days > 0,
      .data$LOS_days <= max_los
    ) %>%
    dplyr::distinct(
      .data[[patient_id_col]],
      .data$.adm_date,
      .data$.out_date,
      .data$LOS_days
    ) %>%
    dplyr::select(-".adm_date", -".out_date")
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
#' @param los_data Data frame from \code{prepare_los_data()}.
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
    dplyr::filter(.data[[organism_col]] == organism)

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

  los_R <- org_los %>%
    dplyr::filter(.data$R_count > 0) %>%
    tidyr::uncount(.data$R_count) %>%
    dplyr::pull("LOS_days")

  los_S <- org_los %>%
    dplyr::filter(.data$S_count > 0) %>%
    tidyr::uncount(.data$S_count) %>%
    dplyr::pull("LOS_days")

  list(R = los_R, S = los_S)
}


# ===========================================================================
# Relative risk reference data, pathogen-drug mappings, and profile-level
# RR assignment.
#
# Loads GBD-derived relative risk (RR) estimates for pathogen-antibiotic class
# combinations, maps organism and drug class names to RR categories, and
# assigns profile-level RR values to resistance profiles using the GBD max
# rule: the profile RR equals the highest class-level RR among all resistant
# classes in that profile. Re-normalises profile probabilities after dropping
# profiles whose resistant classes carry no RR estimate, and computes the
# fatal resistance prevalence R_k and expected odds ratio E[OR_death] for
# each pathogen.
#
#   daly_load_rr_reference()              loads rr_list_gbd.xlsx
#   daly_add_rr_mappings()                appends rr_pathogen and rr_drug columns
#   daly_assign_rr_to_profiles()          assigns RR_LOS_profile via max rule
#   daly_filter_profiles_to_rr_classes()  drops and re-normalises unmatched profiles
#   daly_calc_resistance_prevalence_fatal() computes R_k and E[OR_death]
# ===========================================================================

#' Load RR (Relative Risk) Reference Data
#'
#' Loads the rr_list_gbd.xlsx file containing pathogen-drug combinations
#' with relative risk values for burden estimation.
#'
#' @return Data frame with Pathogen, Drug, and RR columns
#' @export
#'
#' @examples
#' \dontrun{
#' rr_data <- daly_load_rr_reference()
#' head(rr_data)
#' }
daly_load_rr_reference <- function() {
  # Find Excel file
  xlsx_path <- find_extdata_file("rr_list_gbd.xlsx")

  if (xlsx_path == "" || !file.exists(xlsx_path)) {
    stop("rr_list_gbd.xlsx not found in inst/extdata/")
  }

  # Load Excel file
  rr_data <- readxl::read_excel(xlsx_path, sheet = "Sheet2")

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
                                 class_col    = "antibiotic_class") {
  # -- Pathogen mapping -------------------------------------------------------
  if (!is.null(organism_col) && organism_col %in% names(data)) {
    taxonomy <- get_organism_taxonomy()
    if (nrow(taxonomy) > 0) {
      rr_map <- unique(data.frame(
        organism_name = taxonomy$organism_name,
        rr_pathogen   = taxonomy$organism_name,
        stringsAsFactors = FALSE
      ))
      data <- dplyr::left_join(
        data, rr_map,
        by = stats::setNames("organism_name", organism_col)
      )
      n_mapped <- sum(!is.na(data$rr_pathogen))
      message(sprintf(
        "Mapped to RR pathogen: %d/%d (%.1f%%)",
        n_mapped, nrow(data), 100 * n_mapped / nrow(data)
      ))
      unmapped <- unique(data[[organism_col]][is.na(data$rr_pathogen)])
      if (length(unmapped) > 0)
        message("Unmapped organisms: ", paste(head(unmapped, 10), collapse = ", "))
    } else {
      warning("Organism taxonomy is empty. Skipping pathogen mapping.")
    }
  }

  # -- Drug class mapping -----------------------------------------------------
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
#' For each resistance profile delta (from compute_resistance_profiles()),
#' determines the profile-level RR_kd_LOS using the GBD max rule:
#'   RR_kd_LOS = max over c in C_R(d) of RR_kc_LOS   [if C_R(d) non-empty]
#'             = 1                                      [if d = all-susceptible]
#' where C_R(d) = \{c : d_c = 1\}.
#' The CI reported for each profile is that of its dominant (max-RR) class.
#'
#' @param profiles_output Named list from compute_resistance_profiles().
#' @param rr_table Data frame from daly_fit_los_rr() or daly_fit_los_rr_nima().
#'   Must have columns pathogen_col, class_col, rr_col, and optionally
#'   CI_lower / CI_upper.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param rr_col Character. Default \code{"RR_LOS"}.
#' @param fallback_rr Numeric. RR for resistant classes with no match.
#'   Default \code{1} (no attributable effect).
#'
#' @return Named list (one entry per pathogen): original profiles data frame
#'   augmented with RR_LOS_profile, dominant_class, and (if available)
#'   CI_lower_profile / CI_upper_profile.
#' @export
daly_assign_rr_to_profiles <- function(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1
) {
  if (!is.list(profiles_output) ||
    !all(sapply(
      profiles_output,
      function(x) all(c("profiles", "classes") %in% names(x))
    ))) {
    stop("profiles_output must be the list returned by compute_resistance_profiles().")
  }
  missing_rr <- setdiff(c(pathogen_col, class_col, rr_col), names(rr_table))
  if (length(missing_rr) > 0L) {
    stop(sprintf(
      "Column(s) not found in rr_table: %s",
      paste(missing_rr, collapse = ", ")
    ))
  }

  has_ci <- all(c("CI_lower", "CI_upper") %in% names(rr_table))
  out <- list()

  for (path in names(profiles_output)) {
    prof_list <- profiles_output[[path]]
    profiles <- prof_list$profiles
    classes <- prof_list$classes # alphabetical; matches profile cols

    rr_k <- rr_table[rr_table[[pathogen_col]] == path, ]
    rr_lookup <- setNames(rr_k[[rr_col]], rr_k[[class_col]])
    ci_lo_lookup <- if (has_ci) setNames(rr_k$CI_lower, rr_k[[class_col]]) else NULL
    ci_hi_lookup <- if (has_ci) setNames(rr_k$CI_upper, rr_k[[class_col]]) else NULL

    n_prof <- nrow(profiles)
    rr_profile <- numeric(n_prof)
    dom_class <- character(n_prof)
    ci_lo_prof <- if (has_ci) numeric(n_prof) else NULL
    ci_hi_prof <- if (has_ci) numeric(n_prof) else NULL

    for (i in seq_len(n_prof)) {
      resist_cls <- classes[as.integer(profiles[i, classes]) == 1L]

      if (length(resist_cls) == 0L) {
        rr_profile[i] <- 1.0
        dom_class[i] <- "all_susceptible"
        if (has_ci) {
          ci_lo_prof[i] <- 1.0
          ci_hi_prof[i] <- 1.0
        }
        next
      }

      rrs_c <- sapply(resist_cls, function(cc) {
        ifelse(cc %in% names(rr_lookup), rr_lookup[[cc]], fallback_rr)
      })
      max_idx <- which.max(rrs_c)
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
    }

    profiles$RR_LOS_profile <- round(rr_profile, 4L)
    profiles$dominant_class <- dom_class
    if (has_ci) {
      profiles$CI_lower_profile <- round(ci_lo_prof, 4L)
      profiles$CI_upper_profile <- round(ci_hi_prof, 4L)
    }

    out[[path]] <- profiles
    message(sprintf(
      "'%s': relative LOS assigned to %d profiles. Max relative LOS = %.4f (dominant class: %s).",
      path, n_prof, max(rr_profile), dom_class[which.max(rr_profile)]
    ))
  }

  return(out)
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
