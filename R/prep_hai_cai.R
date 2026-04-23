# prep_hai_cai.R
# Layer 6b: HAI/CAI derivation
#
# Functions:
#   - prep_derive_hai_cai
#   - prep_flag_hai_inferred
#   - prep_reconcile_hai_observed_inferred
#   - prep_derive_icu_flag
#
# Note: prep_check_hai_inputs removed — covered by prep_validate_date_logic().


#' Derive HAI/CAI Infection Type
#'
#' Infers Community-Acquired (CAI) vs Hospital-Acquired (HAI) infection
#' using the admission-to-culture date gap (default cutoff: 2 days / 48 hours).
#'
#' @param data Data frame.
#' @param infection_type_col Character. Infection type column.
#'   Default "infection_type".
#' @param admission_col Character. Admission date column. Default "date_of_admission".
#' @param culture_col Character. Culture date column. Default "date_of_culture".
#' @param hai_cutoff Numeric. Days after admission to classify as HAI. Default 2.
#' @param overwrite Logical. Recalculate even when already present. Default FALSE.
#'
#' @return Data frame with \code{infection_type} enriched.
#' @export
prep_derive_hai_cai <- function(data,
                                 infection_type_col = "infection_type",
                                 admission_col      = "date_of_admission",
                                 culture_col        = "date_of_culture",
                                 hai_cutoff         = 2,
                                 overwrite          = FALSE) {
  has_infection_type <- infection_type_col %in% names(data)
  has_admission      <- admission_col      %in% names(data)
  has_culture        <- culture_col        %in% names(data)

  if (!has_infection_type) data[[infection_type_col]] <- NA_character_

  if (!has_admission || !has_culture) {
    message(sprintf("[!] Cannot infer infection type: missing '%s' or '%s'",
                    admission_col, culture_col))
    return(data)
  }

  n_before_missing <- sum(is.na(data[[infection_type_col]]))
  message(sprintf("Inferring infection type using %d-day HAI cutoff...", hai_cutoff))

  data <- data %>%
    dplyr::mutate(
      days_to_culture = dplyr::case_when(
        !is.na(!!rlang::sym(admission_col)) & !is.na(!!rlang::sym(culture_col)) ~
          as.numeric(difftime(!!rlang::sym(culture_col), !!rlang::sym(admission_col), units = "days")),
        TRUE ~ NA_real_
      ),
      inferred_type = dplyr::case_when(
        !is.na(days_to_culture) & days_to_culture >= hai_cutoff ~ "HAI",
        !is.na(days_to_culture) & days_to_culture <  hai_cutoff ~ "CAI",
        TRUE                                                      ~ NA_character_
      ),
      !!infection_type_col := dplyr::case_when(
        overwrite & !is.na(inferred_type)                                           ~ inferred_type,
        is.na(!!rlang::sym(infection_type_col)) & !is.na(inferred_type)             ~ inferred_type,
        TRUE                                                                          ~ !!rlang::sym(infection_type_col)
      ),
      infection_type_method = dplyr::case_when(
        !is.na(inferred_type)                          ~ sprintf("inferred_%dday_cutoff", hai_cutoff),
        !is.na(!!rlang::sym(infection_type_col))       ~ "provided",
        TRUE                                            ~ NA_character_
      ),
      infection_type_confidence = dplyr::case_when(
        infection_type_method == "provided" ~ "high",
        !is.na(inferred_type)               ~ "medium",
        TRUE                                 ~ NA_character_
      )
    ) %>%
    dplyr::select(-inferred_type, -days_to_culture)

  n_enriched <- n_before_missing - sum(is.na(data[[infection_type_col]]))
  if (n_enriched > 0)
    message(sprintf("Enriched infection_type: %d rows filled", n_enriched))

  message("\nInfection type distribution:")
  print(dplyr::arrange(
    dplyr::count(data, !!rlang::sym(infection_type_col), infection_type_method),
    dplyr::desc(n)
  ))

  return(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 6b)
# ---------------------------------------------------------------------------

# prep_check_hai_inputs removed: its admission <= culture check is now covered
# by prep_validate_date_logic() in prep_dates.R.

#' Flag HAI Inferred vs Observed
#'
#' Adds an \code{infection_type_src} column indicating whether the infection type
#' was present in the raw data ("observed"), derived from dates ("inferred"), or
#' could not be determined ("unknown").
#'
#' @param data Data frame.
#' @param infection_type_col Character. Raw/observed infection type column.
#'   Default "infection_type".
#' @param infection_type_method_col Character. Method column added by
#'   \code{prep_derive_hai_cai()}. Default "infection_type_method".
#'
#' @return Data frame with \code{infection_type_src} column added.
#' @export
prep_flag_hai_inferred <- function(data,
                                    infection_type_col        = "infection_type",
                                    infection_type_method_col = "infection_type_method") {
  if (!infection_type_col %in% names(data)) {
    warning(sprintf("[prep_flag_hai_inferred] Column '%s' not found.", infection_type_col))
    data$infection_type_src <- NA_character_
    return(data)
  }

  method_present <- infection_type_method_col %in% names(data)

  data$infection_type_src <- dplyr::case_when(
    method_present & !is.na(data[[infection_type_method_col]]) &
      data[[infection_type_method_col]] == "provided"    ~ "observed",
    method_present & !is.na(data[[infection_type_method_col]]) &
      grepl("inferred", data[[infection_type_method_col]])  ~ "inferred",
    !is.na(data[[infection_type_col]])                   ~ "observed",
    TRUE                                                  ~ "unknown"
  )

  dist <- table(data$infection_type_src, useNA = "ifany")
  message("[prep_flag_hai_inferred] infection_type_src distribution:")
  print(dist)

  return(data)
}


#' Reconcile Observed and Inferred HAI/CAI Classification
#'
#' Where both an observed (raw) infection type and an inferred (date-derived)
#' type exist, flags rows where they disagree. Retains observed classification
#' by default and marks discordant records.
#'
#' @param data Data frame.
#' @param observed_col Character. Observed infection type column. Default "infection_type_observed".
#' @param inferred_col Character. Inferred infection type column. Default "infection_type_inferred".
#' @param output_col Character. Reconciled output column. Default "infection_type".
#'
#' @return Data frame with reconciled \code{infection_type} and
#'   \code{hai_discordant} logical flag.
#' @export
prep_reconcile_hai_observed_inferred <- function(data,
                                                  observed_col = "infection_type_observed",
                                                  inferred_col = "infection_type_inferred",
                                                  output_col   = "infection_type") {
  if (!any(c(observed_col, inferred_col) %in% names(data))) {
    message("[prep_reconcile_hai_observed_inferred] Neither observed nor inferred column found. Skipping.")
    return(data)
  }

  has_obs <- observed_col %in% names(data)
  has_inf <- inferred_col %in% names(data)

  if (!has_obs) {
    data[[output_col]] <- data[[inferred_col]]
    data$hai_discordant <- FALSE
    return(data)
  }
  if (!has_inf) {
    data[[output_col]] <- data[[observed_col]]
    data$hai_discordant <- FALSE
    return(data)
  }

  # Observed takes precedence; flag disagreements
  data[[output_col]] <- dplyr::coalesce(data[[observed_col]], data[[inferred_col]])

  data$hai_discordant <- !is.na(data[[observed_col]]) &
    !is.na(data[[inferred_col]]) &
    data[[observed_col]] != data[[inferred_col]]

  n_discordant <- sum(data$hai_discordant, na.rm = TRUE)
  if (n_discordant > 0)
    warning(sprintf("[prep_reconcile_hai_observed_inferred] %d row(s) have discordant observed vs inferred HAI/CAI.",
                    n_discordant))
  else
    message("[prep_reconcile_hai_observed_inferred] No discordant HAI/CAI classifications found.")

  return(data)
}


#' Derive ICU Flag
#'
#' Sets \code{icu_flag = TRUE} where ward/department information indicates
#' the patient was in an ICU at the time of culture.
#'
#' @param data Data frame.
#' @param ward_col Character. Ward/unit column. Default "ward_icu".
#' @param dept_col Character. Department column. Default "hospital_department".
#'
#' @return Data frame with \code{icu_flag} logical column added.
#' @export
prep_derive_icu_flag <- function(data,
                                  ward_col = "ward_icu",
                                  dept_col = "hospital_department") {
  icu_pattern <- "icu|intensive care|critical care|intensive therapy|itcu|nicu|picu|sicu|cicu"

  icu_from_ward <- if (ward_col %in% names(data))
    !is.na(data[[ward_col]]) &
      (grepl(icu_pattern, tolower(data[[ward_col]]), ignore.case = TRUE) |
         toupper(trimws(data[[ward_col]])) == "ICU")
  else
    rep(FALSE, nrow(data))

  icu_from_dept <- if (dept_col %in% names(data))
    !is.na(data[[dept_col]]) &
      grepl(icu_pattern, tolower(data[[dept_col]]), ignore.case = TRUE)
  else
    rep(FALSE, nrow(data))

  data$icu_flag <- icu_from_ward | icu_from_dept

  n_icu <- sum(data$icu_flag, na.rm = TRUE)
  message(sprintf("[prep_derive_icu_flag] %d row(s) flagged as ICU.", n_icu))

  return(data)
}
