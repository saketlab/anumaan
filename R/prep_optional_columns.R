# prep_types.R
# Layer 3b: Type coercion - sex, outcome, infection type
#
# Functions:
#   - prep_standardize_sex
#   - prep_standardize_final_outcome  (canonical outcome standardizer)
#   - prep_standardize_infection_type
#   - prep_standardize_outcome        (deprecated wrapper -> prep_standardize_final_outcome)


#' Standardize Sex Values
#'
#' Maps various gender/sex representations to standard "M" or "F".
#'
#' @param data Data frame containing gender column.
#' @param col Character. Name of gender column. Default "gender".
#'
#' @return Data frame with standardized gender values.
#' @export
prep_standardize_sex <- function(data, col = "gender") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  data[[col]] <- dplyr::case_when(
    toupper(data[[col]]) %in% c("M", "MALE", "MAN")        ~ "M",
    toupper(data[[col]]) %in% c("F", "FEMALE", "WOMAN")    ~ "F",
    TRUE                                                     ~ NA_character_
  )

  na_count <- sum(is.na(data[[col]]))
  if (na_count > 0)
    warning(sprintf("Gender column: %d values could not be standardized to M/F", na_count))

  message(sprintf("Standardized gender: M=%d, F=%d, NA=%d",
                  sum(data[[col]] == "M", na.rm = TRUE),
                  sum(data[[col]] == "F", na.rm = TRUE),
                  na_count))
  return(data)
}


#' Standardize Outcome Values (deprecated)
#'
#' Deprecated. Use \code{prep_standardize_final_outcome()} instead, which
#' adds \code{outcome_std} ("Survived"/"Died"/NA) without overwriting the
#' original column.
#'
#' @param data Data frame containing outcome column.
#' @param col Character. Name of outcome column. Default "final_outcome".
#'
#' @return Data frame with \code{outcome_std} column added.
#' @export
prep_standardize_outcome <- function(data, col = "final_outcome") {
  .Deprecated("prep_standardize_final_outcome")
  prep_standardize_final_outcome(data, col = col)
}


# ---------------------------------------------------------------------------
# New functions (Layer 3b)
# ---------------------------------------------------------------------------

#' Standardize Final Outcome Column
#'
#' Normalizes outcome variants to a two-level controlled vocabulary:
#'   "Survived", "Died", or NA.
#'
#' \describe{
#'   \item{Survived}{Survived / Alive / Discharge / Discharged / Recovered}
#'   \item{Died}{Died / Death / Expired / Deceased / Dead}
#'   \item{NA}{Unknown / Missing / Absconded / LAMA / DAMA}
#' }
#'
#' @param data Data frame.
#' @param col Character. Outcome column. Default "final_outcome".
#'
#' @return Data frame with \code{outcome_std} column added.
#' @export
prep_standardize_final_outcome <- function(data, col = "final_outcome") {
  if (!col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping outcome standardization.", col))
    data$outcome_std <- NA_character_
    return(data)
  }

  val_up <- toupper(trimws(as.character(data[[col]])))

  data$outcome_std <- dplyr::case_when(
    val_up %in% c("SURVIVED", "ALIVE", "DISCHARGE", "DISCHARGED", "RECOVERED",
                  "DISCHARGED ALIVE", "SURVIVED/DISCHARGED") ~ "Survived",
    val_up %in% c("DIED", "DEATH", "EXPIRED", "DECEASED", "DEAD",
                  "DIED/EXPIRED")                             ~ "Died",
    val_up %in% c("UNKNOWN", "MISSING", "ABSCONDED", "LAMA", "DAMA",
                  "LEFT AGAINST MEDICAL ADVICE", "NA", "")    ~ NA_character_,
    TRUE                                                       ~ NA_character_
  )

  dist <- table(data$outcome_std, useNA = "ifany")
  message("Standardized final outcome distribution:")
  print(dist)

  return(data)
}


#' Standardize Infection Type Column
#'
#' Normalizes HAI/CAI variants to a controlled vocabulary: "HAI", "CAI", or NA.
#'
#' @param data Data frame.
#' @param col Character. Infection type column. Default "infection_type".
#'
#' @return Data frame with standardized \code{infection_type} column.
#' @export
prep_standardize_infection_type <- function(data, col = "infection_type") {
  if (!col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping infection type standardization.", col))
    return(data)
  }

  val_up <- toupper(trimws(as.character(data[[col]])))

  data[[col]] <- dplyr::case_when(
    val_up %in% c("HAI", "HOSPITAL", "HOSPITAL-ACQUIRED", "HOSPITAL ACQUIRED",
                  "NOSOCOMIAL", "HEALTHCARE ASSOCIATED",
                  "HEALTHCARE-ASSOCIATED")                    ~ "HAI",
    val_up %in% c("CAI", "COMMUNITY", "COMMUNITY-ACQUIRED",
                  "COMMUNITY ACQUIRED")                       ~ "CAI",
    TRUE                                                       ~ NA_character_
  )

  dist <- table(data[[col]], useNA = "ifany")
  message("Standardized infection type distribution:")
  print(dist)

  return(data)
}
