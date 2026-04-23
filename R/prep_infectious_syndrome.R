# prep_infectious_syndrome.R
# Wide-format syndrome assignment helper.
#
# Functions removed (merged into prep_diagnosis.R):
#   normalize_syndrome_name     -> internal helper in prep_diagnosis.R
#   read_syndrome_hierarchy     -> replaced by load_syndrome_hierarchy()
#   infer_patient_syndrome_long -> merged into prep_assign_patient_syndrome()
#
# Remaining:
#   infer_patient_syndrome_wide -> pivots wide syndrome flags to long, then
#                                  calls prep_assign_patient_syndrome()


#' Assign Syndrome from Wide-Format Syndrome Flags
#'
#' Converts a wide-format data frame (one column per syndrome, values indicating
#' presence) to long format and assigns the highest-priority syndrome per patient
#' using the infectious syndrome hierarchy.
#'
#' This is a convenience wrapper around \code{prep_assign_patient_syndrome()}.
#' All syndrome-selection logic and parameters (hierarchy, burden filter,
#' respiratory collapsing) are handled there.
#'
#' @param data Data frame in wide format where syndrome columns contain presence
#'   indicators (\code{1}/\code{TRUE}/\code{"Yes"} etc.).
#' @param patient_col Character. Patient ID column. Default \code{"patient_id"}.
#' @param syndrome_cols Character vector. Columns to treat as syndrome flags.
#'   If NULL, all columns except \code{patient_col} are used. Default NULL.
#' @param positive_values Character/logical/numeric vector. Values treated as
#'   "syndrome present". Default covers common representations of TRUE/1/Yes.
#' @param collapse_unspecified_respiratory Logical. Passed to
#'   \code{prep_assign_patient_syndrome()}. Default TRUE.
#' @param keep_only_burden_syndromes Logical. Passed to
#'   \code{prep_assign_patient_syndrome()}. Default FALSE.
#'
#' @return Data frame with one row per patient and the selected syndrome.
#' @export
infer_patient_syndrome_wide <- function(data,
                                        patient_col                      = "patient_id",
                                        syndrome_cols                    = NULL,
                                        positive_values                  = c(1, "1", TRUE, "TRUE",
                                                                             "True", "true",
                                                                             "Yes", "YES", "yes"),
                                        collapse_unspecified_respiratory = TRUE,
                                        keep_only_burden_syndromes       = FALSE) {
  if (!patient_col %in% names(data))
    stop(sprintf("[infer_patient_syndrome_wide] Column '%s' not found.", patient_col))

  if (is.null(syndrome_cols))
    syndrome_cols <- setdiff(names(data), patient_col)

  missing_cols <- setdiff(syndrome_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("[infer_patient_syndrome_wide] Missing syndrome columns: %s",
                 paste(missing_cols, collapse = ", ")))

  long_df <- data %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(syndrome_cols),
      names_to  = "syndrome",
      values_to = ".present"
    ) %>%
    dplyr::filter(as.character(.data$.present) %in% as.character(positive_values)) %>%
    dplyr::select(-.present)

  message(sprintf("[infer_patient_syndrome_wide] %d positive syndrome rows after pivot.",
                  nrow(long_df)))

  prep_assign_patient_syndrome(
    data                             = long_df,
    patient_col                      = patient_col,
    syndrome_col                     = "syndrome",
    score_col                        = NULL,
    keep_all_candidates              = FALSE,
    collapse_unspecified_respiratory = collapse_unspecified_respiratory,
    keep_only_burden_syndromes       = keep_only_burden_syndromes
  )
}
