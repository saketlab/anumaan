# prep_outputs.R
# Layer 7: Analysis-ready output construction
#
# Functions:
#   - prep_filter_minimally_usable
#   - prep_filter_analysis_ready
#   - prep_build_fatal_cohort
#   - prep_build_nonfatal_cohort
#   - prep_attrition_flow
#   - prep_missingness_report
#   - prep_validate_analysis_ready


#' Filter Minimally Usable Records
#'
#' Keeps records that have the bare minimum fields needed for any analysis:
#' a patient identifier, a culture date, an organism name, and at least one
#' non-NA AST result. This is a looser filter than \code{prep_filter_analysis_ready()}.
#'
#' @param data Data frame after preprocessing.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param culture_date_col Character. Culture date column. Default "culture_date".
#' @param organism_col Character. Organism column. Default "organism_name".
#' @param ast_col Character. Harmonized AST column. Default "ast_value_harmonized".
#'
#' @return Filtered data frame with only minimally usable rows.
#' @export
prep_filter_minimally_usable <- function(data,
                                          patient_col      = "patient_id",
                                          culture_date_col = "culture_date",
                                          organism_col     = "organism_name",
                                          ast_col          = "ast_value_harmonized") {
  n_before <- nrow(data)

  keep <- rep(TRUE, n_before)

  if (patient_col %in% names(data))
    keep <- keep & !is.na(data[[patient_col]]) & trimws(data[[patient_col]]) != ""

  if (culture_date_col %in% names(data))
    keep <- keep & !is.na(data[[culture_date_col]])

  if (organism_col %in% names(data))
    keep <- keep & !is.na(data[[organism_col]]) & trimws(data[[organism_col]]) != ""

  if (ast_col %in% names(data))
    keep <- keep & !is.na(data[[ast_col]]) & data[[ast_col]] %in% c("S", "I", "R")

  result   <- data[keep, , drop = FALSE]
  n_after  <- nrow(result)
  n_removed <- n_before - n_after

  message(sprintf("[prep_filter_minimally_usable] %d rows retained; %d removed (%.1f%%).",
                  n_after, n_removed, 100 * n_removed / max(n_before, 1)))

  return(result)
}


#' Filter Analysis-Ready Records
#'
#' Keeps records that meet the full burden-analysis readiness criteria:
#' \itemize{
#'   \item Patient ID present
#'   \item Culture date present
#'   \item Organism name present and standardized
#'   \item Antibiotic name present and standardized
#'   \item Harmonized AST value in S/I/R
#'   \item Not flagged as a contaminant (if column present)
#' }
#'
#' @param data Data frame after preprocessing.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param culture_date_col Character. Culture date column. Default "culture_date".
#' @param organism_col Character. Organism column. Default "organism_name".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_name".
#' @param ast_col Character. Harmonized AST column. Default "ast_value_harmonized".
#' @param exclude_contaminants Logical. Exclude rows flagged as contaminants.
#'   Default TRUE.
#' @param contaminant_col Character. Contaminant flag column. Default "contaminant_flag".
#'
#' @return Filtered data frame.
#' @export
prep_filter_analysis_ready <- function(data,
                                        patient_col          = "patient_id",
                                        culture_date_col     = "culture_date",
                                        organism_col         = "organism_name",
                                        antibiotic_col       = "antibiotic_name",
                                        ast_col              = "ast_value_harmonized",
                                        exclude_contaminants = TRUE,
                                        contaminant_col      = "contaminant_flag") {
  n_before <- nrow(data)
  keep     <- rep(TRUE, n_before)

  check_col <- function(col, not_na = TRUE, valid_vals = NULL) {
    if (!col %in% names(data)) return()
    if (not_na)
      keep <<- keep & !is.na(data[[col]]) & trimws(as.character(data[[col]])) != ""
    if (!is.null(valid_vals))
      keep <<- keep & !is.na(data[[col]]) & data[[col]] %in% valid_vals
  }

  check_col(patient_col)
  check_col(culture_date_col)
  check_col(organism_col)
  check_col(antibiotic_col)
  check_col(ast_col, valid_vals = c("S", "I", "R"))

  if (exclude_contaminants && contaminant_col %in% names(data))
    keep <- keep & (is.na(data[[contaminant_col]]) | !data[[contaminant_col]])

  result    <- data[keep, , drop = FALSE]
  n_after   <- nrow(result)
  n_removed <- n_before - n_after

  message(sprintf("[prep_filter_analysis_ready] %d rows retained; %d removed (%.1f%%).",
                  n_after, n_removed, 100 * n_removed / max(n_before, 1)))

  return(result)
}


#' Build Fatal Cohort
#'
#' Extracts records for patients who died, for use in YLL calculation.
#'
#' @param data Analysis-ready data frame.
#' @param outcome_col Character. Final outcome column. Default "final_outcome".
#' @param died_value Character. Value indicating death. Default "Died".
#'
#' @return Data frame with only fatal episodes.
#' @export
prep_build_fatal_cohort <- function(data,
                                     outcome_col = "final_outcome",
                                     died_value  = "Died") {
  if (!outcome_col %in% names(data)) {
    warning(sprintf("[prep_build_fatal_cohort] Column '%s' not found.", outcome_col))
    return(data[integer(0), , drop = FALSE])
  }

  result <- data[!is.na(data[[outcome_col]]) & data[[outcome_col]] == died_value, , drop = FALSE]

  message(sprintf("[prep_build_fatal_cohort] %d fatal records (%.1f%% of input).",
                  nrow(result), 100 * nrow(result) / max(nrow(data), 1)))
  return(result)
}


#' Build Non-Fatal Cohort
#'
#' Extracts records for patients who survived, for use in YLD calculation.
#'
#' @param data Analysis-ready data frame.
#' @param outcome_col Character. Final outcome column. Default "final_outcome".
#' @param survived_value Character. Value indicating survival. Default "Survived".
#'
#' @return Data frame with only non-fatal episodes.
#' @export
prep_build_nonfatal_cohort <- function(data,
                                        outcome_col     = "final_outcome",
                                        survived_value  = "Survived") {
  if (!outcome_col %in% names(data)) {
    warning(sprintf("[prep_build_nonfatal_cohort] Column '%s' not found.", outcome_col))
    return(data[integer(0), , drop = FALSE])
  }

  result <- data[!is.na(data[[outcome_col]]) & data[[outcome_col]] == survived_value, , drop = FALSE]

  message(sprintf("[prep_build_nonfatal_cohort] %d non-fatal records (%.1f%% of input).",
                  nrow(result), 100 * nrow(result) / max(nrow(data), 1)))
  return(result)
}


#' Track Attrition Through the Pipeline
#'
#' Records a stage entry in an attrition flow table. Call once per pipeline
#' stage to build up a cumulative record of how many rows and patients are
#' retained after each step.
#'
#' Usage pattern:
#' \preformatted{
#'   flow <- NULL
#'   flow <- prep_attrition_flow(flow, data_raw,    "raw_intake",    "All raw records")
#'   flow <- prep_attrition_flow(flow, data_dated,  "dates_parsed",  "After date coercion")
#'   flow <- prep_attrition_flow(flow, data_ready,  "analysis_ready","After all filters")
#'   print(flow)
#' }
#'
#' @param flow Data frame of previous attrition stages (or NULL to start).
#' @param data Current data frame at this stage.
#' @param stage_name Character. Short label for this stage.
#' @param reason Character. Description of what changed at this stage.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param event_col Character. Event ID column (optional). Default "event_id".
#'
#' @return Updated attrition flow data frame.
#' @export
prep_attrition_flow <- function(flow,
                                 data,
                                 stage_name,
                                 reason     = "",
                                 patient_col = "patient_id",
                                 event_col   = "event_id") {
  n_rows     <- nrow(data)
  n_patients <- if (patient_col %in% names(data))
    dplyr::n_distinct(data[[patient_col]], na.rm = TRUE) else NA_integer_
  n_events   <- if (event_col %in% names(data))
    dplyr::n_distinct(data[[event_col]], na.rm = TRUE) else NA_integer_

  n_removed <- if (!is.null(flow) && nrow(flow) > 0)
    flow$n_rows[nrow(flow)] - n_rows else 0L

  new_row <- data.frame(
    stage      = stage_name,
    n_rows     = n_rows,
    n_patients = n_patients,
    n_events   = n_events,
    n_removed  = n_removed,
    reason     = reason,
    stringsAsFactors = FALSE
  )

  if (is.null(flow)) new_row else rbind(flow, new_row)
}


#' Missingness Report
#'
#' Generates a per-column summary of missing values. Flags columns where
#' missingness exceeds a threshold.
#'
#' @param data Data frame.
#' @param threshold Numeric. Percentage threshold (0-100) above which a column
#'   is flagged as high-missingness. Default 20.
#' @param cols Character vector. Subset of columns to report on. NULL = all.
#'
#' @return Data frame with columns: col_name, n_total, n_missing, pct_missing,
#'   is_high_missing.
#' @export
prep_missingness_report <- function(data, threshold = 20, cols = NULL) {
  target_cols <- if (!is.null(cols)) intersect(cols, names(data)) else names(data)

  report <- do.call(rbind, lapply(target_cols, function(col) {
    n_total   <- nrow(data)
    n_missing <- sum(is.na(data[[col]]) | trimws(as.character(data[[col]])) %in%
                       c("", "NA", "NULL", "N/A", "None"))
    pct       <- 100 * n_missing / max(n_total, 1)
    data.frame(
      col_name        = col,
      n_total         = n_total,
      n_missing       = n_missing,
      pct_missing     = round(pct, 1),
      is_high_missing = pct > threshold,
      stringsAsFactors = FALSE
    )
  }))

  n_high <- sum(report$is_high_missing)
  if (n_high > 0)
    message(sprintf("[prep_missingness_report] %d column(s) exceed %.0f%% missingness threshold.",
                    n_high, threshold))

  report[order(report$pct_missing, decreasing = TRUE), ]
}


#' Validate Analysis-Ready Dataset
#'
#' Final gate check before burden estimation. Confirms that the minimum standard
#' AMR schema fields are present and have acceptable completeness.
#'
#' @param data Analysis-ready data frame.
#' @param min_rows Integer. Minimum acceptable row count. Default 10.
#' @param max_missing_pct Numeric. Maximum percent missing allowed for core fields.
#'   Default 30.
#' @param stop_on_failure Logical. Stop if validation fails. Default FALSE.
#'
#' @return Invisibly returns a list: passes (logical), issues (character vector).
#' @export
prep_validate_analysis_ready <- function(data,
                                          min_rows         = 10L,
                                          max_missing_pct  = 30,
                                          stop_on_failure  = FALSE) {
  issues <- character(0)

  # Row count
  if (nrow(data) < min_rows)
    issues <- c(issues, sprintf("Only %d rows present (minimum: %d).", nrow(data), min_rows))

  # Core required fields for burden analysis
  core_fields <- c("patient_id", "culture_date", "organism_name",
                   "antibiotic_name", "ast_value_harmonized")

  missing_fields <- setdiff(core_fields, names(data))
  if (length(missing_fields) > 0)
    issues <- c(issues, sprintf("Missing core fields: %s.",
                                paste(missing_fields, collapse = ", ")))

  # Missingness in present core fields
  present_core <- intersect(core_fields, names(data))
  for (col in present_core) {
    pct <- 100 * sum(is.na(data[[col]])) / max(nrow(data), 1)
    if (pct > max_missing_pct)
      issues <- c(issues, sprintf("'%s' has %.1f%% missing (threshold: %.0f%%).",
                                  col, pct, max_missing_pct))
  }

  passes <- length(issues) == 0L

  if (passes) {
    message("[prep_validate_analysis_ready] Dataset passes all validation checks.")
  } else {
    msg <- paste("[prep_validate_analysis_ready] Validation issues found:\n",
                 paste(" -", issues, collapse = "\n"))
    if (stop_on_failure) stop(msg) else warning(msg)
  }

  invisible(list(passes = passes, issues = issues))
}
