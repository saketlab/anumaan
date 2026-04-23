# prep_intake.R
# Layer 1: File/table intake checks
# Purpose: confirm the file loaded, has usable structure, and columns exist.
# No data transformation happens here. Only checks.
#
# Functions moved here from:
#   - prep_stewardship_join.R : prep_check_columns, prep_check_keys, prep_validate_table
#   - prep_clean_and_standardize.R : validate_required_fields, validate_data_quality
#
# New functions: prep_log_source, prep_inventory_columns, prep_detect_schema_drift


#' Check required columns exist and report types
#'
#' Runs before any table is merged. Reports missing columns, type mismatches,
#' and blank/duplicate column names.  Stops on hard failures; warns on soft ones.
#'
#' @param data Data frame to inspect.
#' @param required Character vector of column names that must be present.
#' @param expected_types Named character vector mapping column names to expected
#'   R classes (e.g. \code{c(PatientInformation_id = "character")}). Only
#'   checked when the column exists.
#' @param table_label Character. Label used in messages (e.g. "patient_KGMU").
#' @param stop_on_missing Logical. Stop if required columns are absent (default TRUE).
#'
#' @return Invisibly returns a tibble summarising every checked column.
#' @export
prep_check_columns <- function(data,
                               required = character(),
                               expected_types = character(),
                               table_label = "table",
                               stop_on_missing = TRUE) {
  if (!is.data.frame(data)) {
    stop(sprintf("[%s] `data` must be a data frame.", table_label))
  }

  blank_names <- sum(is.na(names(data)) | trimws(names(data)) == "")
  dup_names   <- sum(duplicated(names(data)))
  if (blank_names > 0L)
    warning(sprintf("[%s] %d blank column name(s) found.", table_label, blank_names))
  if (dup_names > 0L)
    warning(sprintf("[%s] %d duplicated column name(s): %s",
                    table_label, dup_names,
                    paste(names(data)[duplicated(names(data))], collapse = ", ")))

  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    msg <- sprintf("[%s] Missing required column(s): %s",
                   table_label, paste(missing, collapse = ", "))
    if (stop_on_missing) stop(msg) else warning(msg)
  }

  report_rows <- list()
  for (col in names(data)) {
    actual_class   <- paste(class(data[[col]]), collapse = "/")
    expected_class <- if (col %in% names(expected_types)) expected_types[[col]] else NA_character_
    type_ok <- if (!is.na(expected_class)) inherits(data[[col]], expected_class) else TRUE
    n_na    <- sum(is.na(data[[col]]))
    n_total <- nrow(data)
    pct_na  <- if (n_total > 0) round(100 * n_na / n_total, 1) else NA_real_

    report_rows[[length(report_rows) + 1L]] <- data.frame(
      table         = table_label,
      column        = col,
      required      = col %in% required,
      present       = TRUE,
      actual_type   = actual_class,
      expected_type = if (!is.na(expected_class)) expected_class else "",
      type_ok       = type_ok,
      n_total       = n_total,
      n_na          = n_na,
      pct_na        = pct_na,
      stringsAsFactors = FALSE
    )

    if (!type_ok)
      warning(sprintf("[%s] Column '%s': expected class '%s', found '%s'.",
                      table_label, col, expected_class, actual_class))
  }

  for (col in missing) {
    report_rows[[length(report_rows) + 1L]] <- data.frame(
      table         = table_label,
      column        = col,
      required      = TRUE,
      present       = FALSE,
      actual_type   = NA_character_,
      expected_type = if (col %in% names(expected_types)) expected_types[[col]] else "",
      type_ok       = FALSE,
      n_total       = NA_integer_,
      n_na          = NA_integer_,
      pct_na        = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  report <- dplyr::bind_rows(report_rows)

  message(sprintf(
    "[%s] Column check: %d cols | %d required | %d missing | %d type warnings",
    table_label, ncol(data), length(required), length(missing),
    sum(!report$type_ok, na.rm = TRUE)
  ))

  invisible(report)
}


#' Check join key quality
#'
#' Reports missing, blank, duplicate, and placeholder keys in a column before
#' a merge. Run this on both sides of every planned join.
#'
#' @param data Data frame.
#' @param key_col Character. Name of the key column.
#' @param table_label Character. Label used in messages.
#' @param warn_missing_pct Numeric. Warn when proportion of missing keys exceeds
#'   this threshold (0-100). Default 5.
#'
#' @return Invisibly returns a one-row summary tibble.
#' @export
prep_check_keys <- function(data,
                            key_col,
                            table_label = "table",
                            warn_missing_pct = 5) {
  if (!key_col %in% names(data)) {
    warning(sprintf("[%s] Key column '%s' not found.", table_label, key_col))
    return(invisible(NULL))
  }

  raw         <- as.character(data[[key_col]])
  placeholders <- c("", "NULL", "null", "NA", "N/A", "None", "none", "nan", "NaN")
  is_missing  <- is.na(raw) | trimws(raw) %in% placeholders
  n_total     <- length(raw)
  n_missing   <- sum(is_missing)
  pct_missing <- round(100 * n_missing / max(n_total, 1), 1)
  n_distinct  <- length(unique(raw[!is_missing]))
  n_dup       <- sum(duplicated(raw[!is_missing]))

  if (pct_missing > warn_missing_pct)
    warning(sprintf("[%s] Key '%s': %.1f%% missing/placeholder (%d of %d rows).",
                    table_label, key_col, pct_missing, n_missing, n_total))

  message(sprintf(
    "[%s] Key '%s': %d rows | %d missing (%.1f%%) | %d distinct | %d duplicated",
    table_label, key_col, n_total, n_missing, pct_missing, n_distinct, n_dup
  ))

  invisible(data.frame(
    table        = table_label,
    key_col      = key_col,
    n_total      = n_total,
    n_missing    = n_missing,
    pct_missing  = pct_missing,
    n_distinct   = n_distinct,
    n_duplicated = n_dup,
    stringsAsFactors = FALSE
  ))
}


#' Run all pre-join sanity checks for one table
#'
#' Convenience wrapper that calls \code{prep_check_columns()},
#' \code{prep_check_keys()}, and \code{prep_coerce_dates()} in sequence.
#' Returns the date-coerced data along with a check report.
#'
#' @param data Data frame.
#' @param required_cols Character vector of required column names.
#' @param key_col Character. Primary/foreign key column to quality-check.
#' @param expected_types Named character vector of expected column classes.
#' @param date_cols Character vector of date column names to coerce. NULL
#'   triggers auto-detection.
#' @param table_label Character. Label used in all messages.
#' @param stop_on_missing Logical. Passed to \code{prep_check_columns()}.
#'
#' @return A list with \code{data} (date-coerced) and \code{report} (check summary).
#' @export
prep_validate_table <- function(data,
                                required_cols   = character(),
                                key_col         = NULL,
                                expected_types  = character(),
                                date_cols       = NULL,
                                table_label     = "table",
                                stop_on_missing = TRUE) {
  col_report <- prep_check_columns(
    data            = data,
    required        = required_cols,
    expected_types  = expected_types,
    table_label     = table_label,
    stop_on_missing = stop_on_missing
  )

  key_report <- NULL
  if (!is.null(key_col)) {
    key_report <- prep_check_keys(
      data        = data,
      key_col     = key_col,
      table_label = table_label
    )
  }

  data <- prep_coerce_dates(data, cols = date_cols, table_label = table_label)

  list(data = data, col_report = col_report, key_report = key_report)
}


#' Validate required fields exist and meet completeness threshold
#'
#' @param data Data frame.
#' @param required_cols Character vector of required column names.
#' @param stop_on_failure Logical. Stop if validation fails. Default TRUE.
#' @param allow_na Logical. Skip completeness checks if TRUE. Default FALSE.
#' @param min_completeness Numeric (0-1). Minimum completeness required. Default 0.8.
#'
#' @return List with validation result details.
#' @export
validate_required_fields <- function(data,
                                     required_cols,
                                     stop_on_failure = TRUE,
                                     allow_na = FALSE,
                                     min_completeness = 0.8) {
  validation_msgs <- character()
  missing_cols    <- character()
  incomplete_cols <- data.frame()

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    msg <- sprintf("Missing required columns: %s",
                   paste(missing_cols, collapse = ", "))
    validation_msgs <- c(validation_msgs, msg)
  }

  existing_cols <- intersect(required_cols, names(data))

  if (length(existing_cols) > 0 && !allow_na) {
    completeness <- sapply(existing_cols, function(col) {
      sum(!is.na(data[[col]])) / nrow(data)
    })

    incomplete <- completeness < min_completeness

    if (any(incomplete)) {
      incomplete_cols <- data.frame(
        column       = names(completeness)[incomplete],
        completeness = completeness[incomplete],
        n_missing    = sapply(names(completeness)[incomplete], function(col) {
          sum(is.na(data[[col]]))
        }),
        stringsAsFactors = FALSE
      )

      msg <- sprintf("Columns below %.0f%% completeness: %s",
                     min_completeness * 100,
                     paste(incomplete_cols$column, collapse = ", "))
      validation_msgs <- c(validation_msgs, msg)
    }
  }

  is_valid <- length(missing_cols) == 0 && nrow(incomplete_cols) == 0

  result <- list(
    valid          = is_valid,
    missing_cols   = missing_cols,
    incomplete_cols = incomplete_cols,
    messages       = validation_msgs,
    n_rows         = nrow(data),
    n_cols_checked = length(required_cols)
  )

  if (is_valid) {
    message(sprintf("[v] Validation passed: All %d required columns present and complete",
                    length(required_cols)))
  } else {
    message("[x] Validation failed:")
    for (msg in validation_msgs) message(sprintf("  - %s", msg))
    if (nrow(incomplete_cols) > 0) {
      message("\nCompleteness details:")
      print(incomplete_cols)
    }
  }

  if (!is_valid && stop_on_failure)
    stop("Data validation failed. See messages above.")

  return(result)
}


#' Validate Data Quality
#'
#' Runs quality checks: minimum row count, required column presence,
#' and maximum missing-value percentage per column.
#'
#' @param data Data frame to validate.
#' @param min_rows Integer. Minimum acceptable number of rows. Default 10.
#' @param max_missing_pct Numeric. Maximum acceptable percent missing per column (0-100). Default 50.
#' @param required_cols Character vector. Columns that must be present.
#' @param stop_on_failure Logical. Stop with error on failure. Default FALSE.
#'
#' @return List with quality assessment.
#' @export
validate_data_quality <- function(data,
                                  min_rows        = 10,
                                  max_missing_pct = 50,
                                  required_cols   = c("patient_id", "organism_normalized"),
                                  stop_on_failure = FALSE) {
  quality_issues <- character()

  n_rows <- nrow(data)
  if (n_rows < min_rows)
    quality_issues <- c(quality_issues,
                        sprintf("Dataset too small: %d rows (minimum: %d)", n_rows, min_rows))

  missing_req <- setdiff(required_cols, names(data))
  if (length(missing_req) > 0)
    quality_issues <- c(quality_issues,
                        sprintf("Missing required columns: %s", paste(missing_req, collapse = ", ")))

  completeness <- sapply(names(data), function(col) {
    100 * sum(!is.na(data[[col]])) / n_rows
  })

  col_completeness <- data.frame(
    column          = names(completeness),
    completeness_pct = as.numeric(completeness),
    n_missing       = sapply(names(data), function(col) sum(is.na(data[[col]]))),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(completeness_pct)

  poor_cols <- col_completeness %>%
    dplyr::filter(completeness_pct < (100 - max_missing_pct))

  if (nrow(poor_cols) > 0)
    quality_issues <- c(quality_issues,
                        sprintf("%d columns exceed %.0f%% missing threshold: %s",
                                nrow(poor_cols), max_missing_pct,
                                paste(poor_cols$column[1:min(5, nrow(poor_cols))], collapse = ", ")))

  overall_completeness <- sum(!is.na(data)) / (nrow(data) * ncol(data))
  if (overall_completeness < 0.5)
    quality_issues <- c(quality_issues,
                        sprintf("Overall completeness too low: %.1f%% (target: >=50%%)",
                                overall_completeness * 100))

  passes_quality <- length(quality_issues) == 0

  result <- list(
    passes_quality       = passes_quality,
    n_rows               = n_rows,
    n_cols               = ncol(data),
    overall_completeness = overall_completeness,
    column_completeness  = col_completeness,
    quality_issues       = quality_issues
  )

  if (passes_quality) {
    message(sprintf("[v] Quality check passed: %d rows x %d cols, %.1f%% complete",
                    n_rows, ncol(data), overall_completeness * 100))
  } else {
    message("[x] Quality issues detected:")
    for (issue in quality_issues) message(sprintf("  - %s", issue))
  }

  if (!passes_quality && stop_on_failure) stop("Data quality validation failed")

  return(result)
}


# ---------------------------------------------------------------------------
# New functions (Layer 1)
# ---------------------------------------------------------------------------

#' Log Data Source Provenance
#'
#' Records metadata about the source of data entering the pipeline.
#' Returns a provenance metadata list attached as an attribute.
#'
#' @param data Data frame.
#' @param study_type Character. One of "stewardship", "surveillance", "aiims_icu_bsi", "generic".
#' @param centre_name Character. Name of the contributing centre.
#' @param file_path Character. Path to the source file (optional).
#' @param sheet_name Character. Sheet name within Excel file (optional).
#'
#' @return Data frame with \code{.provenance} attribute attached.
#' @export
prep_log_source <- function(data,
                            study_type  = "generic",
                            centre_name = NULL,
                            file_path   = NULL,
                            sheet_name  = NULL) {
  provenance <- list(
    centre_name   = centre_name %||% NA_character_,
    study_type    = study_type,
    file_path     = file_path %||% NA_character_,
    sheet_name    = sheet_name %||% NA_character_,
    extracted_at  = Sys.time(),
    n_rows        = nrow(data),
    n_cols        = ncol(data)
  )

  message(sprintf(
    "[prep_log_source] Centre: %s | Study: %s | Rows: %d | Cols: %d",
    provenance$centre_name, study_type, nrow(data), ncol(data)
  ))

  attr(data, ".provenance") <- provenance
  data
}


#' Inventory Columns
#'
#' Returns a tibble describing every column: class, distinct values, and
#' missingness. Used to detect completely empty columns and schema drift.
#'
#' @param data Data frame.
#'
#' @return Tibble: col_name | class | n_distinct | pct_missing | has_any_value
#' @export
prep_inventory_columns <- function(data) {
  n <- nrow(data)
  tibble::tibble(
    col_name      = names(data),
    class         = sapply(data, function(x) paste(class(x), collapse = "/")),
    n_distinct    = sapply(data, function(x) dplyr::n_distinct(x, na.rm = TRUE)),
    n_missing     = sapply(data, function(x) sum(is.na(x))),
    pct_missing   = round(100 * sapply(data, function(x) sum(is.na(x))) / max(n, 1), 1),
    has_any_value = sapply(data, function(x) any(!is.na(x)))
  )
}


#' Detect Schema Drift Across Centres
#'
#' Compares column sets across a list of centre data frames. Flags columns
#' present in some centres but not others.
#'
#' @param data_list Named list of data frames (one per centre).
#' @param reference_centre Character or NULL. Name of the centre to treat as
#'   the reference schema. If NULL, the union of all columns is used.
#'
#' @return Tibble: column | present_in (comma-separated centres) | missing_from | n_centres_present
#' @export
prep_detect_schema_drift <- function(data_list, reference_centre = NULL) {
  if (!is.list(data_list) || length(data_list) == 0L)
    stop("`data_list` must be a non-empty named list of data frames.")

  centre_names <- names(data_list)
  if (is.null(centre_names)) centre_names <- paste0("centre_", seq_along(data_list))

  all_cols <- unique(unlist(lapply(data_list, names)))

  if (!is.null(reference_centre) && reference_centre %in% centre_names) {
    ref_cols <- names(data_list[[reference_centre]])
  } else {
    ref_cols <- all_cols
  }

  presence_matrix <- sapply(centre_names, function(ctr) {
    all_cols %in% names(data_list[[ctr]])
  })

  if (is.vector(presence_matrix))
    presence_matrix <- matrix(presence_matrix, nrow = length(all_cols),
                              dimnames = list(all_cols, centre_names))

  result <- tibble::tibble(
    column           = all_cols,
    n_centres_present = rowSums(presence_matrix),
    present_in       = apply(presence_matrix, 1, function(x)
      paste(centre_names[x], collapse = ", ")),
    missing_from     = apply(presence_matrix, 1, function(x)
      paste(centre_names[!x], collapse = ", ")),
    is_universal     = rowSums(presence_matrix) == length(centre_names),
    is_reference_col = all_cols %in% ref_cols
  ) %>%
    dplyr::arrange(n_centres_present, column)

  n_drift <- sum(!result$is_universal)
  if (n_drift > 0) {
    message(sprintf(
      "[prep_detect_schema_drift] %d column(s) not present in all %d centres.",
      n_drift, length(centre_names)
    ))
  } else {
    message("[prep_detect_schema_drift] Schema is consistent across all centres.")
  }

  result
}
