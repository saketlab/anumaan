# validate.R
# Data validation and quality control functions

#' Validate Required Fields
#'
#' Checks that all required columns are present and contain valid data.
#' Returns validation results and optionally stops on failure.
#'
#' @param data Data frame to validate
#' @param required_cols Character vector. Required column names.
#' @param stop_on_failure Logical. If TRUE, stops execution on validation failure.
#'   If FALSE, returns validation report. Default TRUE.
#' @param allow_na Logical. If TRUE, allows NA values in required columns.
#'   Default FALSE.
#' @param min_completeness Numeric. Minimum proportion of non-NA values required
#'   (0-1). Default 0.8 (80 percent completeness).
#'
#' @return List with validation results:
#' \itemize{
#'   \item \code{valid}: Logical. Overall validation status
#'   \item \code{missing_cols}: Character vector of missing columns
#'   \item \code{incomplete_cols}: Data frame of columns below min_completeness
#'   \item \code{messages}: Character vector of validation messages
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Stop on failure (default)
#' validate_required_fields(
#'   data,
#'   required_cols = c("patient_id", "organism_normalized", "antibiotic_normalized")
#' )
#'
#' # Get validation report without stopping
#' validation <- validate_required_fields(
#'   data,
#'   required_cols = c("patient_id", "organism_normalized"),
#'   stop_on_failure = FALSE
#' )
#' }
validate_required_fields <- function(data,
                                     required_cols,
                                     stop_on_failure = TRUE,
                                     allow_na = FALSE,
                                     min_completeness = 0.8) {
  validation_msgs <- character()
  missing_cols <- character()
  incomplete_cols <- data.frame()

  # Check 1: Column existence
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    msg <- sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
    validation_msgs <- c(validation_msgs, msg)
  }

  # Check 2: Completeness (for existing columns)
  existing_cols <- intersect(required_cols, names(data))

  if (length(existing_cols) > 0 && !allow_na) {
    completeness <- sapply(existing_cols, function(col) {
      sum(!is.na(data[[col]])) / nrow(data)
    })

    incomplete <- completeness < min_completeness

    if (any(incomplete)) {
      incomplete_cols <- data.frame(
        column = names(completeness)[incomplete],
        completeness = completeness[incomplete],
        n_missing = sapply(names(completeness)[incomplete], function(col) {
          sum(is.na(data[[col]]))
        }),
        stringsAsFactors = FALSE
      )

      msg <- sprintf(
        "Columns below %.0f%% completeness: %s",
        min_completeness * 100,
        paste(incomplete_cols$column, collapse = ", ")
      )
      validation_msgs <- c(validation_msgs, msg)
    }
  }

  # Determine overall validity
  is_valid <- length(missing_cols) == 0 && nrow(incomplete_cols) == 0

  # Create validation result
  result <- list(
    valid = is_valid,
    missing_cols = missing_cols,
    incomplete_cols = incomplete_cols,
    messages = validation_msgs,
    n_rows = nrow(data),
    n_cols_checked = length(required_cols)
  )

  # Print messages
  if (is_valid) {
    message(sprintf(
      "[v] Validation passed: All %d required columns present and complete",
      length(required_cols)
    ))
  } else {
    message("[x] Validation failed:")
    for (msg in validation_msgs) {
      message(sprintf("  - %s", msg))
    }

    if (nrow(incomplete_cols) > 0) {
      message("\nCompleteness details:")
      print(incomplete_cols)
    }
  }

  # Stop or return
  if (!is_valid && stop_on_failure) {
    stop("Data validation failed. See messages above.")
  }

  return(result)
}


#' Remove Duplicate Rows
#'
#' Identifies and removes exact duplicate rows from the dataset.
#' Optionally keeps first or last occurrence.
#'
#' @param data Data frame
#' @param keep Character. Which duplicate to keep: "first" (default), "last", or "none".
#' @param subset Character vector. Column names to check for duplicates.
#'   If NULL, checks all columns. Default NULL.
#' @param report Logical. If TRUE, prints detailed duplicate report. Default TRUE.
#'
#' @return Data frame with duplicates removed
#' @export
#'
#' @examples
#' \dontrun{
#' # Remove exact duplicates (all columns)
#' clean_data <- remove_duplicate_rows(data)
#'
#' # Remove duplicates based on specific columns
#' clean_data <- remove_duplicate_rows(
#'   data,
#'   subset = c("patient_id", "date_of_culture", "organism_normalized")
#' )
#' }
remove_duplicate_rows <- function(data,
                                  keep = "first",
                                  subset = NULL,
                                  report = TRUE) {
  # Validate keep parameter
  valid_keep <- c("first", "last", "none")
  if (!keep %in% valid_keep) {
    stop(sprintf("keep must be one of: %s", paste(valid_keep, collapse = ", ")))
  }

  # Validate subset columns
  if (!is.null(subset)) {
    missing_cols <- setdiff(subset, names(data))
    if (length(missing_cols) > 0) {
      stop(sprintf("Subset columns not found: %s", paste(missing_cols, collapse = ", ")))
    }
  }

  n_before <- nrow(data)

  # Identify duplicates
  if (is.null(subset)) {
    is_dup <- duplicated(data, fromLast = (keep == "last"))
  } else {
    is_dup <- duplicated(data[, subset], fromLast = (keep == "last"))
  }

  n_duplicates <- sum(is_dup)

  # Report duplicates before removal
  if (report && n_duplicates > 0) {
    message(sprintf(
      "Found %d duplicate rows (%.1f%%)",
      n_duplicates,
      100 * n_duplicates / n_before
    ))

    # Show example duplicates
    dup_indices <- which(is_dup | duplicated(data[, subset], fromLast = !is_dup))
    example_dups <- data[dup_indices[1:min(6, length(dup_indices))], ]

    if (!is.null(subset)) {
      example_dups <- example_dups %>%
        dplyr::select(dplyr::all_of(subset), dplyr::everything())
    }

    message("\nExample duplicates:")
    print(utils::head(example_dups))
  }

  # Remove duplicates
  if (keep == "none") {
    # Remove all duplicates (keep no occurrences)
    if (is.null(subset)) {
      data <- data[!duplicated(data) & !duplicated(data, fromLast = TRUE), ]
    } else {
      data <- data[!duplicated(data[, subset]) & !duplicated(data[, subset], fromLast = TRUE), ]
    }
  } else {
    # Keep first or last
    data <- data[!is_dup, ]
  }

  n_after <- nrow(data)
  n_removed <- n_before - n_after

  message(sprintf(
    "Removed %d duplicate rows: %d -> %d (%.1f%% retained)",
    n_removed, n_before, n_after, 100 * n_after / n_before
  ))

  return(data)
}


#' Validate Data Quality
#'
#' Assesses overall data quality metrics and flags datasets that don't
#' meet minimum quality thresholds.
#'
#' @param data Data frame
#' @param min_rows Numeric. Minimum number of rows required. Default 10.
#' @param max_missing_pct Numeric. Maximum percent of missing values allowed
#'   per column (0-100). Default 50.
#' @param required_cols Character vector. Columns that must be present.
#'   Default c("patient_id", "organism_normalized").
#' @param stop_on_failure Logical. If TRUE, stops on quality failure. Default FALSE.
#'
#' @return List with quality assessment:
#'   - passes_quality: Logical
#'   - n_rows: Number of rows
#'   - n_cols: Number of columns
#'   - overall_completeness: Overall proportion non-missing (0-1)
#'   - column_completeness: Data frame of per-column completeness
#'   - quality_issues: Character vector of identified issues
#'
#' @export
#'
#' @examples
#' \dontrun{
#' quality <- validate_data_quality(data, min_rows = 100, max_missing_pct = 30)
#'
#' if (!quality$passes_quality) {
#'   print(quality$quality_issues)
#' }
#' }
validate_data_quality <- function(data,
                                  min_rows = 10,
                                  max_missing_pct = 50,
                                  required_cols = c("patient_id", "organism_normalized"),
                                  stop_on_failure = FALSE) {
  quality_issues <- character()

  # Check 1: Minimum rows
  n_rows <- nrow(data)
  if (n_rows < min_rows) {
    quality_issues <- c(
      quality_issues,
      sprintf("Dataset too small: %d rows (minimum: %d)", n_rows, min_rows)
    )
  }

  # Check 2: Required columns
  missing_req <- setdiff(required_cols, names(data))
  if (length(missing_req) > 0) {
    quality_issues <- c(
      quality_issues,
      sprintf("Missing required columns: %s", paste(missing_req, collapse = ", "))
    )
  }

  # Check 3: Completeness per column
  completeness <- sapply(names(data), function(col) {
    100 * sum(!is.na(data[[col]])) / n_rows
  })

  col_completeness <- data.frame(
    column = names(completeness),
    completeness_pct = as.numeric(completeness),
    n_missing = sapply(names(data), function(col) sum(is.na(data[[col]]))),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(completeness_pct)

  # Identify columns with excessive missingness
  poor_cols <- col_completeness %>%
    dplyr::filter(completeness_pct < (100 - max_missing_pct))

  if (nrow(poor_cols) > 0) {
    quality_issues <- c(
      quality_issues,
      sprintf(
        "%d columns exceed %.0f%% missing threshold: %s",
        nrow(poor_cols),
        max_missing_pct,
        paste(poor_cols$column[1:min(5, nrow(poor_cols))], collapse = ", ")
      )
    )
  }

  # Check 4: Overall completeness
  overall_completeness <- sum(!is.na(data)) / (nrow(data) * ncol(data))

  if (overall_completeness < 0.5) {
    quality_issues <- c(
      quality_issues,
      sprintf(
        "Overall completeness too low: %.1f%% (target: >=50%%)",
        overall_completeness * 100
      )
    )
  }

  # Determine pass/fail
  passes_quality <- length(quality_issues) == 0

  # Create result
  result <- list(
    passes_quality = passes_quality,
    n_rows = n_rows,
    n_cols = ncol(data),
    overall_completeness = overall_completeness,
    column_completeness = col_completeness,
    quality_issues = quality_issues
  )

  # Print results
  if (passes_quality) {
    message(sprintf(
      "[v] Quality check passed: %d rows x %d cols, %.1f%% complete",
      n_rows, ncol(data), overall_completeness * 100
    ))
  } else {
    message("[x] Quality issues detected:")
    for (issue in quality_issues) {
      message(sprintf("  - %s", issue))
    }
  }

  # Stop if requested
  if (!passes_quality && stop_on_failure) {
    stop("Data quality validation failed")
  }

  return(result)
}


#' Check Logical Consistency
#'
#' Validates logical relationships in the data (e.g., date sequences,
#' age consistency, valid value ranges).
#'
#' @param data Data frame
#' @param checks Character vector. Which checks to perform:
#'   - "date_sequence": admission < culture < outcome
#'   - "age_range": Age between 0-120
#'   - "age_dob_match": Age matches DOB
#'   - "outcome_consistency": Died patients have outcome date
#'   - "all": All checks (default)
#' @param stop_on_failure Logical. Stop on inconsistency. Default FALSE.
#'
#' @return List with consistency check results:
#'   - consistent: Logical
#'   - issues_found: Data frame of inconsistent rows
#'   - summary: Character vector of issue summaries
#'
#' @export
#'
#' @examples
#' \dontrun{
#' consistency <- check_logical_consistency(data, checks = "all")
#'
#' if (!consistency$consistent) {
#'   print(consistency$summary)
#'   View(consistency$issues_found)
#' }
#' }
check_logical_consistency <- function(data,
                                      checks = "all",
                                      stop_on_failure = FALSE) {
  # Expand "all" to all check types
  all_checks <- c("date_sequence", "age_range", "age_dob_match", "outcome_consistency")
  if ("all" %in% checks) {
    checks <- all_checks
  }

  issues <- list()
  summary_msgs <- character()

  # Check 1: Date sequence
  if ("date_sequence" %in% checks) {
    date_cols <- c("date_of_admission", "date_of_culture", "date_of_final_outcome")
    if (all(date_cols %in% names(data))) {
      date_issues <- data %>%
        dplyr::filter(
          !is.na(date_of_admission) & !is.na(date_of_culture) &
            (date_of_admission > date_of_culture |
              (!is.na(date_of_final_outcome) & date_of_culture > date_of_final_outcome))
        )

      if (nrow(date_issues) > 0) {
        issues$date_sequence <- date_issues
        summary_msgs <- c(
          summary_msgs,
          sprintf("Date sequence violations: %d rows", nrow(date_issues))
        )
      }
    }
  }

  # Check 2: Age range
  if ("age_range" %in% checks && "Age" %in% names(data)) {
    age_issues <- data %>%
      dplyr::filter(!is.na(Age) & (Age < 0 | Age > 120))

    if (nrow(age_issues) > 0) {
      issues$age_range <- age_issues
      summary_msgs <- c(
        summary_msgs,
        sprintf("Age out of range (0-120): %d rows", nrow(age_issues))
      )
    }
  }

  # Check 3: Age-DOB match
  if ("age_dob_match" %in% checks && all(c("Age", "DOB", "date_of_culture") %in% names(data))) {
    age_dob_issues <- data %>%
      dplyr::filter(
        !is.na(Age) & !is.na(DOB) & !is.na(date_of_culture)
      ) %>%
      dplyr::mutate(
        calculated_age = as.numeric(difftime(date_of_culture, DOB, units = "days")) / 365.25,
        age_diff = abs(Age - calculated_age)
      ) %>%
      dplyr::filter(age_diff > 2) # Allow 2-year tolerance

    if (nrow(age_dob_issues) > 0) {
      issues$age_dob_match <- age_dob_issues
      summary_msgs <- c(
        summary_msgs,
        sprintf("Age-DOB mismatch: %d rows (>2 year difference)", nrow(age_dob_issues))
      )
    }
  }

  # Check 4: Outcome consistency
  if ("outcome_consistency" %in% checks && all(c("final_outcome", "date_of_final_outcome") %in% names(data))) {
    outcome_issues <- data %>%
      dplyr::filter(
        final_outcome == "Died" & is.na(date_of_final_outcome)
      )

    if (nrow(outcome_issues) > 0) {
      issues$outcome_consistency <- outcome_issues
      summary_msgs <- c(
        summary_msgs,
        sprintf("Died patients without outcome date: %d rows", nrow(outcome_issues))
      )
    }
  }

  # Combine all issues
  is_consistent <- length(issues) == 0

  if (length(issues) > 0) {
    all_issues <- dplyr::bind_rows(lapply(names(issues), function(check_name) {
      issues[[check_name]] %>%
        dplyr::mutate(consistency_check = check_name)
    }))
  } else {
    all_issues <- data.frame()
  }

  # Create result
  result <- list(
    consistent = is_consistent,
    issues_found = all_issues,
    summary = summary_msgs,
    n_checks_performed = length(checks),
    n_issues = length(issues)
  )

  # Print results
  if (is_consistent) {
    message(sprintf(
      "[v] Logical consistency: All %d checks passed",
      length(checks)
    ))
  } else {
    message("[x] Logical inconsistencies found:")
    for (msg in summary_msgs) {
      message(sprintf("  - %s", msg))
    }
  }

  # Stop if requested
  if (!is_consistent && stop_on_failure) {
    stop("Logical consistency validation failed")
  }

  return(result)
}
