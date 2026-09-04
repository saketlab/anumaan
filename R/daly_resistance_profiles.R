# daly_resistance_profiles.R
#
# Estimates antimicrobial resistance profiles for two data sources that both
# feed the DALY (YLL/YLD) pipeline: Pathway 1 solves a simplex-constrained QP
# over aggregate/surveillance marginals (estimate_profiles_convex()); Pathway 2
# fits a Bayesian multivariate probit model to facility-level AST records.
#
# compute_marginal_resistance()/compute_pairwise_coresistance()/
# compute_resistance_profiles() is a second, independent implementation of the
# Pathway 1 QP, taking isolate-level input; only bootstrap_profiles_convex()
# uses it. Merging the two needs numerical equivalence tests first -- neither
# path has test coverage.


# Internal helpers

.null_default <- function(x, default) if (is.null(x)) default else x

.log_step <- function(log_list, step, n_in, n_out, message, status = "ok") {
  log_list[[length(log_list) + 1L]] <- tibble::tibble(
    step       = step,
    n_rows_in  = as.integer(n_in),
    n_rows_out = as.integer(n_out),
    message    = message,
    status     = status
  )
  log_list
}

.add_check <- function(results, check_name, status, message, n_affected = NA_integer_) {
  results[[length(results) + 1L]] <- tibble::tibble(
    check_name = check_name,
    status     = status,
    message    = message,
    n_affected = as.integer(n_affected)
  )
  results
}

# Resolve antibiotic-class and metadata columns from a wide data frame.
# Used by compute_marginals_from_data() and compute_pairwise_from_data().
.resolve_class_cols <- function(data_wide, col_map, panel_map, outcome_col = NULL) {
  pathogen_col <- .null_default(col_map$pathogen_col, "pathogen")
  geography_col <- .null_default(col_map$geography_col, "state")
  isolate_col <- .null_default(col_map$isolate_col, "isolate_id")

  meta_cols <- unique(c(
    isolate_col, pathogen_col,
    col_map$patient_col, col_map$date_col, geography_col,
    col_map$specimen_col, col_map$age_col, col_map$dob_col,
    col_map$location_col, outcome_col, "year"
  ))
  meta_cols <- meta_cols[!is.null(meta_cols) & meta_cols %in% names(data_wide)]
  class_cols <- setdiff(names(data_wide), meta_cols)
  class_cols <- class_cols[class_cols %in% unlist(panel_map)]

  list(
    pathogen_col  = pathogen_col,
    geography_col = geography_col,
    isolate_col   = isolate_col,
    meta_cols     = meta_cols,
    class_cols    = class_cols
  )
}

# Build the stratification column vector from stratify_by flags.
# Used by compute_marginals_from_data() and compute_pairwise_from_data().
.build_strat_cols <- function(stratify_by = NULL,
                              outcome_col = NULL,
                              geography_col = "state",
                              data) {
  strat_cols <- character(0)
  if ("geography" %in% stratify_by && geography_col %in% names(data)) {
    strat_cols <- c(strat_cols, geography_col)
  }
  if ("year" %in% stratify_by && "year" %in% names(data)) {
    strat_cols <- c(strat_cols, "year")
  }
  if (!is.null(outcome_col) && outcome_col %in% names(data)) {
    strat_cols <- c(strat_cols, outcome_col)
  }
  unique(strat_cols)
}

# Validate that paired arguments are both supplied or both NULL.
.check_paired_args <- function(col, val, col_arg, val_arg) {
  if (xor(is.null(col), is.null(val))) {
    stop(sprintf(
      "Both %s and %s must be provided together, or both NULL.",
      col_arg, val_arg
    ), call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate Inputs for Resistance Profile Estimation (Pathway 1)
#'
#' Checks that all mandatory columns exist, detects long vs wide format,
#' validates AST content, and verifies that columns required for any requested
#' stratification are present. Returns a tibble of check results and stops on
#' any failed mandatory check.
#'
#' \strong{Mandatory columns (always):}
#' \itemize{
#'   \item \code{isolate_col}   -- unique isolate identifier; the counting unit
#'     for all resistance statistics. One patient may have multiple isolates.
#'   \item \code{pathogen_col}  -- organism / pathogen name.
#'   \item \code{ast_col}       -- susceptibility result; must contain S / I / R.
#'   \item \code{patient_col}   -- patient identifier (retained for context, not counted).
#'   \item \code{date_col}      -- culture date; used to derive \code{year} when
#'     \code{"year"} is in \code{stratify_by}.
#'   \item \code{geography_col} -- geographic unit (state / district).
#'   \item \code{specimen_col}  -- specimen type (blood, urine, etc.).
#'   \item \code{age_col} OR \code{dob_col} -- at least one must be present.
#' }
#'
#' \strong{Mandatory in long format only:}
#' \itemize{
#'   \item \code{antibiotic_col} -- antibiotic / drug name column.
#' }
#'
#' \strong{Optional (warn if absent, do not stop):}
#' \itemize{
#'   \item \code{class_col}    -- antibiotic class; derived from \code{antibiotic_col}
#'     via WHO AWaRe reference if absent.
#'   \item \code{location_col} -- ward / ICU / OPD; retained for subgroup outputs.
#'   \item \code{outcome_col}  -- patient outcome; required for outcome stratification.
#' }
#'
#' @param data Data frame to validate.
#' @param col_map Named list mapping logical roles to actual column names.
#'   Defaults use standard package column names. Set any entry to \code{NULL}
#'   to skip that check (only valid for optional roles).
#' @param stratify_by Character vector or \code{NULL}. Requested stratification
#'   dimensions. Accepted values: \code{"geography"}, \code{"year"},
#'   \code{"outcome"}. Triggers additional column checks for each dimension.
#'   Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Outcome column name. Required
#'   when \code{"outcome"} is in \code{stratify_by}. Can also be supplied via
#'   \code{col_map$outcome_col}. Default \code{NULL}.
#'
#' @return A tibble with columns \code{check_name}, \code{status}
#'   (\code{"pass"} / \code{"warn"} / \code{"fail"}), \code{message}, and
#'   \code{n_affected}. The attribute \code{detected_format} (\code{"long"} or
#'   \code{"wide"}) is attached. Stops with an informative message if any
#'   mandatory check fails.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- validate_profile_inputs(
#'   data = my_ast_data,
#'   col_map = list(
#'     isolate_col = "isolate_id",
#'     pathogen_col = "organism",
#'     ast_col = "result",
#'     patient_col = "pid",
#'     date_col = "culture_date",
#'     geography_col = "state",
#'     specimen_col = "specimen",
#'     age_col = "age_years",
#'     dob_col = NULL,
#'     antibiotic_col = "drug_name",
#'     class_col = NULL,
#'     location_col = "ward",
#'     outcome_col = "discharge_status"
#'   ),
#'   stratify_by = c("geography", "year"),
#'   outcome_col = "discharge_status"
#' )
#' }
validate_profile_inputs <- function(
  data,
  col_map = list(
    isolate_col    = "isolate_id",
    pathogen_col   = "pathogen",
    ast_col        = "ast_value",
    patient_col    = "patient_id",
    date_col       = "date_of_culture",
    geography_col  = "state",
    specimen_col   = "specimen_type",
    age_col        = "age",
    dob_col        = "dob",
    antibiotic_col = "antibiotic_name",
    class_col      = "antibiotic_class",
    location_col   = NULL,
    outcome_col    = NULL
  ),
  stratify_by = NULL,
  outcome_col = NULL
) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (nrow(data) == 0L) stop("`data` has zero rows.")

  results <- list()

  # 1. Format detection
  abx_col <- .null_default(col_map$antibiotic_col, "antibiotic_name")
  iso_col <- .null_default(col_map$isolate_col, "isolate_id")

  has_abx_col <- abx_col %in% names(data)
  has_iso_col <- iso_col %in% names(data)

  if (has_abx_col && has_iso_col) {
    # Long format: each isolate appears in multiple rows (one per antibiotic)
    max_rows_per_iso <- max(table(data[[iso_col]]), na.rm = TRUE)
    detected_format <- if (max_rows_per_iso > 1L) "long" else "wide"
  } else if (has_abx_col) {
    detected_format <- "long"
  } else {
    detected_format <- "wide"
  }

  results <- .add_check(
    results, "format_detection", "pass",
    sprintf(
      "Data detected as '%s' format. %d rows x %d columns.",
      detected_format, nrow(data), ncol(data)
    )
  )

  # 2. Mandatory columns (all formats)
  mandatory_always <- c(
    isolate_col   = .null_default(col_map$isolate_col, "isolate_id"),
    pathogen_col  = .null_default(col_map$pathogen_col, "pathogen"),
    ast_col       = .null_default(col_map$ast_col, "ast_value"),
    patient_col   = .null_default(col_map$patient_col, "patient_id"),
    date_col      = .null_default(col_map$date_col, "date_of_culture"),
    geography_col = .null_default(col_map$geography_col, "state"),
    specimen_col  = .null_default(col_map$specimen_col, "specimen_type")
  )

  for (role in names(mandatory_always)) {
    col <- mandatory_always[[role]]
    if (is.na(col) || !col %in% names(data)) {
      results <- .add_check(
        results, paste0("mandatory_", role), "fail",
        sprintf("Mandatory column '%s' [role: %s] not found.", col, role)
      )
    } else {
      n_na <- sum(is.na(data[[col]]))
      results <- .add_check(
        results, paste0("mandatory_", role), "pass",
        sprintf("Column '%s' found. %d NA value(s).", col, n_na),
        n_na
      )
    }
  }

  # 3. Mandatory in long format: antibiotic_name column
  if (detected_format == "long") {
    if (!has_abx_col) {
      results <- .add_check(
        results, "mandatory_antibiotic_col", "fail",
        sprintf(
          "Long format detected but antibiotic column '%s' not found. ",
          abx_col
        )
      )
    } else {
      n_na <- sum(is.na(data[[abx_col]]))
      results <- .add_check(
        results, "mandatory_antibiotic_col", "pass",
        sprintf("Antibiotic column '%s' found. %d NA value(s).", abx_col, n_na),
        n_na
      )
    }
  }

  # 4. Age OR DOB -- at least one mandatory
  age_col <- col_map$age_col
  dob_col <- col_map$dob_col
  has_age <- !is.null(age_col) && age_col %in% names(data)
  has_dob <- !is.null(dob_col) && dob_col %in% names(data)

  if (!has_age && !has_dob) {
    results <- .add_check(
      results, "mandatory_age_or_dob", "fail",
      sprintf(
        "Neither age column ('%s') nor DOB column ('%s') found. At least one is mandatory.",
        .null_default(age_col, "age"), .null_default(dob_col, "dob")
      )
    )
  } else {
    found_desc <- paste(
      c(
        if (has_age) sprintf("age ('%s')", age_col),
        if (has_dob) sprintf("dob ('%s')", dob_col)
      ),
      collapse = " and "
    )
    results <- .add_check(
      results, "mandatory_age_or_dob", "pass",
      sprintf("Age/DOB check passed. Found: %s.", found_desc)
    )
  }

  # 5. Optional columns -- warn but do not stop
  eff_outcome_col <- .null_default(outcome_col, col_map$outcome_col)

  optional_checks <- list(
    class_col = list(
      col = col_map$class_col,
      msg_absent = "Antibiotic class column '%s' absent. Will derive classes from antibiotic names via WHO AWaRe reference."
    ),
    location_col = list(
      col = col_map$location_col,
      msg_absent = "Location column '%s' absent. Ward/ICU/OPD subgroup outputs will not be available."
    ),
    outcome_col = list(
      col = eff_outcome_col,
      msg_absent = "Outcome column '%s' absent. Outcome stratification will not be possible."
    )
  )

  for (role in names(optional_checks)) {
    col <- optional_checks[[role]]$col
    if (is.null(col)) next
    if (!col %in% names(data)) {
      results <- .add_check(
        results, paste0("optional_", role), "warn",
        sprintf(optional_checks[[role]]$msg_absent, col)
      )
    } else {
      n_na <- sum(is.na(data[[col]]))
      results <- .add_check(
        results, paste0("optional_", role), "pass",
        sprintf("Optional column '%s' found. %d NA value(s).", col, n_na),
        n_na
      )
    }
  }

  # 6. AST value content check
  ast_col_name <- .null_default(col_map$ast_col, "ast_value")
  if (ast_col_name %in% names(data)) {
    ast_raw <- toupper(trimws(as.character(data[[ast_col_name]])))
    valid_sir <- c("S", "I", "R")
    n_valid <- sum(ast_raw %in% valid_sir, na.rm = TRUE)
    n_missing <- sum(is.na(data[[ast_col_name]]))
    n_invalid <- nrow(data) - n_valid - n_missing
    pct_valid <- 100 * n_valid / nrow(data)

    check_status <- if (pct_valid < 50) "fail" else if (n_invalid > 0) "warn" else "pass"
    check_msg <- sprintf(
      "AST column '%s': %.1f%% valid S/I/R (%d valid, %d missing, %d non-standard).",
      ast_col_name, pct_valid, n_valid, n_missing, n_invalid
    )
    if (n_invalid > 0 && check_status != "fail") {
      check_msg <- paste(check_msg, "Non-standard values will be treated as NA.")
    }
    results <- .add_check(results, "ast_value_content", check_status, check_msg, n_invalid)
  }

  # 7. Isolate ID uniqueness (wide format only)
  if (detected_format == "wide" && has_iso_col) {
    n_dup <- sum(duplicated(data[[iso_col]], incomparables = NA))
    if (n_dup > 0L) {
      results <- .add_check(
        results, "isolate_id_uniqueness", "warn",
        sprintf(
          "%d duplicate isolate IDs in '%s' (wide format expects one row per isolate).",
          n_dup, iso_col
        ),
        n_dup
      )
    } else {
      results <- .add_check(
        results, "isolate_id_uniqueness", "pass",
        sprintf("All isolate IDs in '%s' are unique.", iso_col)
      )
    }
  }

  # 8. Stratification-specific checks
  # outcome_col is a separate split variable, not a stratify_by dimension.
  valid_strat <- c("geography", "year")
  if (!is.null(stratify_by)) {
    unknown_strat <- setdiff(stratify_by, valid_strat)
    if (length(unknown_strat) > 0L) {
      warning(sprintf(
        "Unknown stratify_by value(s) ignored: %s. Valid: %s.",
        paste(unknown_strat, collapse = ", "),
        paste(valid_strat, collapse = ", ")
      ), call. = FALSE)
    }
  }

  strat_checks <- list(
    geography = list(
      col      = .null_default(col_map$geography_col, "state"),
      fail_msg = "Geography stratification requested but geography column '%s' not found."
    ),
    year = list(
      col      = .null_default(col_map$date_col, "date_of_culture"),
      fail_msg = "Year stratification requested but date column '%s' not found. Year is derived from the date column."
    )
  )

  # Separate check: outcome column required when outcome_col is supplied
  if (!is.null(eff_outcome_col)) {
    if (!eff_outcome_col %in% names(data)) {
      results <- .add_check(
        results, "outcome_col", "fail",
        sprintf("outcome_col '%s' requested but not found in data.", eff_outcome_col)
      )
    } else {
      n_na_out <- sum(is.na(data[[eff_outcome_col]]))
      results <- .add_check(
        results, "outcome_col", "pass",
        sprintf("Outcome column '%s' found. %d NA value(s).", eff_outcome_col, n_na_out), n_na_out
      )
    }
  }

  for (dim in intersect(stratify_by, valid_strat)) {
    col <- strat_checks[[dim]]$col
    if (is.null(col) || !col %in% names(data)) {
      results <- .add_check(
        results, paste0("stratify_", dim), "fail",
        sprintf(strat_checks[[dim]]$fail_msg, .null_default(col, paste0(dim, "_col")))
      )
    } else {
      n_strata <- dplyr::n_distinct(data[[col]], na.rm = TRUE)
      extra <- if (dim == "year") {
        " Year will be extracted as integer from this column."
      } else {
        ""
      }
      results <- .add_check(
        results, paste0("stratify_", dim), "pass",
        sprintf(
          "Stratification by %s: column '%s' found with %d unique value(s).%s",
          dim, col, n_strata, extra
        ),
        n_strata
      )
    }
  }

  # 9. Collate, report, and stop/warn
  results_tbl <- dplyr::bind_rows(results)
  n_fail <- sum(results_tbl$status == "fail")
  n_warn <- sum(results_tbl$status == "warn")
  n_pass <- sum(results_tbl$status == "pass")

  message(sprintf(
    "\n[validate_profile_inputs] %d passed | %d warnings | %d failed",
    n_pass, n_warn, n_fail
  ))

  if (n_warn > 0L) {
    warn_msgs <- results_tbl$message[results_tbl$status == "warn"]
    for (msg in warn_msgs) message(sprintf("  [!] %s", msg))
  }

  if (n_fail > 0L) {
    fail_msgs <- results_tbl$message[results_tbl$status == "fail"]
    stop(
      sprintf(
        "%d mandatory check(s) failed:\n%s",
        n_fail, paste(sprintf("  - %s", fail_msgs), collapse = "\n")
      ),
      call. = FALSE
    )
  }

  attr(results_tbl, "detected_format") <- detected_format
  invisible(results_tbl)
}


#' Preprocess AST Data for Resistance Profile Estimation (Pathway 1)
#'
#' Orchestrates the full preprocessing pipeline required before marginal
#' resistance computation and convex optimisation. Calls existing
#' \code{prep_*} package functions in sequence and returns a wide tibble
#' (one row per isolate) with antibiotic-class columns ready for
#' \code{compute_marginals_from_data()}.
#'
#' \strong{Pipeline steps (in order):}
#' \enumerate{
#'   \item Input validation via \code{validate_profile_inputs()}.
#'   \item Wide -> long conversion if wide format detected
#'     (\code{prep_pivot_ast_wide_to_long()}).
#'   \item AST value harmonisation -- clean messy strings, resolve I values
#'     (\code{prep_harmonize_ast()}).
#'   \item Antibiotic name standardisation if no class column
#'     (\code{prep_standardize_antibiotics()}).
#'   \item Class assignment if still absent
#'     (\code{prep_classify_antibiotic_class()}).
#'   \item Class-level collapse: any R in class => class = R
#'     (\code{prep_collapse_class_level()}). Counting unit is \code{isolate_id}.
#'   \item Year derivation from \code{date_col} when \code{"year"} is in
#'     \code{stratify_by}.
#'   \item Organism-specific panel filter: retain only user-specified classes
#'     per pathogen.
#'   \item Long -> wide pivot with antibiotic-class columns
#'     (\code{prep_create_wide_ast_matrix()}).
#' }
#'
#' \strong{Counting unit:} All statistics are computed per unique
#' \code{isolate_id}. A patient with multiple isolates contributes one count
#' per isolate, not one count per patient.
#'
#' @param data Data frame. Raw AST data in long or wide format.
#' @param col_map Named list mapping logical roles to actual column names.
#'   Same structure as in \code{validate_profile_inputs()}. Default names
#'   follow standard package conventions.
#' @param panel_map Named list. \strong{Mandatory.} Maps each pathogen name to
#'   a character vector of antibiotic class names to include in profile
#'   estimation. Classes not in the panel are excluded and logged with reason
#'   code \code{"not_in_user_panel"}. Example:
#'   \code{list("Klebsiella pneumoniae" = c("3GC", "Carbapenems", "Fluoroquinolones"))}.
#' @param stratify_by Character vector or \code{NULL}. Stratification dimensions
#'   to carry through to output. Valid values: \code{"geography"},
#'   \code{"year"}, \code{"outcome"}. Each requested dimension adds the
#'   corresponding column to the wide output. Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Column name for patient outcome
#'   (e.g. \code{"final_outcome"}). Required when \code{"outcome"} is in
#'   \code{stratify_by}. Default \code{NULL}.
#' @param who_table Data frame or \code{NULL}. Custom WHO AWaRe classification
#'   table passed to \code{prep_standardize_antibiotics()}. If \code{NULL},
#'   uses the built-in \code{inst/extdata/WHO_aware_class.csv}. Default
#'   \code{NULL}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{data_wide}}{Tibble. One row per \code{isolate_id}. Columns:
#'     all metadata columns + one column per antibiotic class (values S / R / NA).}
#'   \item{\code{preprocessing_log}}{Tibble. One row per pipeline step with
#'     columns \code{step}, \code{n_rows_in}, \code{n_rows_out}, \code{message},
#'     \code{status}.}
#'   \item{\code{panel_exclusions}}{Tibble. All organism-class combinations
#'     removed by the panel filter, with \code{reason_code} and
#'     \code{reason_text}.}
#'   \item{\code{col_map_resolved}}{Named list. Updated \code{col_map} with
#'     any column names derived during preprocessing (e.g. the harmonised AST
#'     column name, the resolved class column name).}
#'   \item{\code{detected_format}}{Character. \code{"long"} or \code{"wide"}.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- preprocess_for_profiles(
#'   data = raw_ast,
#'   col_map = list(
#'     isolate_col = "isolate_id",
#'     pathogen_col = "organism",
#'     ast_col = "result",
#'     patient_col = "pid",
#'     date_col = "culture_date",
#'     geography_col = "state",
#'     specimen_col = "specimen",
#'     age_col = "age",
#'     antibiotic_col = "drug_name"
#'   ),
#'   panel_map = list(
#'     "Klebsiella pneumoniae" = c("3GC", "Carbapenems", "Fluoroquinolones"),
#'     "Escherichia coli"      = c("3GC", "Fluoroquinolones", "Aminoglycosides")
#'   ),
#'   stratify_by = c("geography", "year"),
#'   outcome_col = "final_outcome"
#' )
#'
#' result$data_wide # ready for compute_marginals_from_data()
#' result$preprocessing_log # step-by-step row counts
#' result$panel_exclusions # classes dropped per organism
#' }
preprocess_for_profiles <- function(
  data,
  col_map = list(
    isolate_col    = "isolate_id",
    pathogen_col   = "pathogen",
    ast_col        = "ast_value",
    patient_col    = "patient_id",
    date_col       = "date_of_culture",
    geography_col  = "state",
    specimen_col   = "specimen_type",
    age_col        = "age",
    dob_col        = "dob",
    antibiotic_col = "antibiotic_name",
    class_col      = "antibiotic_class",
    location_col   = NULL,
    outcome_col    = NULL
  ),
  panel_map,
  stratify_by = NULL,
  outcome_col = NULL,
  who_table = NULL
) {
  # Argument checks
  if (missing(panel_map) || is.null(panel_map) || length(panel_map) == 0L) {
    stop(
      "`panel_map` is mandatory. Supply a named list mapping each pathogen ",
      "to its antibiotic class panel, e.g.:\n",
      "  list('Klebsiella pneumoniae' = c('3GC', 'Carbapenems', 'Fluoroquinolones'))",
      call. = FALSE
    )
  }
  if (!is.list(panel_map) || is.null(names(panel_map))) {
    stop("`panel_map` must be a named list (pathogen names as names, class vectors as values).",
      call. = FALSE
    )
  }

  eff_outcome_col <- .null_default(outcome_col, col_map$outcome_col)

  # Resolve column names with defaults
  isolate_col <- .null_default(col_map$isolate_col, "isolate_id")
  pathogen_col <- .null_default(col_map$pathogen_col, "pathogen")
  ast_col <- .null_default(col_map$ast_col, "ast_value")
  patient_col <- .null_default(col_map$patient_col, "patient_id")
  date_col <- .null_default(col_map$date_col, "date_of_culture")
  geography_col <- .null_default(col_map$geography_col, "state")
  specimen_col <- .null_default(col_map$specimen_col, "specimen_type")
  age_col <- col_map$age_col
  dob_col <- col_map$dob_col
  antibiotic_col <- .null_default(col_map$antibiotic_col, "antibiotic_name")
  class_col <- col_map$class_col
  location_col <- col_map$location_col

  plog <- list()
  exclusions <- list()

  # Step 0: Input validation
  message("\n[preprocess_for_profiles] Step 0: Validating inputs...")
  validation_report <- validate_profile_inputs(
    data        = data,
    col_map     = col_map,
    stratify_by = stratify_by,
    outcome_col = eff_outcome_col
  )
  detected_format <- attr(validation_report, "detected_format")
  plog <- .log_step(
    plog, "0_validate", nrow(data), nrow(data),
    sprintf("Validation passed. Format: %s.", detected_format)
  )

  # Step 1: Wide -> long conversion
  message("\n[preprocess_for_profiles] Step 1: Format handling...")
  n_in <- nrow(data)

  if (detected_format == "wide") {
    # Identify antibiotic / class columns: all columns not in the known
    # metadata set are candidates for antibiotic columns.
    known_meta <- unique(c(
      isolate_col, pathogen_col, patient_col, date_col,
      geography_col, specimen_col,
      age_col, dob_col, location_col, eff_outcome_col
    ))
    known_meta <- known_meta[!is.null(known_meta) & known_meta %in% names(data)]
    candidate_abx <- setdiff(names(data), known_meta)

    if (length(candidate_abx) == 0L) {
      stop("Wide format detected but no antibiotic/class columns found after excluding metadata columns.",
        call. = FALSE
      )
    }
    message(sprintf("  Wide -> long: pivoting %d antibiotic column(s).", length(candidate_abx)))

    data <- prep_pivot_ast_wide_to_long(
      data = data,
      antibiotic_cols = candidate_abx,
      id_cols = known_meta,
      antibiotic_name_col = antibiotic_col,
      antibiotic_value_col = ast_col,
      remove_missing = FALSE # keep NA rows; will be handled in harmonisation
    )
    plog <- .log_step(
      plog, "1_wide_to_long", n_in, nrow(data),
      sprintf("Wide -> long: pivoted %d antibiotic column(s).", length(candidate_abx))
    )
  } else {
    plog <- .log_step(plog, "1_format", n_in, nrow(data), "Long format confirmed. No pivot needed.")
  }

  # Step 2: AST harmonisation
  message("\n[preprocess_for_profiles] Step 2: Harmonising AST values...")
  n_in <- nrow(data)

  # prep_harmonize_ast expects ast_col to exist; rename temporarily if needed
  abx_std_col <- if (antibiotic_col %in% names(data)) antibiotic_col else "antibiotic_name"
  data <- prep_harmonize_ast(
    data           = data,
    antibiotic_col = abx_std_col,
    ast_col        = ast_col,
    colistin_to_s  = TRUE,
    others_to_r    = TRUE
  )
  # After harmonisation, the clean values are in ast_value_harmonized
  harmonised_ast_col <- "ast_value_harmonized"
  plog <- .log_step(
    plog, "2_harmonise_ast", n_in, nrow(data),
    sprintf("AST values harmonised into '%s'. I values recoded.", harmonised_ast_col)
  )

  # Step 3: Antibiotic name standardisation + class assignment
  message("\n[preprocess_for_profiles] Step 3: Antibiotic standardisation and class assignment...")
  n_in <- nrow(data)

  has_class_col <- !is.null(class_col) && class_col %in% names(data)

  if (!has_class_col) {
    message("  No antibiotic class column found. Deriving from antibiotic names via WHO AWaRe reference.")

    if (antibiotic_col %in% names(data)) {
      data <- prep_standardize_antibiotics(
        data           = data,
        antibiotic_col = antibiotic_col,
        who_table      = who_table,
        add_class      = TRUE,
        add_aware      = FALSE
      )
      # prep_standardize_antibiotics adds antibiotic_normalized and antibiotic_class
      class_col <- "antibiotic_class"
      plog <- .log_step(
        plog, "3a_standardize_abx", n_in, nrow(data),
        "Antibiotic names standardised via WHO AWaRe; antibiotic_class derived."
      )
    } else {
      stop(sprintf(
        "No antibiotic class column ('%s') and no antibiotic name column ('%s') found. Cannot derive classes.",
        .null_default(col_map$class_col, "antibiotic_class"), antibiotic_col
      ), call. = FALSE)
    }
  } else {
    plog <- .log_step(
      plog, "3a_class_col", n_in, nrow(data),
      sprintf("Antibiotic class column '%s' already present. Standardisation skipped.", class_col)
    )
  }

  # Ensure class_col is set
  if (!"antibiotic_class" %in% names(data) && class_col %in% names(data)) {
    data$antibiotic_class <- data[[class_col]]
    class_col <- "antibiotic_class"
  }

  plog <- .log_step(
    plog, "3b_class_assignment", n_in, nrow(data),
    sprintf("Antibiotic class column resolved to '%s'.", class_col)
  )

  # Step 4: Class-level collapse
  #   Rule: any R in class => class = R   (ertapenem NA + imipenem S + meropenem R -> Carbapenems R)
  #   Counting unit: isolate_id (not patient_id)
  message("\n[preprocess_for_profiles] Step 4: Collapsing to class level (any R -> class R)...")
  n_in <- nrow(data)

  # prep_collapse_class_level() hardcodes a reference to "antibiotic_normalized"
  # for the drugs_tested summary column. Ensure it exists so the call never
  # errors when the user supplied antibiotic_class directly (skipping Step 3).
  if (!"antibiotic_normalized" %in% names(data)) {
    data$antibiotic_normalized <- if (antibiotic_col %in% names(data)) {
      data[[antibiotic_col]]
    } else {
      NA_character_
    }
  }

  # Determine which metadata columns to carry through collapse
  # (prep_collapse_class_level uses extra_cols for passthrough)
  carry_cols <- unique(c(
    patient_col, date_col, geography_col, specimen_col,
    age_col, dob_col, location_col, eff_outcome_col,
    if ("antibiotic_normalized" %in% names(data)) "antibiotic_normalized"
  ))
  carry_cols <- carry_cols[!is.null(carry_cols) & carry_cols %in% names(data)]
  # Remove columns that are already part of the grouping keys in collapse
  grouping_keys <- c(isolate_col, pathogen_col, class_col)
  carry_cols <- setdiff(carry_cols, grouping_keys)

  data_collapsed <- prep_collapse_class_level(
    data = data,
    event_col = isolate_col, # isolate_id is the counting unit
    organism_col = pathogen_col,
    class_col = class_col,
    susceptibility_col = harmonised_ast_col,
    extra_cols = if (length(carry_cols) > 0L) carry_cols else NULL
  )

  plog <- .log_step(
    plog, "4_collapse_class", n_in, nrow(data_collapsed),
    sprintf(
      "Class-level collapse: %d drug rows -> %d isolate-class rows. Unit: %s.",
      n_in, nrow(data_collapsed), isolate_col
    )
  )
  data <- data_collapsed

  # Step 4b: Isolate deduplication
  #   Same isolate_id x antibiotic combination can appear more than once
  #   (e.g. duplicate data entry, multi-join artefact). Keep the worst
  #   phenotype per isolate x drug (R > I > S > NA).
  message("\n[preprocess_for_profiles] Step 4b: Deduplicating isolate x drug rows...")
  n_in <- nrow(data)

  abx_dedup_col <- if ("antibiotic_normalized" %in% names(data)) {
    "antibiotic_normalized"
  } else if (antibiotic_col %in% names(data)) {
    antibiotic_col
  } else {
    NULL
  }

  if (!is.null(abx_dedup_col)) {
    sir_rank <- c("R" = 3L, "I" = 2L, "S" = 1L)
    data <- data %>%
      dplyr::mutate(
        .sir_rank = dplyr::coalesce(
          sir_rank[.data[[harmonised_ast_col]]],
          0L
        )
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(isolate_col, abx_dedup_col)))) %>%
      dplyr::arrange(dplyr::desc(.sir_rank), .by_group = TRUE) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.sir_rank)
    n_removed <- n_in - nrow(data)
    plog <- .log_step(
      plog, "4b_dedup", n_in, nrow(data),
      sprintf(
        "Dedup: %d duplicate isolate-drug row(s) removed (kept worst phenotype).",
        n_removed
      )
    )
  } else {
    plog <- .log_step(plog, "4b_dedup", n_in, nrow(data),
      "Dedup skipped: no antibiotic name column identified.",
      status = "warn"
    )
  }

  # Step 5: Year derivation -- always performed when date_col is present.
  #   Year is always added to the output so downstream stratification
  #   by year does not require re-running preprocessing.
  message("\n[preprocess_for_profiles] Step 5: Deriving year from date column...")
  n_in <- nrow(data)

  if (date_col %in% names(data)) {
    parsed_dates <- suppressWarnings(as.Date(as.character(data[[date_col]])))
    data$year <- as.integer(format(parsed_dates, "%Y"))
    n_missing_yr <- sum(is.na(data$year))
    plog <- .log_step(
      plog, "5_year_derivation", n_in, nrow(data),
      sprintf("Year derived from '%s'. %d missing year value(s).", date_col, n_missing_yr)
    )
  } else {
    warning(sprintf("Date column '%s' not found. 'year' column not added.", date_col),
      call. = FALSE
    )
    plog <- .log_step(plog, "5_year_derivation", n_in, nrow(data),
      sprintf("Skipped: date column '%s' not present.", date_col),
      status = "warn"
    )
  }

  # Step 6: Panel filter
  #   For each pathogen, retain only classes listed in panel_map.
  #   Everything else is logged as excluded with reason_code.
  message("\n[preprocess_for_profiles] Step 6: Applying organism-specific class panels...")
  n_in <- nrow(data)

  panel_tbl <- dplyr::bind_rows(lapply(names(panel_map), function(org) {
    tibble::tibble(
      !!pathogen_col := org,
      antibiotic_class = panel_map[[org]]
    )
  }))

  # Validate that panel class names actually appear in the collapsed data.
  # Misspelled class names produce silently empty panels otherwise.
  classes_in_data <- unique(data[["antibiotic_class"]])
  for (org in names(panel_map)) {
    bad_classes <- setdiff(panel_map[[org]], classes_in_data)
    if (length(bad_classes) > 0L) {
      warning(sprintf(
        "panel_map['%s']: %d class name(s) not found in data and will produce no rows: %s",
        org, length(bad_classes), paste(bad_classes, collapse = ", ")
      ), call. = FALSE)
    }
  }

  # Pathogens not in panel_map at all -- flag but keep (no panel defined for them)
  pathogens_in_data <- unique(data[[pathogen_col]])
  pathogens_no_panel <- setdiff(pathogens_in_data, names(panel_map))
  if (length(pathogens_no_panel) > 0L) {
    warning(sprintf(
      "%d pathogen(s) in data have no entry in panel_map and will be excluded from profiling: %s",
      length(pathogens_no_panel),
      paste(head(pathogens_no_panel, 5L), collapse = ", ")
    ), call. = FALSE)
    exclusions <- c(exclusions, list(tibble::tibble(
      !!pathogen_col := pathogens_no_panel,
      antibiotic_class = NA_character_,
      reason_code = "no_panel_defined",
      reason_text = "Pathogen has no entry in panel_map; excluded from profiling."
    )))
  }

  # Classes present in data for panel pathogens but not in their panel
  data_panel_pathogens <- data[data[[pathogen_col]] %in% names(panel_map), ]
  excluded_rows <- dplyr::anti_join(
    data_panel_pathogens[, c(pathogen_col, "antibiotic_class")],
    panel_tbl,
    by = c(pathogen_col, "antibiotic_class")
  ) %>%
    dplyr::distinct()

  if (nrow(excluded_rows) > 0L) {
    exclusions <- c(exclusions, list(
      excluded_rows %>%
        dplyr::mutate(
          reason_code = "not_in_user_panel",
          reason_text = "Antibiotic class is not in the user-specified panel for this pathogen."
        )
    ))
    message(sprintf(
      "  Panel filter: %d pathogen-class combination(s) excluded (not in panel).",
      nrow(excluded_rows)
    ))
  }

  # Apply filter: keep only panel-matching rows
  data_filtered <- dplyr::semi_join(data, panel_tbl, by = c(pathogen_col, "antibiotic_class"))

  plog <- .log_step(
    plog, "6_panel_filter", n_in, nrow(data_filtered),
    sprintf(
      "Panel filter: %d -> %d rows. %d exclusion(s) logged.",
      n_in, nrow(data_filtered), nrow(excluded_rows)
    )
  )
  data <- data_filtered

  # Step 7: Long -> wide with antibiotic classes as columns
  #   One row per isolate_id. Class columns: S / R / NA.
  message("\n[preprocess_for_profiles] Step 7: Pivoting to wide format (classes as columns)...")
  n_in <- nrow(data)

  # Assemble the metadata columns to carry into the wide output
  wide_keep_cols <- unique(c(
    patient_col, pathogen_col, date_col, geography_col, specimen_col,
    age_col, dob_col, location_col, eff_outcome_col,
    if ("year" %in% names(data)) "year"
  ))
  wide_keep_cols <- wide_keep_cols[!is.null(wide_keep_cols) & wide_keep_cols %in% names(data)]

  data_wide <- prep_create_wide_ast_matrix(
    data = data,
    event_col = isolate_col,
    antibiotic_col = "antibiotic_class",
    susceptibility_col = "class_resistance",
    prefix = "", # class names become column names directly
    keep_cols = wide_keep_cols
  )

  # Sanitise column names (prep_create_wide_ast_matrix already does this, but be explicit)
  names(data_wide) <- gsub("[^A-Za-z0-9_]", "_", names(data_wide))
  names(data_wide) <- gsub("_{2,}", "_", names(data_wide))

  plog <- .log_step(
    plog, "7_wide_pivot", n_in, nrow(data_wide),
    sprintf(
      "Long -> wide: %d isolates x %d class column(s).",
      nrow(data_wide),
      ncol(data_wide) - length(wide_keep_cols) - 1L
    )
  )

  # Collate outputs
  col_map_resolved <- col_map
  col_map_resolved$ast_col_harmonised <- harmonised_ast_col
  col_map_resolved$class_col <- "antibiotic_class"
  col_map_resolved$outcome_col <- eff_outcome_col

  exclusions_tbl <- if (length(exclusions) > 0L) {
    dplyr::bind_rows(exclusions)
  } else {
    tibble::tibble(
      pathogen         = character(),
      antibiotic_class = character(),
      reason_code      = character(),
      reason_text      = character()
    )
  }

  log_tbl <- dplyr::bind_rows(plog)
  message(sprintf(
    "\n[preprocess_for_profiles] Complete. Output: %d isolate(s) x %d column(s). Steps logged: %d.",
    nrow(data_wide), ncol(data_wide), nrow(log_tbl)
  ))

  list(
    data_wide = tibble::as_tibble(data_wide),
    preprocessing_log = log_tbl,
    panel_exclusions = exclusions_tbl,
    col_map_resolved = col_map_resolved,
    detected_format = detected_format
  )
}


#' Formally Check Resistance Profile Probability Constraints
#'
#' Validates that the profile probability distributions returned by
#' \code{compute_resistance_profiles()} satisfy all mathematical requirements:
#' non-negativity, sum-to-one, and reconstruction of the input marginal and
#' pairwise resistance rates within a specified tolerance.
#'
#' The function uses two complementary sources of evidence:
#' \enumerate{
#'   \item \strong{Constraint residuals} already stored in
#'     \code{profiles_output[[pathogen]]$constraint_residuals} -- these are
#'     \eqn{M \hat{p} - v} computed at solve time and are available without
#'     re-supplying the original rates.
#'   \item \strong{Direct reconstruction} from the binary class-indicator
#'     columns in the profiles data frame when \code{marginals} and/or
#'     \code{pairwise} are supplied -- independently verifies the residuals
#'     and catches any post-solve modifications to the probability vector.
#' }
#'
#' @param profiles_output Named list returned by
#'   \code{compute_resistance_profiles()}. Each element must contain
#'   \code{$profiles} (data frame with \code{probability} and binary class
#'   columns) and optionally \code{$constraint_residuals}.
#' @param marginals Data frame or \code{NULL}. Original marginal resistance
#'   rates used as QP constraints. When provided, reconstructed marginals are
#'   compared against these values. Must have columns \code{pathogen_col},
#'   \code{class_col}, and \code{rate_col}. Default \code{NULL}.
#' @param pairwise Data frame or \code{NULL}. Original pairwise co-resistance
#'   rates. When provided, reconstructed pairwise values are also checked.
#'   Must have columns \code{pathogen_col}, \code{class1_col},
#'   \code{class2_col}, and \code{pairwise_rate_col}. Default \code{NULL}.
#' @param tolerance Numeric. Maximum absolute residual considered a pass.
#'   Default \code{1e-6}.
#' @param pathogen_col Character. Pathogen column in \code{marginals} /
#'   \code{pairwise}. Default \code{"pathogen"}.
#' @param class_col Character. Class column in \code{marginals}.
#'   Default \code{"antibiotic_class"}.
#' @param rate_col Character. Marginal rate column in \code{marginals}.
#'   Default \code{"marginal_resistance"}.
#' @param class1_col Character. First class column in \code{pairwise}.
#'   Default \code{"antibiotic_class_1"}.
#' @param class2_col Character. Second class column in \code{pairwise}.
#'   Default \code{"antibiotic_class_2"}.
#' @param pairwise_rate_col Character. Pairwise rate column.
#'   Default \code{"pairwise_resistance_prevalence"}.
#'
#' @return A tibble with one row per constraint per pathogen:
#' \describe{
#'   \item{\code{pathogen}}{Pathogen name.}
#'   \item{\code{constraint_type}}{One of \code{"nonneg"}, \code{"sum_to_one"},
#'     \code{"marginal"}, \code{"pairwise"}.}
#'   \item{\code{constraint_name}}{Specific label (e.g. \code{"marg_3GC"},
#'     \code{"pair_3GC_Carbapenems"}).}
#'   \item{\code{target}}{The constraint target value (\code{NA} for
#'     \code{nonneg}).}
#'   \item{\code{reconstructed}}{The value implied by \eqn{\hat{p}}.}
#'   \item{\code{abs_residual}}{Absolute difference between target and
#'     reconstructed.}
#'   \item{\code{pass}}{Logical. \code{TRUE} if \code{abs_residual < tolerance}
#'     or the structural check passed.}
#'   \item{\code{source}}{Either \code{"stored_residuals"} (from solve-time
#'     values) or \code{"recomputed"} (from direct reconstruction).}
#' }
#' Also prints a per-pathogen summary of passed and failed checks.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(amr_clean)
#' co_res <- compute_pairwise_coresistance(marg)
#' rp <- compute_resistance_profiles(marg, co_res)
#'
#' # Using stored constraint residuals only
#' checks <- check_profile_constraints(rp)
#'
#' # Full re-verification with original rates
#' checks <- check_profile_constraints(
#'   rp,
#'   marginals = marg$marginal,
#'   pairwise  = NULL,
#'   tolerance = 1e-5
#' )
#'
#' # Inspect failures
#' checks[!checks$pass, ]
#' }
check_profile_constraints <- function(
  profiles_output,
  marginals = NULL,
  pairwise = NULL,
  tolerance = 1e-6,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rate_col = "marginal_resistance",
  class1_col = "antibiotic_class_1",
  class2_col = "antibiotic_class_2",
  pairwise_rate_col = "pairwise_resistance_prevalence"
) {
  if (!is.list(profiles_output) || length(profiles_output) == 0L) {
    stop("`profiles_output` must be the non-empty named list from compute_resistance_profiles().",
      call. = FALSE
    )
  }

  # Non-class columns that should not be treated as class indicator columns
  non_class_cols <- c(
    "profile", "probability", "RR_LOS_profile", "dominant_class",
    "CI_lower_profile", "CI_upper_profile", "numerator", "PAF_LOS",
    "denominator", "numerator_R_kd", "R_kd"
  )

  all_rows <- list()

  for (path in names(profiles_output)) {
    entry <- profiles_output[[path]]

    if (!is.list(entry) || !"profiles" %in% names(entry)) {
      warning(sprintf("'%s': no $profiles element found; skipping.", path), call. = FALSE)
      next
    }

    prof_df <- entry$profiles
    p_hat <- prof_df$probability
    n_prof <- length(p_hat)
    classes <- if ("classes" %in% names(entry)) {
      entry$classes
    } else {
      setdiff(names(prof_df), non_class_cols)
    }

    .row <- function(type, name, target, reconstructed, source) {
      abs_res <- if (is.na(target)) NA_real_ else abs(reconstructed - target)
      pass <- if (is.na(target)) reconstructed >= -tolerance else abs_res < tolerance
      tibble::tibble(
        pathogen        = path,
        constraint_type = type,
        constraint_name = name,
        target          = target,
        reconstructed   = reconstructed,
        abs_residual    = abs_res,
        pass            = pass,
        source          = source
      )
    }

    # Check 1: Non-negativity -- all p_hat >= 0
    min_p <- min(p_hat, na.rm = TRUE)
    all_rows <- c(all_rows, list(.row(
      "nonneg", "min_probability",
      target = NA_real_,
      reconstructed = min_p,
      source = "recomputed"
    )))

    # Check 2: Sum-to-one
    p_sum <- sum(p_hat, na.rm = TRUE)
    all_rows <- c(all_rows, list(.row(
      "sum_to_one", "sum_probability",
      target = 1.0,
      reconstructed = p_sum,
      source = "recomputed"
    )))

    # Check 3a: Marginal constraints from stored constraint_residuals
    if ("constraint_residuals" %in% names(entry)) {
      resid <- entry$constraint_residuals
      targets <- if ("constraint_targets" %in% names(entry)) {
        entry$constraint_targets
      } else {
        setNames(rep(NA_real_, length(resid)), names(resid))
      }

      marg_nms <- names(resid)[grepl("^marg_", names(resid))]
      pair_nms <- names(resid)[grepl("^pair_", names(resid))]

      for (nm in marg_nms) {
        tgt <- targets[[nm]]
        recon <- if (!is.na(tgt)) tgt + resid[[nm]] else NA_real_
        all_rows <- c(all_rows, list(tibble::tibble(
          pathogen        = path,
          constraint_type = "marginal",
          constraint_name = nm,
          target          = tgt,
          reconstructed   = recon,
          abs_residual    = abs(resid[[nm]]),
          pass            = abs(resid[[nm]]) < tolerance,
          source          = "stored_residuals"
        )))
      }
      for (nm in pair_nms) {
        tgt <- targets[[nm]]
        recon <- if (!is.na(tgt)) tgt + resid[[nm]] else NA_real_
        all_rows <- c(all_rows, list(tibble::tibble(
          pathogen        = path,
          constraint_type = "pairwise",
          constraint_name = nm,
          target          = tgt,
          reconstructed   = recon,
          abs_residual    = abs(resid[[nm]]),
          pass            = abs(resid[[nm]]) < tolerance,
          source          = "stored_residuals"
        )))
      }
    }

    # Check 3b: Marginal constraints recomputed from supplied marginals
    if (!is.null(marginals) && pathogen_col %in% names(marginals)) {
      marg_k <- marginals[marginals[[pathogen_col]] == path, ]
      marg_k <- marg_k[marg_k[[class_col]] %in% classes, ]

      if (nrow(marg_k) > 0L) {
        bin_mat <- as.matrix(prof_df[, classes, drop = FALSE])
        storage.mode(bin_mat) <- "double"

        for (i in seq_len(nrow(marg_k))) {
          cls <- marg_k[[class_col]][i]
          target <- marg_k[[rate_col]][i]
          if (!cls %in% colnames(bin_mat)) next
          reconstructed <- sum(p_hat * bin_mat[, cls], na.rm = TRUE)
          all_rows <- c(all_rows, list(.row(
            "marginal", paste0("marg_", cls),
            target, reconstructed, "recomputed"
          )))
        }
      }
    }

    # Check 3c: Pairwise constraints recomputed from supplied pairwise
    if (!is.null(pairwise) && pathogen_col %in% names(pairwise)) {
      pair_k <- pairwise[pairwise[[pathogen_col]] == path, ]

      if (nrow(pair_k) > 0L && class1_col %in% names(pair_k) &&
        class2_col %in% names(pair_k) && pairwise_rate_col %in% names(pair_k)) {
        bin_mat <- as.matrix(prof_df[, classes, drop = FALSE])
        storage.mode(bin_mat) <- "double"

        for (i in seq_len(nrow(pair_k))) {
          c1 <- pair_k[[class1_col]][i]
          c2 <- pair_k[[class2_col]][i]
          target <- pair_k[[pairwise_rate_col]][i]
          if (is.na(target) || !c1 %in% colnames(bin_mat) || !c2 %in% colnames(bin_mat)) next
          reconstructed <- sum(p_hat * bin_mat[, c1] * bin_mat[, c2], na.rm = TRUE)
          all_rows <- c(all_rows, list(.row(
            "pairwise", paste0("pair_", c1, "_", c2),
            target, reconstructed, "recomputed"
          )))
        }
      }
    }
  }

  result_tbl <- dplyr::bind_rows(all_rows)

  # Per-pathogen summary printed to console
  if (nrow(result_tbl) > 0L) {
    summary_tbl <- result_tbl %>%
      dplyr::group_by(pathogen) %>%
      dplyr::summarise(
        n_checks = dplyr::n(),
        n_pass = sum(pass, na.rm = TRUE),
        n_fail = sum(!pass, na.rm = TRUE),
        max_abs_residual = max(abs_residual, na.rm = TRUE),
        .groups = "drop"
      )
    message("\n[check_profile_constraints] Results (tolerance = ", tolerance, "):")
    for (i in seq_len(nrow(summary_tbl))) {
      r <- summary_tbl[i, ]
      status <- if (r$n_fail == 0L) "[OK]" else "[!!]"
      message(sprintf(
        "  %s %-40s  %d/%d passed | max |residual| = %.2e",
        status, r$pathogen, r$n_pass, r$n_checks, r$max_abs_residual
      ))
    }
  }

  result_tbl
}


#' Bootstrap Uncertainty Intervals for Resistance Profile Probabilities
#'
#' Quantifies uncertainty in resistance-profile probability estimates by
#' repeatedly resampling isolate counts from a binomial distribution implied
#' by the observed \code{n_tested} and \code{n_resistant} values, refitting
#' the convex optimisation QP for each replicate, and returning percentile
#' confidence intervals across replicates.
#'
#' \strong{Resampling mechanism:}
#' For each bootstrap replicate \eqn{b} and each (pathogen, class) cell:
#' \deqn{n_{\text{resistant}}^{(b)} \sim \text{Binomial}(n_{\text{tested}},\;
#' \hat{r}_{kd})}
#' The bootstrap marginal is then
#' \eqn{\hat{r}_{kd}^{(b)} = n_{\text{resistant}}^{(b)} / n_{\text{tested}}}.
#' Pairwise co-resistance rates (if supplied) are resampled analogously.
#' The QP is re-solved on each bootstrap marginal using the existing
#' \code{compute_resistance_profiles()} engine.
#'
#' \strong{Performance:}
#' For \eqn{n \leq 12} classes the QP solves in milliseconds; B = 500
#' replicates will complete in seconds. For \eqn{n \geq 14} classes use
#' \code{n_cores > 1} or reduce B.
#'
#' @param marginals Data frame. Marginal resistance rates with at minimum
#'   columns \code{pathogen_col}, \code{class_col}, \code{n_tested_col},
#'   \code{n_resistant_col}. This is the \code{$marginal} element from
#'   \code{compute_marginal_resistance()}, or a flat tibble with the same
#'   structure from \code{compute_marginals_from_data()}.
#' @param coresistance_output Named list or \code{NULL}. Output of
#'   \code{compute_pairwise_coresistance()} (one entry per pathogen with
#'   \code{$T_matrix}, \code{$R_matrix}, \code{$prevalence}). When supplied,
#'   pairwise co-resistance is resampled per replicate. When \code{NULL},
#'   independence fallback \eqn{P(A \cap B) = P(A) P(B)} is used in every
#'   replicate. Default \code{NULL}.
#' @param B Integer. Number of bootstrap replicates. Default \code{500}.
#' @param seed Integer. Random seed for reproducibility. Default \code{123}.
#' @param alpha Numeric. Two-sided interval width. Default \code{0.05}
#'   (i.e., 95\% intervals).
#' @param n_cores Integer. Parallel cores via \code{parallel::mclapply}.
#'   Default \code{1L}.
#' @param exclude_near_zero Logical. Passed to \code{compute_resistance_profiles()}.
#'   Default \code{TRUE}.
#' @param top_n_classes Integer or \code{NULL}. Cap on number of classes per
#'   pathogen. Default \code{NULL}.
#' @param sigma_sq Numeric. QP constraint variance. Default \code{1}.
#' @param ridge Numeric. QP Hessian ridge term. Default \code{1e-8}.
#' @param pathogen_col Character. Pathogen column in \code{marginals}.
#'   Default \code{"organism_name"}.
#' @param class_col Character. Class column in \code{marginals}.
#'   Default \code{"antibiotic_class"}.
#' @param n_tested_col Character. Tested-count column. Default \code{"n_tested"}.
#' @param n_resistant_col Character. Resistant-count column. Default
#'   \code{"n_resistant"}.
#' @param org_group_col Character. Organism-group column (passed through to
#'   QP engine). Default \code{"org_group"}.
#'
#' @return A named list, one entry per pathogen. Each entry is a tibble with
#'   columns:
#' \describe{
#'   \item{\code{profile}}{Profile label (e.g. \code{"RSS"}).}
#'   \item{\code{probability_mean}}{Mean profile probability across B replicates.}
#'   \item{\code{probability_median}}{Median across replicates.}
#'   \item{\code{lower}}{Lower percentile (\code{alpha / 2}).}
#'   \item{\code{upper}}{Upper percentile (\code{1 - alpha / 2}).}
#'   \item{\code{n_replicates_converged}}{Number of replicates for which the
#'     QP converged (non-uniform solution).}
#'   \item{\code{convergence_rate}}{Proportion of replicates that converged.}
#' }
#' The point-estimate \code{profiles} tibble (from a single run on the
#' original marginals) is stored as the attribute \code{"point_estimate"}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(amr_clean)
#' co_res <- compute_pairwise_coresistance(marg)
#'
#' boot <- bootstrap_profiles_convex(
#'   marginals = marg$marginal,
#'   coresistance_output = co_res,
#'   B = 500,
#'   seed = 42
#' )
#'
#' # 95% intervals for K. pneumoniae
#' boot[["Klebsiella pneumoniae"]]
#'
#' # Point estimate
#' attr(boot[["Klebsiella pneumoniae"]], "point_estimate")
#' }
bootstrap_profiles_convex <- function(
  marginals,
  coresistance_output = NULL,
  B = 500L,
  seed = 123L,
  alpha = 0.05,
  n_cores = 1L,
  exclude_near_zero = TRUE,
  top_n_classes = NULL,
  sigma_sq = 1,
  ridge = 1e-8,
  pathogen_col = "organism_name",
  class_col = "antibiotic_class",
  n_tested_col = "n_tested",
  n_resistant_col = "n_resistant",
  org_group_col = "org_group"
) {
  # Input validation
  required_marg <- c(pathogen_col, class_col, n_tested_col, n_resistant_col)
  missing_marg <- setdiff(required_marg, names(marginals))
  if (length(missing_marg) > 0L) {
    stop(sprintf(
      "Columns not found in `marginals`: %s",
      paste(missing_marg, collapse = ", ")
    ), call. = FALSE)
  }

  B <- as.integer(B)
  n_cores <- as.integer(n_cores)
  if (is.na(B) || B < 1L) stop("`B` must be a positive integer.", call. = FALSE)
  if (is.na(n_cores) || n_cores < 1L) stop("`n_cores` must be a positive integer.", call. = FALSE)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  # Add org_group if absent (required by compute_resistance_profiles engine)
  if (!org_group_col %in% names(marginals)) {
    marginals[[org_group_col]] <- marginals[[pathogen_col]]
  }

  # Add marginal_resistance if absent
  if (!"marginal_resistance" %in% names(marginals)) {
    marginals$marginal_resistance <- marginals[[n_resistant_col]] / marginals[[n_tested_col]]
  }

  # Point estimate on original marginals (stored as attribute)
  message(sprintf(
    "[bootstrap_profiles_convex] Computing point estimate before bootstrap (B=%d, seed=%d)...",
    B, seed
  ))

  # Build marginal_output list expected by compute_resistance_profiles()
  .build_marginal_output <- function(marg_df, zero_thr = 0) {
    near_zero <- marg_df[marg_df$marginal_resistance <= zero_thr, , drop = FALSE]
    # class_long is not used by the QP step; supply a minimal placeholder
    list(
      marginal   = marg_df,
      near_zero  = near_zero,
      class_long = data.frame()
    )
  }

  point_out <- compute_resistance_profiles(
    marginal_output = .build_marginal_output(marginals),
    coresistance_output = if (is.null(coresistance_output)) list() else coresistance_output,
    exclude_near_zero = exclude_near_zero,
    top_n_classes = top_n_classes,
    sigma_sq = sigma_sq,
    ridge = ridge,
    pathogen_col = pathogen_col,
    antibiotic_class_col = class_col,
    n_cores = 1L
  )

  if (length(point_out) == 0L) {
    stop("Point estimate returned no profiles. Check marginals and class filters.", call. = FALSE)
  }

  # Bootstrap replicates
  set.seed(seed)
  # Pre-draw all random seeds for replicates (ensures reproducibility
  # regardless of parallel vs sequential execution order)
  rep_seeds <- sample.int(.Machine$integer.max, B, replace = FALSE)

  message(sprintf(
    "[bootstrap_profiles_convex] Running %d bootstrap replicates on %d pathogen(s) with %d core(s)...",
    B, length(point_out), n_cores
  ))

  .one_replicate <- function(b) {
    set.seed(rep_seeds[b])

    # Resample resistant counts from Binomial(n_tested, observed_prevalence)
    marg_b <- marginals
    n_t <- marg_b[[n_tested_col]]
    r_obs <- pmin(pmax(marg_b$marginal_resistance, 0), 1) # guard [0,1]
    r_obs[is.na(r_obs)] <- 0

    n_r_b <- stats::rbinom(nrow(marg_b), size = n_t, prob = r_obs)
    marg_b$marginal_resistance <- n_r_b / n_t

    # Resample pairwise T and R matrices if available
    co_b <- coresistance_output
    if (!is.null(co_b) && length(co_b) > 0L) {
      co_b <- lapply(co_b, function(co_path) {
        T_mat <- co_path$T_matrix
        R_mat <- co_path$R_matrix
        prev <- co_path$prevalence
        n_c <- nrow(T_mat)
        R_b <- matrix(0L, n_c, n_c, dimnames = dimnames(T_mat))
        prev_b <- matrix(NA_real_, n_c, n_c, dimnames = dimnames(T_mat))

        if (n_c >= 2L) {
          # Same i-outer/j-inner pair order as the old nested loop, so this
          # keeps rbinom() draw-for-draw reproducible for a given seed.
          idx <- do.call(rbind, lapply(seq_len(n_c - 1L), function(i) cbind(i, (i + 1L):n_c)))
          t_vec <- T_mat[idx]
          valid <- !is.na(t_vec) & t_vec != 0L
          idx <- idx[valid, , drop = FALSE]
          t_vec <- t_vec[valid]

          if (length(t_vec) > 0L) {
            p_vec <- prev[idx]
            p_vec[is.na(p_vec)] <- 0
            p_vec <- pmin(pmax(p_vec, 0), 1)
            r_vec <- stats::rbinom(length(t_vec), size = t_vec, prob = p_vec)

            idx_lo <- idx[, c(2L, 1L), drop = FALSE]
            R_b[idx] <- r_vec
            R_b[idx_lo] <- r_vec
            pv <- r_vec / t_vec
            prev_b[idx] <- pv
            prev_b[idx_lo] <- pv
          }
        }
        diag(prev_b) <- NA_real_
        list(
          prevalence = prev_b,
          T_matrix   = T_mat,
          R_matrix   = R_b,
          classes    = co_path$classes
        )
      })
    }

    tryCatch(
      compute_resistance_profiles(
        marginal_output      = .build_marginal_output(marg_b),
        coresistance_output  = if (is.null(co_b)) list() else co_b,
        exclude_near_zero    = exclude_near_zero,
        top_n_classes        = top_n_classes,
        sigma_sq             = sigma_sq,
        ridge                = ridge,
        pathogen_col         = pathogen_col,
        antibiotic_class_col = class_col,
        n_cores              = 1L
      ),
      error = function(e) NULL # failed replicate returns NULL
    )
  }

  if (n_cores > 1L) {
    all_reps <- parallel::mclapply(seq_len(B), .one_replicate,
      mc.cores = n_cores, mc.preschedule = FALSE
    )
  } else {
    all_reps <- lapply(seq_len(B), .one_replicate)
  }

  # Remove NULL (failed) replicates
  valid_reps <- Filter(Negate(is.null), all_reps)
  n_valid <- length(valid_reps)

  if (n_valid == 0L) stop("All bootstrap replicates failed.", call. = FALSE)
  if (n_valid < B) {
    warning(sprintf("%d/%d bootstrap replicates failed and were dropped.", B - n_valid, B),
      call. = FALSE
    )
  }

  # Aggregate: collect p_hat vectors per pathogen x profile
  lo_q <- alpha / 2
  hi_q <- 1 - alpha / 2

  out <- list()

  for (path in names(point_out)) {
    point_prof <- point_out[[path]]$profiles

    # Collect probability matrix: rows = profiles, cols = replicates
    prob_matrix <- vapply(valid_reps, function(rep_result) {
      if (is.null(rep_result) || !path %in% names(rep_result)) {
        return(rep(NA_real_, nrow(point_prof)))
      }
      rep_profs <- rep_result[[path]]$profiles
      # Align on profile label
      matched <- rep_profs$probability[match(point_prof$profile, rep_profs$profile)]
      ifelse(is.na(matched), 0, matched)
    }, numeric(nrow(point_prof)))

    # Each row of prob_matrix is one profile across B replicates
    # Detect uniform (failed) replicates: all probs equal
    n_prof <- nrow(point_prof)
    uniform_p <- 1 / n_prof
    converged <- apply(prob_matrix, 2, function(col) !all(abs(col - uniform_p) < 1e-10))
    n_converged <- sum(converged, na.rm = TRUE)
    conv_matrix <- prob_matrix[, converged, drop = FALSE]

    if (ncol(conv_matrix) == 0L) {
      warning(sprintf("'%s': no converged bootstrap replicates. Intervals will be NA.", path),
        call. = FALSE
      )
      lower <- rep(NA_real_, n_prof)
      upper <- rep(NA_real_, n_prof)
      pmean <- rep(NA_real_, n_prof)
      pmed <- rep(NA_real_, n_prof)
    } else {
      lower <- apply(conv_matrix, 1, stats::quantile, probs = lo_q, na.rm = TRUE)
      upper <- apply(conv_matrix, 1, stats::quantile, probs = hi_q, na.rm = TRUE)
      pmean <- apply(conv_matrix, 1, mean, na.rm = TRUE)
      pmed <- apply(conv_matrix, 1, stats::median, na.rm = TRUE)
    }

    result_df <- tibble::tibble(
      profile                  = point_prof$profile,
      probability_mean         = round(pmean, 6L),
      probability_median       = round(pmed, 6L),
      lower                    = round(lower, 6L),
      upper                    = round(upper, 6L),
      n_replicates_converged   = n_converged,
      convergence_rate         = round(n_converged / n_valid, 4L)
    )

    attr(result_df, "point_estimate") <- point_prof
    out[[path]] <- result_df

    message(sprintf(
      "[bootstrap] '%s': %d profiles | %d/%d replicates converged | interval width: %.1f%%",
      path, n_prof, n_converged, n_valid, 100 * (1 - alpha)
    ))
  }

  message(sprintf(
    "\n[bootstrap_profiles_convex] Done. %d pathogen(s) | %d/%d valid replicates.",
    length(out), n_valid, B
  ))

  out
}


# enumerate_binary_profiles()

#' Enumerate All Binary Resistance Profiles for a Set of Antibiotic Classes
#'
#' Generates all \eqn{2^n} binary resistance profiles for \code{n} antibiotic
#' classes using vectorised bit manipulation. Each row is one profile; each
#' class column contains 1 (resistant) or 0 (susceptible).
#'
#' @param classes Character vector. Ordered antibiotic class names. The order
#'   determines which bit position each class occupies and must be consistent
#'   across all calls within a single analysis.
#'
#' @return Tibble with \code{2^n} rows and \code{n + 1} columns:
#'   \code{profile_delta} (character label, e.g. \code{"RSS"}) and one integer
#'   column per class (0 / 1).
#'
#' @export
enumerate_binary_profiles <- function(classes) {
  if (length(classes) == 0L) stop("`classes` must be a non-empty character vector.", call. = FALSE)
  if (length(classes) > 20L) {
    warning(sprintf("2^%d = %d profiles may be very slow. Consider top_n_classes.", length(classes), 2L^length(classes)), call. = FALSE)
  }

  n <- length(classes)
  n_profiles <- 2L^n

  profiles_mat <- matrix(
    as.integer(outer(0L:(n_profiles - 1L), 2L^(0L:(n - 1L)), bitwAnd) > 0L),
    nrow = n_profiles, ncol = n,
    dimnames = list(NULL, classes)
  )

  char_mat <- matrix("S", nrow = n_profiles, ncol = n)
  char_mat[profiles_mat == 1L] <- "R"

  # check.names = FALSE preserves class names that contain hyphens, spaces, or
  # other characters that base R would otherwise sanitise to dots.
  label_df <- as.data.frame(char_mat, check.names = FALSE)
  colnames(label_df) <- classes
  binary_df <- as.data.frame(profiles_mat, check.names = FALSE)

  tibble::as_tibble(
    cbind(
      data.frame(
        profile_delta = do.call(paste0, label_df),
        stringsAsFactors = FALSE, check.names = FALSE
      ),
      binary_df
    ),
    .name_repair = "minimal"
  )
}


# build_constraint_matrix()

#' Build QP Constraint Matrix and Target Vector
#'
#' Constructs the constraint matrix \strong{M} and target vector \strong{v}
#' used in the simplex-constrained weighted least-squares QP:
#'
#' \deqn{\min_p \|Mp - v\|^2_W \quad \text{s.t.} \quad p \ge 0,\; \sum p = 1}
#'
#' Marginal rows: \eqn{M_{d,\delta} = 1} iff class \eqn{d} is resistant in
#' profile \eqn{\delta}.
#' Pairwise rows: \eqn{M_{d_1 d_2,\delta} = 1} iff both classes are resistant.
#' When a pairwise value is unavailable (too few co-tested isolates), the
#' product of marginals (independence assumption) is substituted.
#' Pairwise values that exceed \eqn{\min(P(A), P(B))} are capped.
#'
#' @param profiles_enum Tibble from \code{enumerate_binary_profiles()}. The
#'   \code{profile_delta} column is ignored; only the binary class columns are
#'   used.
#' @param r_marg Named numeric vector. Marginal resistance rates, one per
#'   class. Names must match column names in \code{profiles_enum}.
#' @param co_mat Numeric matrix or \code{NULL}. Square symmetric pairwise
#'   co-resistance prevalence matrix with row/column names matching
#'   \code{names(r_marg)}. \code{NA} cells trigger independence fallback.
#'   Default \code{NULL}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{M}}{Numeric matrix \eqn{(n + n(n-1)/2) \times 2^n}.}
#'   \item{\code{v}}{Numeric vector of constraint targets.}
#'   \item{\code{constraint_names}}{Character vector labelling each row of M.}
#'   \item{\code{fallback_pairs}}{Character vector of class pairs that used
#'     the independence fallback.}
#'   \item{\code{capped_pairs}}{Named numeric vector of pairs whose pairwise
#'     value was capped, showing original and capped values.}
#' }
#'
#' @export
build_constraint_matrix <- function(profiles_enum, r_marg, co_mat = NULL) {
  classes <- names(r_marg)
  n <- length(classes)

  # Binary indicator matrix (2^n x n)
  profiles_mat <- as.matrix(profiles_enum[, classes, drop = FALSE])
  storage.mode(profiles_mat) <- "double"

  # Marginal rows
  M_marg <- t(profiles_mat) # n x 2^n
  v_marg <- r_marg[classes]
  marg_names <- paste0("marg_", classes)

  # Pairwise rows
  pairs_mat <- utils::combn(n, 2L)
  n_pair <- ncol(pairs_mat)
  d1_idx <- pairs_mat[1L, ]
  d2_idx <- pairs_mat[2L, ]
  c1_names <- classes[d1_idx]
  c2_names <- classes[d2_idx]
  pair_names <- paste0("pair_", c1_names, "_", c2_names)

  M_pair <- t(
    profiles_mat[, d1_idx, drop = FALSE] *
      profiles_mat[, d2_idx, drop = FALSE]
  ) # n(n-1)/2 x 2^n

  r1 <- r_marg[c1_names]
  r2 <- r_marg[c2_names]

  # Look up observed pairwise values
  co_vals <- rep(NA_real_, n_pair)
  if (!is.null(co_mat) && !is.null(rownames(co_mat))) {
    in_r <- c1_names %in% rownames(co_mat)
    in_c <- c2_names %in% colnames(co_mat)
    ok <- in_r & in_c
    if (any(ok)) co_vals[ok] <- co_mat[cbind(c1_names[ok], c2_names[ok])]
  }

  # Cap pairwise to min(P(A), P(B))
  cap_vals <- pmin(r1, r2)
  has_co <- !is.na(co_vals)
  capped_co <- pmin(co_vals, cap_vals)
  was_capped <- has_co & !is.na(capped_co) & (capped_co < co_vals)
  v_pair <- ifelse(has_co, capped_co, r1 * r2)

  fallback_pairs <- pair_names[!has_co]
  capped_pairs <- if (any(was_capped)) {
    stats::setNames(
      paste0(round(co_vals[was_capped], 4L), "->", round(capped_co[was_capped], 4L)),
      pair_names[was_capped]
    )
  } else {
    character(0L)
  }

  M <- rbind(M_marg, M_pair)
  v <- c(v_marg, v_pair)
  storage.mode(M) <- "double"

  list(
    M                = M,
    v                = v,
    constraint_names = c(marg_names, pair_names),
    fallback_pairs   = fallback_pairs,
    capped_pairs     = capped_pairs
  )
}


# validate_aggregate_inputs()

#' Validate Pre-computed Aggregate Marginal Inputs
#'
#' Validates a pre-computed marginal resistance table (e.g. from GBD or a
#' surveillance network) before it is passed to \code{estimate_profiles_convex()}.
#' Checks prevalence bounds, duplicate rows, zero-denominator cells, impossible
#' pairwise values, and minimum class counts per pathogen.
#'
#' @param marginals Data frame with at minimum columns \code{pathogen_col},
#'   \code{class_col}, and \code{rate_col}. Optionally also \code{n_tested_col}
#'   and \code{n_resistant_col}.
#' @param pairwise Data frame or \code{NULL}. Pairwise co-resistance table with
#'   columns \code{pathogen_col}, \code{class1_col}, \code{class2_col}, and
#'   \code{pairwise_rate_col}. Default \code{NULL}.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param rate_col Character. Marginal resistance rate column.
#'   Default \code{"marginal_resistance"}.
#' @param n_tested_col Character or \code{NULL}. Default \code{"n_tested"}.
#' @param n_resistant_col Character or \code{NULL}. Default \code{"n_resistant"}.
#' @param class1_col Character. Default \code{"antibiotic_class_1"}.
#' @param class2_col Character. Default \code{"antibiotic_class_2"}.
#' @param pairwise_rate_col Character. Default \code{"pairwise_resistance_prevalence"}.
#' @param min_classes_per_pathogen Integer. Pathogens with fewer classes are
#'   flagged. Default \code{2L}.
#'
#' @return Tibble of check results (same schema as
#'   \code{validate_profile_inputs()}). Stops on any \code{"fail"}.
#'
#' @export
validate_aggregate_inputs <- function(
  marginals,
  pairwise = NULL,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rate_col = "marginal_resistance",
  n_tested_col = "n_tested",
  n_resistant_col = "n_resistant",
  class1_col = "antibiotic_class_1",
  class2_col = "antibiotic_class_2",
  pairwise_rate_col = "pairwise_resistance_prevalence",
  min_classes_per_pathogen = 2L
) {
  results <- list()

  # 1. Required columns
  required <- c(pathogen_col, class_col, rate_col)
  missing <- setdiff(required, names(marginals))
  if (length(missing) > 0L) {
    results <- .add_check(
      results, "required_cols", "fail",
      sprintf("Required column(s) missing from marginals: %s", paste(missing, collapse = ", "))
    )
  } else {
    results <- .add_check(
      results, "required_cols", "pass",
      sprintf("All required columns present (%s).", paste(required, collapse = ", "))
    )
  }

  if ("fail" %in% sapply(results, `[[`, "status")) {
    stop(dplyr::bind_rows(results)$message[dplyr::bind_rows(results)$status == "fail"][[1]],
      call. = FALSE
    )
  }

  # 2. Prevalence bounds [0, 1]
  rates <- marginals[[rate_col]]
  n_oob <- sum(!is.na(rates) & (rates < 0 | rates > 1))
  if (n_oob > 0L) {
    results <- .add_check(
      results, "prevalence_bounds", "fail",
      sprintf("%d row(s) have %s outside [0, 1].", n_oob, rate_col), n_oob
    )
  } else {
    results <- .add_check(
      results, "prevalence_bounds", "pass",
      sprintf("All %s values in [0, 1].", rate_col)
    )
  }

  # 3. n_tested > 0 (if column present)
  if (!is.null(n_tested_col) && n_tested_col %in% names(marginals)) {
    n_zero <- sum(!is.na(marginals[[n_tested_col]]) & marginals[[n_tested_col]] <= 0L)
    if (n_zero > 0L) {
      results <- .add_check(
        results, "n_tested_positive", "fail",
        sprintf("%d row(s) have %s <= 0.", n_zero, n_tested_col), n_zero
      )
    } else {
      results <- .add_check(
        results, "n_tested_positive", "pass",
        sprintf("All %s > 0.", n_tested_col)
      )
    }

    # n_resistant <= n_tested
    if (!is.null(n_resistant_col) && n_resistant_col %in% names(marginals)) {
      n_bad <- sum(!is.na(marginals[[n_resistant_col]]) &
        !is.na(marginals[[n_tested_col]]) &
        marginals[[n_resistant_col]] > marginals[[n_tested_col]])
      if (n_bad > 0L) {
        results <- .add_check(
          results, "resistant_le_tested", "fail",
          sprintf("%d row(s) have %s > %s.", n_bad, n_resistant_col, n_tested_col), n_bad
        )
      } else {
        results <- .add_check(
          results, "resistant_le_tested", "pass",
          sprintf("%s <= %s for all rows.", n_resistant_col, n_tested_col)
        )
      }
    }
  }

  # 4. No duplicate pathogen x class rows
  dup_key <- paste(marginals[[pathogen_col]], marginals[[class_col]], sep = "||")
  n_dup <- sum(duplicated(dup_key))
  if (n_dup > 0L) {
    results <- .add_check(
      results, "no_duplicates", "fail",
      sprintf("%d duplicate pathogen x class row(s) found.", n_dup), n_dup
    )
  } else {
    results <- .add_check(
      results, "no_duplicates", "pass",
      "No duplicate pathogen x class rows."
    )
  }

  # 5. Minimum classes per pathogen
  class_counts <- tapply(marginals[[class_col]], marginals[[pathogen_col]], length)
  n_too_few <- sum(class_counts < min_classes_per_pathogen)
  if (n_too_few > 0L) {
    pathogens_few <- names(class_counts)[class_counts < min_classes_per_pathogen]
    results <- .add_check(
      results, "min_classes", "warn",
      sprintf(
        "%d pathogen(s) have < %d classes (profiles trivial): %s",
        n_too_few, min_classes_per_pathogen,
        paste(head(pathogens_few, 5L), collapse = ", ")
      ),
      n_too_few
    )
  } else {
    results <- .add_check(
      results, "min_classes", "pass",
      sprintf("All pathogens have >= %d classes.", min_classes_per_pathogen)
    )
  }

  # 6. Pairwise checks (if supplied)
  if (!is.null(pairwise)) {
    pw_req <- c(pathogen_col, class1_col, class2_col, pairwise_rate_col)
    pw_miss <- setdiff(pw_req, names(pairwise))
    if (length(pw_miss) > 0L) {
      results <- .add_check(
        results, "pairwise_cols", "fail",
        sprintf("Pairwise table missing column(s): %s", paste(pw_miss, collapse = ", "))
      )
    } else {
      # Pairwise rate in [0,1]
      pw_rates <- pairwise[[pairwise_rate_col]]
      n_pw_oob <- sum(!is.na(pw_rates) & (pw_rates < 0 | pw_rates > 1))
      if (n_pw_oob > 0L) {
        results <- .add_check(
          results, "pairwise_bounds", "fail",
          sprintf("%d pairwise row(s) have %s outside [0, 1].", n_pw_oob, pairwise_rate_col), n_pw_oob
        )
      }

      # Pairwise <= min(marginal_A, marginal_B)
      pw_joined <- pairwise %>%
        dplyr::left_join(marginals[, c(pathogen_col, class_col, rate_col)],
          by = stats::setNames(c(pathogen_col, class_col), c(pathogen_col, class1_col))
        ) %>%
        dplyr::rename(rate_A = !!rate_col) %>%
        dplyr::left_join(marginals[, c(pathogen_col, class_col, rate_col)],
          by = stats::setNames(c(pathogen_col, class_col), c(pathogen_col, class2_col))
        ) %>%
        dplyr::rename(rate_B = !!rate_col)

      n_impossible <- sum(!is.na(pw_joined[[pairwise_rate_col]]) &
        !is.na(pw_joined$rate_A) & !is.na(pw_joined$rate_B) &
        pw_joined[[pairwise_rate_col]] > pmin(pw_joined$rate_A, pw_joined$rate_B))
      if (n_impossible > 0L) {
        results <- .add_check(
          results, "pairwise_le_marginals", "warn",
          sprintf("%d pairwise value(s) exceed min(marginal_A, marginal_B) and will be capped.", n_impossible), n_impossible
        )
      } else {
        results <- .add_check(
          results, "pairwise_le_marginals", "pass",
          "All pairwise values <= min(marginal_A, marginal_B)."
        )
      }
    }
  }

  result_tbl <- dplyr::bind_rows(results)
  n_fail <- sum(result_tbl$status == "fail")
  n_warn <- sum(result_tbl$status == "warn")
  n_pass <- sum(result_tbl$status == "pass")

  message(sprintf(
    "[validate_aggregate_inputs] %d passed | %d warnings | %d failed",
    n_pass, n_warn, n_fail
  ))
  if (n_warn > 0L) {
    for (msg in result_tbl$message[result_tbl$status == "warn"]) message("  [!] ", msg)
  }

  if (n_fail > 0L) {
    stop(sprintf(
      "%d check(s) failed:\n%s", n_fail,
      paste(sprintf("  - %s", result_tbl$message[result_tbl$status == "fail"]),
        collapse = "\n"
      )
    ), call. = FALSE)
  }

  invisible(result_tbl)
}


#' Compute Marginal Resistance Rates from Preprocessed Wide Data
#'
#' Takes the wide-format tibble produced by \code{preprocess_for_profiles()}
#' and computes marginal resistance rates per pathogen x antibiotic class,
#' optionally stratified by geography, year, and/or patient outcome. The
#' counting unit is always \code{isolate_id}: a patient with multiple isolates
#' contributes one count per isolate.
#'
#' Optionally, externally modelled marginals (e.g. from GBD ST-GPR) can be
#' supplied to override the locally computed values while retaining the
#' pairwise correlation structure from the local data.
#'
#' @param data_wide Tibble. Output of \code{preprocess_for_profiles()$data_wide}.
#'   One row per isolate; antibiotic-class columns contain \code{"S"},
#'   \code{"R"}, or \code{NA} (not tested).
#' @param col_map Named list. \code{col_map_resolved} from
#'   \code{preprocess_for_profiles()}.
#' @param panel_map Named list. Same panel used in \code{preprocess_for_profiles()}.
#' @param stratify_by Character vector or \code{NULL}. Dimensions to stratify
#'   over. Valid: \code{"geography"}, \code{"year"}. Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Column name for patient outcome.
#'   When supplied, separate marginals are computed for each outcome value
#'   (e.g. Died vs Discharged). Default \code{NULL}.
#' @param min_n_tested Integer. Minimum tested-isolate count per
#'   (stratum x pathogen x class) cell. Cells below this threshold are dropped
#'   with a logged reason. Default \code{30L}.
#' @param external_marginals Data frame or \code{NULL}. Pre-modelled marginal
#'   resistance rates (e.g. GBD ST-GPR estimates). When provided, locally
#'   computed rates are replaced where a match exists on
#'   pathogen x class [x geography x year]. Must have columns matching
#'   \code{ext_col_map}. Default \code{NULL}.
#' @param ext_col_map Named list. Column names in \code{external_marginals}.
#'   Default \code{list(pathogen_col="pathogen", class_col="antibiotic_class",
#'   geography_col="geography", year_col="year", rate_col="resistance_prevalence")}.
#'
#' @return Tibble with columns: \code{pathogen}, \code{antibiotic_class},
#'   [stratum columns], \code{n_tested}, \code{n_resistant},
#'   \code{marginal_resistance}, \code{marginal_source}
#'   (\code{"computed"} or \code{"external"}).
#'
#' @export
compute_marginals_from_data <- function(
  data_wide,
  col_map,
  panel_map,
  stratify_by = NULL,
  outcome_col = NULL,
  min_n_tested = 30L,
  external_marginals = NULL,
  ext_col_map = list(
    pathogen_col  = "pathogen",
    class_col     = "antibiotic_class",
    geography_col = "geography",
    year_col      = "year",
    rate_col      = "resistance_prevalence"
  )
) {
  cols <- .resolve_class_cols(data_wide, col_map, panel_map, outcome_col)
  pathogen_col <- cols$pathogen_col
  geography_col <- cols$geography_col
  isolate_col <- cols$isolate_col
  class_cols <- cols$class_cols

  if (length(class_cols) == 0L) {
    stop("No antibiotic-class columns found in data_wide matching panel_map. Check panel_map class names.", call. = FALSE)
  }

  strat_cols <- .build_strat_cols(stratify_by, outcome_col, geography_col, data_wide)

  group_vars <- c(strat_cols, pathogen_col)

  # Pivot class columns to long, then compute marginals
  data_long <- data_wide %>%
    dplyr::select(dplyr::all_of(c(isolate_col, group_vars, class_cols))) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(class_cols),
      names_to  = "antibiotic_class",
      values_to = "class_result"
    ) %>%
    dplyr::filter(!is.na(class_result)) # only tested isolate-class pairs

  marginals_tbl <- data_long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "antibiotic_class")))) %>%
    dplyr::summarise(
      n_tested    = dplyr::n(),
      n_resistant = sum(class_result == "R", na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      marginal_resistance = n_resistant / n_tested,
      marginal_source     = "computed"
    )

  # Apply min_n_tested filter
  n_before <- nrow(marginals_tbl)
  excluded_min <- marginals_tbl[marginals_tbl$n_tested < min_n_tested, , drop = FALSE]
  marginals_tbl <- marginals_tbl[marginals_tbl$n_tested >= min_n_tested, , drop = FALSE]

  if (nrow(excluded_min) > 0L) {
    message(sprintf(
      "[compute_marginals_from_data] min_n_tested = %d: %d / %d pathogen-class cells excluded.",
      min_n_tested, nrow(excluded_min), n_before
    ))
  }

  # External marginals override
  if (!is.null(external_marginals)) {
    message("[compute_marginals_from_data] Merging external marginals override...")

    ext_path_col <- .null_default(ext_col_map$pathogen_col, "pathogen")
    ext_cls_col <- .null_default(ext_col_map$class_col, "antibiotic_class")
    ext_rate_col <- .null_default(ext_col_map$rate_col, "resistance_prevalence")
    ext_geo_col <- .null_default(ext_col_map$geography_col, "geography")
    ext_yr_col <- .null_default(ext_col_map$year_col, "year")

    join_keys <- c(pathogen_col = ext_path_col, antibiotic_class = ext_cls_col)
    if ("geography" %in% stratify_by && ext_geo_col %in% names(external_marginals)) {
      join_keys <- c(join_keys, stats::setNames(ext_geo_col, geography_col))
    }
    if ("year" %in% stratify_by && ext_yr_col %in% names(external_marginals)) {
      join_keys <- c(join_keys, stats::setNames(ext_yr_col, "year"))
    }

    ext_slim <- external_marginals[, c(unname(join_keys), ext_rate_col), drop = FALSE]
    names(ext_slim)[names(ext_slim) == ext_rate_col] <- ".ext_rate"

    marginals_tbl <- marginals_tbl %>%
      dplyr::left_join(ext_slim, by = join_keys) %>%
      dplyr::mutate(
        marginal_source     = dplyr::if_else(!is.na(.ext_rate), "external", marginal_source),
        marginal_resistance = dplyr::coalesce(.ext_rate, marginal_resistance)
      ) %>%
      dplyr::select(-.ext_rate)

    n_overridden <- sum(marginals_tbl$marginal_source == "external")
    message(sprintf("  %d / %d cells overridden with external marginals.", n_overridden, nrow(marginals_tbl)))
  }

  message(sprintf(
    "[compute_marginals_from_data] Done. %d pathogen(s), %d cells, %d strat col(s).",
    dplyr::n_distinct(marginals_tbl[[pathogen_col]]), nrow(marginals_tbl), length(strat_cols)
  ))

  tibble::as_tibble(marginals_tbl)
}


#' Compute Pairwise Co-resistance Using Pearson Back-calculation
#'
#' For each antibiotic-class pair (A, B) and each pathogen [x stratum], this
#' function:
#' \enumerate{
#'   \item Computes the Pearson correlation \eqn{\rho_{AB}} between binary
#'     resistance indicators across co-tested isolates.
#'   \item Back-calculates pairwise prevalence using the GBD formula:
#'     \deqn{P(A \cap B) = P(A) P(B) + \rho_{AB}
#'       \sqrt{P(A)(1-P(A)) \cdot P(B)(1-P(B))}}
#'     where \eqn{P(A)} and \eqn{P(B)} come from \code{marginals}
#'     (which may include externally overridden values).
#'   \item Caps the result to \eqn{\min(P(A), P(B))} and floors it to 0.
#'     Cells below the floor use the independence product as fallback.
#' }
#'
#' @param data_wide Tibble. Output of \code{preprocess_for_profiles()$data_wide}.
#' @param marginals Tibble. Output of \code{compute_marginals_from_data()}.
#' @param col_map Named list. \code{col_map_resolved}.
#' @param panel_map Named list. Same panel as used in preprocessing.
#' @param stratify_by Character vector or \code{NULL}. Must match what was
#'   used in \code{compute_marginals_from_data()}. Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Default \code{NULL}.
#' @param min_co_tested Integer. Pairs with fewer co-tested isolates are
#'   reported with \code{method = "independence_fallback"}. Default \code{10L}.
#'
#' @return Tibble with one row per pathogen x class-pair [x stratum]:
#' \describe{
#'   \item{\code{pathogen}, \code{antibiotic_class_1}, \code{antibiotic_class_2}}{Keys.}
#'   \item{\code{n_co_tested}}{Isolates tested for both classes.}
#'   \item{\code{rho}}{Pearson correlation between binary resistance indicators.}
#'   \item{\code{pairwise_prevalence}}{Back-calculated P(AnB), capped and floored.}
#'   \item{\code{method}}{\code{"pearson_back_calc"} or \code{"independence_fallback"}.}
#' }
#'
#' @export
compute_pairwise_from_data <- function(
  data_wide,
  marginals,
  col_map,
  panel_map,
  stratify_by = NULL,
  outcome_col = NULL,
  min_co_tested = 10L
) {
  cols <- .resolve_class_cols(data_wide, col_map, panel_map, outcome_col)
  pathogen_col <- cols$pathogen_col
  geography_col <- cols$geography_col
  isolate_col <- cols$isolate_col
  class_cols <- cols$class_cols

  strat_cols <- .build_strat_cols(stratify_by, outcome_col, geography_col, data_wide)

  group_key_cols <- c(strat_cols, pathogen_col)

  # Get unique strata x pathogen combinations
  strata_df <- data_wide %>%
    dplyr::select(dplyr::all_of(c(group_key_cols))) %>%
    dplyr::distinct()

  all_pairs <- list()

  for (i in seq_len(nrow(strata_df))) {
    row_key <- strata_df[i, , drop = FALSE]
    path <- row_key[[pathogen_col]]
    cls_panel <- panel_map[[path]]
    if (is.null(cls_panel)) next
    cls_avail <- intersect(cls_panel, class_cols)
    if (length(cls_avail) < 2L) next

    # Filter data_wide to this stratum
    stratum_filter <- rep(TRUE, nrow(data_wide))
    for (sc in strat_cols) {
      stratum_filter <- stratum_filter & (data_wide[[sc]] == row_key[[sc]])
    }
    stratum_filter <- stratum_filter & (data_wide[[pathogen_col]] == path)

    sub_wide <- data_wide[stratum_filter, cls_avail, drop = FALSE]

    # Get marginals for this stratum x pathogen
    marg_filter <- marginals[[pathogen_col]] == path
    for (sc in strat_cols) {
      marg_filter <- marg_filter & (marginals[[sc]] == row_key[[sc]])
    }
    marg_k <- marginals[marg_filter, , drop = FALSE]
    p_lookup <- stats::setNames(marg_k$marginal_resistance, marg_k$antibiotic_class)

    # Enumerate class pairs
    pairs_mat <- utils::combn(cls_avail, 2L)
    for (j in seq_len(ncol(pairs_mat))) {
      c1 <- pairs_mat[1L, j]
      c2 <- pairs_mat[2L, j]
      if (!c1 %in% names(p_lookup) || !c2 %in% names(p_lookup)) next

      PA <- p_lookup[[c1]]
      PB <- p_lookup[[c2]]
      if (is.na(PA) || is.na(PB)) next

      # Co-tested isolates: non-NA in both columns
      v1 <- sub_wide[[c1]]
      v2 <- sub_wide[[c2]]
      both <- !is.na(v1) & !is.na(v2)
      n_co <- sum(both)

      if (n_co < min_co_tested) {
        pairwise_val <- PA * PB
        method <- "independence_fallback"
        rho <- NA_real_
      } else {
        y1 <- as.integer(v1[both] == "R")
        y2 <- as.integer(v2[both] == "R")
        rho <- tryCatch(
          stats::cor(y1, y2),
          error = function(e) NA_real_,
          warning = function(w) suppressWarnings(stats::cor(y1, y2))
        )

        if (is.na(rho) || !is.finite(rho)) {
          pairwise_val <- PA * PB
          method <- "independence_fallback"
        } else {
          # GBD back-calculation formula
          var_A <- PA * (1 - PA)
          var_B <- PB * (1 - PB)
          p_ab_raw <- PA * PB + rho * sqrt(var_A * var_B)
          # Cap to min(PA, PB) and floor to independence product
          p_ab_capped <- pmin(p_ab_raw, min(PA, PB))
          pairwise_val <- max(p_ab_capped, PA * PB * 0) # floor at 0, not independence
          pairwise_val <- max(pairwise_val, 0)
          method <- "pearson_back_calc"
        }
      }

      pair_row <- c(
        as.list(row_key),
        list(
          antibiotic_class_1   = c1,
          antibiotic_class_2   = c2,
          n_co_tested          = n_co,
          rho                  = round(rho, 6L),
          pairwise_prevalence  = round(pairwise_val, 6L),
          method               = method
        )
      )
      all_pairs <- c(all_pairs, list(tibble::as_tibble(pair_row)))
    }
  }

  if (length(all_pairs) == 0L) {
    warning("No pairwise pairs computed. Check panel_map and data.", call. = FALSE)
    return(tibble::tibble())
  }

  result <- dplyr::bind_rows(all_pairs)
  n_pearson <- sum(result$method == "pearson_back_calc")
  n_fallback <- sum(result$method == "independence_fallback")
  message(sprintf(
    "[compute_pairwise_from_data] %d pairs: %d Pearson back-calc | %d independence fallback.",
    nrow(result), n_pearson, n_fallback
  ))
  result
}


#' Estimate Resistance Profile Probabilities via Convex Optimisation
#'
#' The public-facing solver for Pathway 1. Accepts pre-computed marginal and
#' pairwise resistance rates (from \code{compute_marginals_from_data()} /
#' \code{compute_pairwise_from_data()}, or from \code{validate_aggregate_inputs()}
#' for GBD-style pre-computed inputs), enumerates all \eqn{2^n} binary
#' resistance profiles, builds the constraint matrix, and solves the
#' simplex-constrained weighted least-squares QP:
#'
#' \deqn{\hat{p} = \arg\min_{p \in \Delta} \|Mp - v\|^2 + \lambda \|p\|^2}
#'
#' Stratification columns (geography, year, outcome) present in \code{marginals}
#' are detected automatically and carried through to the output.
#'
#' @param marginals Tibble. Output of \code{compute_marginals_from_data()} or
#'   a validated aggregate table. Must have at minimum \code{pathogen_col},
#'   \code{class_col}, and \code{rate_col}.
#' @param pairwise Tibble or \code{NULL}. Output of
#'   \code{compute_pairwise_from_data()}. When \code{NULL}, independence
#'   \eqn{P(A \cap B) = P(A)P(B)} is assumed for all class pairs.
#'   Default \code{NULL}.
#' @param panel_map Named list or \code{NULL}. If supplied, classes for each
#'   pathogen are restricted to the panel. When \code{NULL}, all classes
#'   present in \code{marginals} are used. Default \code{NULL}.
#' @param exclude_near_zero Logical. Drop classes with
#'   \code{marginal_resistance <= zero_threshold}. Default \code{TRUE}.
#' @param zero_threshold Numeric. Classes at or below this value are
#'   considered near-zero. Default \code{0}.
#' @param top_n_classes Integer or \code{NULL}. Cap classes per pathogen to
#'   top N by \code{n_tested}. Limits 2^n explosion. Default \code{NULL}.
#' @param lambda Numeric. Ridge regularisation added to QP Hessian diagonal.
#'   Default \code{1e-8}.
#' @param sigma_sq Numeric. Uniform constraint variance. Default \code{1}.
#' @param solver Character. \code{"osqp"} (default, preferred) or
#'   \code{"quadprog"} (fallback).
#' @param n_cores Integer. Parallel cores. Default \code{1L}.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param rate_col Character. Default \code{"marginal_resistance"}.
#' @param n_tested_col Character. Default \code{"n_tested"}.
#'
#' @return Tibble with one row per (stratum x) pathogen x profile:
#' \describe{
#'   \item{Stratum columns}{Any geography / year / outcome columns from \code{marginals}.}
#'   \item{\code{pathogen}}{Pathogen name.}
#'   \item{\code{profile_set_type}}{\code{"aggregate_convex"}.}
#'   \item{\code{profile_class_set}}{Ordered class names joined by \code{"|"}.}
#'   \item{\code{profile_delta}}{Binary profile label, e.g. \code{"RSS"}.}
#'   \item{\code{profile_probability}}{Estimated probability \eqn{\hat{p}_\delta}.}
#'   \item{\code{estimator}}{\code{"convex"}.}
#'   \item{\code{convergence_flag}}{Logical; \code{FALSE} if QP failed and
#'     uniform distribution was returned.}
#'   \item{\code{identifiability_flag}}{Logical; \code{TRUE} when the constraint
#'     matrix is rank-deficient (system is underdetermined).}
#'   \item{\code{max_abs_residual}}{Maximum absolute constraint residual.}
#'   \item{\code{notes}}{Semicolon-separated notes: capped pairs, fallback pairs.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginals_from_data(result$data_wide, result$col_map_resolved, panel_map)
#' pw <- compute_pairwise_from_data(result$data_wide, marg, result$col_map_resolved, panel_map)
#' out <- estimate_profiles_convex(marg, pw, panel_map)
#' }
estimate_profiles_convex <- function(
  marginals,
  pairwise = NULL,
  panel_map = NULL,
  exclude_near_zero = TRUE,
  zero_threshold = 0,
  top_n_classes = NULL,
  lambda = 1e-8,
  sigma_sq = 1,
  solver = c("osqp", "quadprog"),
  n_cores = 1L,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rate_col = "marginal_resistance",
  n_tested_col = "n_tested"
) {
  solver <- match.arg(solver)
  n_cores <- as.integer(n_cores)

  has_osqp <- requireNamespace("osqp", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)
  has_quadprog <- requireNamespace("quadprog", quietly = TRUE)
  if (!has_osqp && !has_quadprog) {
    stop("Install 'osqp' + 'Matrix' (recommended) or 'quadprog'.", call. = FALSE)
  }
  use_osqp <- has_osqp && solver != "quadprog"

  req_cols <- c(pathogen_col, class_col, rate_col)
  miss <- setdiff(req_cols, names(marginals))
  if (length(miss) > 0L) {
    stop(sprintf("Column(s) not found in `marginals`: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }

  # Detect stratum columns automatically
  known_non_strat <- c(
    pathogen_col, class_col, rate_col, n_tested_col,
    "n_resistant", "marginal_source", "org_group"
  )
  strat_cols <- setdiff(names(marginals), known_non_strat)

  # Build unique (stratum x pathogen) combinations
  key_cols <- c(strat_cols, pathogen_col)
  strata_df <- marginals %>%
    dplyr::select(dplyr::all_of(key_cols)) %>%
    dplyr::distinct()

  message(sprintf(
    "[estimate_profiles_convex] %d stratum x pathogen combination(s) | solver: %s | n_cores: %d",
    nrow(strata_df), if (use_osqp) "osqp" else "quadprog", n_cores
  ))

  # Build a look-up structure for pairwise prevalence
  .get_co_mat <- function(path, row_key) {
    if (is.null(pairwise) || nrow(pairwise) == 0L) {
      return(NULL)
    }
    pw_f <- pairwise[[pathogen_col]] == path
    for (sc in strat_cols) {
      pw_f <- pw_f & (pairwise[[sc]] == row_key[[sc]])
    }
    pw_k <- pairwise[pw_f, , drop = FALSE]
    if (nrow(pw_k) == 0L) {
      return(NULL)
    }

    all_cls <- unique(c(pw_k$antibiotic_class_1, pw_k$antibiotic_class_2))
    mat <- matrix(NA_real_, length(all_cls), length(all_cls),
      dimnames = list(all_cls, all_cls)
    )
    for (ri in seq_len(nrow(pw_k))) {
      c1 <- pw_k$antibiotic_class_1[ri]
      c2 <- pw_k$antibiotic_class_2[ri]
      pv <- pw_k$pairwise_prevalence[ri]
      mat[c1, c2] <- pv
      mat[c2, c1] <- pv
    }
    mat
  }

  .solve_one <- function(idx) {
    row_key <- strata_df[idx, , drop = FALSE]
    path <- row_key[[pathogen_col]]

    # Filter marginals for this stratum x pathogen
    marg_f <- marginals[[pathogen_col]] == path
    for (sc in strat_cols) marg_f <- marg_f & (marginals[[sc]] == row_key[[sc]])
    marg_k <- marginals[marg_f, , drop = FALSE]

    # Determine classes
    if (!is.null(panel_map) && path %in% names(panel_map)) {
      cls <- sort(intersect(panel_map[[path]], marg_k[[class_col]]))
    } else {
      cls <- sort(marg_k[[class_col]])
    }

    if (exclude_near_zero) {
      near_z <- marg_k[[class_col]][!is.na(marg_k[[rate_col]]) &
        marg_k[[rate_col]] <= zero_threshold]
      cls <- setdiff(cls, near_z)
    }

    if (!is.null(top_n_classes) && length(cls) > top_n_classes &&
      n_tested_col %in% names(marg_k)) {
      n_ord <- marg_k[marg_k[[class_col]] %in% cls, ]
      n_ord <- n_ord[order(n_ord[[n_tested_col]], decreasing = TRUE), ]
      cls <- sort(n_ord[[class_col]][seq_len(top_n_classes)])
    }

    if (length(cls) < 2L) {
      message(sprintf("  '%s': fewer than 2 classes after filtering -- skipped.", path))
      return(NULL)
    }

    r_marg <- stats::setNames(marg_k[[rate_col]][match(cls, marg_k[[class_col]])], cls)
    co_mat_k <- .get_co_mat(path, row_key)

    # Enumerate profiles and build constraint matrix
    profiles_enum <- enumerate_binary_profiles(cls)
    cm <- build_constraint_matrix(profiles_enum, r_marg, co_mat_k)
    M <- cm$M
    v <- cm$v
    n_profiles <- nrow(profiles_enum)

    # Identifiability: rank(M) vs 2^n
    rank_M <- qr(M)$rank
    identif_flag <- rank_M < n_profiles

    # QP
    coef <- 2.0 / sigma_sq
    H_mat <- coef * crossprod(M)
    diag(H_mat) <- diag(H_mat) + lambda
    d_qp <- coef * drop(crossprod(M, v))

    converged <- TRUE
    p_hat <- tryCatch(
      {
        if (use_osqp) {
          A_sp <- rbind(
            Matrix::Matrix(1.0, nrow = 1L, ncol = n_profiles, sparse = TRUE),
            Matrix::Diagonal(n_profiles)
          )
          prob <- osqp::osqp(
            P = Matrix::forceSymmetric(Matrix::Matrix(H_mat)),
            q = -d_qp,
            A = A_sp,
            l = c(1.0, rep(0.0, n_profiles)),
            u = c(1.0, rep(Inf, n_profiles)),
            pars = osqp::osqpSettings(
              verbose = FALSE, eps_abs = 1e-8,
              eps_rel = 1e-8, max_iter = 10000L, polish = TRUE
            )
          )
          res <- prob$solve()
          if (!(res$info$status %in% c("solved", "solved_inaccurate"))) stop(res$info$status)
          pmax(res$x, 0.0)
        } else {
          Amat <- cbind(rep(1.0, n_profiles), diag(n_profiles))
          bvec <- c(1.0, rep(0.0, n_profiles))
          pmax(quadprog::solve.QP(H_mat, d_qp, Amat, bvec, meq = 1L)$solution, 0.0)
        }
      },
      error = function(e) {
        converged <<- FALSE
        rep(1.0 / n_profiles, n_profiles)
      }
    )
    p_hat <- p_hat / sum(p_hat)

    residuals <- drop(M %*% p_hat) - v
    max_abs_res <- max(abs(residuals))
    notes_parts <- character(0)
    if (length(cm$capped_pairs) > 0L) {
      notes_parts <- c(notes_parts, paste0("capped: ", paste(names(cm$capped_pairs), collapse = ", ")))
    }
    if (length(cm$fallback_pairs) > 0L) {
      notes_parts <- c(notes_parts, paste0("indep_fallback: ", paste(cm$fallback_pairs, collapse = ", ")))
    }
    if (!converged) notes_parts <- c(notes_parts, "QP_failed_uniform_returned")

    message(sprintf(
      "  [%d/%d] '%s'%s: n=%d classes -> %d profiles | max|resid|=%.2e%s%s",
      idx, nrow(strata_df), path,
      if (length(strat_cols) > 0L) paste0(" [", paste(unlist(row_key[strat_cols]), collapse = "/"), "]") else "",
      length(cls), n_profiles, max_abs_res,
      if (!converged) " [QP FAILED]" else "",
      if (identif_flag) " [UNDERDETERMINED]" else ""
    ))

    out_df <- profiles_enum %>%
      dplyr::mutate(
        profile_probability  = round(p_hat, 8L),
        convergence_flag     = converged,
        identifiability_flag = identif_flag,
        max_abs_residual     = round(max_abs_res, 8L),
        notes                = if (length(notes_parts) > 0L) paste(notes_parts, collapse = "; ") else NA_character_
      ) %>%
      dplyr::mutate(
        pathogen = path,
        profile_set_type = "aggregate_convex",
        profile_class_set = paste(cls, collapse = "|"),
        estimator = "convex"
      )

    # Attach stratum columns
    for (sc in strat_cols) out_df[[sc]] <- row_key[[sc]]

    # Reorder columns: strat | pathogen | profile meta | probability | flags
    front_cols <- c(
      strat_cols, "pathogen", "profile_set_type", "profile_class_set",
      "profile_delta", "profile_probability", "estimator",
      "convergence_flag", "identifiability_flag", "max_abs_residual", "notes"
    )
    class_indicator_cols <- setdiff(names(out_df), front_cols)
    out_df <- out_df[, c(front_cols, class_indicator_cols), drop = FALSE]
    tibble::as_tibble(out_df)
  }

  if (n_cores > 1L) {
    results_list <- parallel::mclapply(seq_len(nrow(strata_df)), .solve_one,
      mc.cores = n_cores, mc.preschedule = FALSE
    )
  } else {
    results_list <- lapply(seq_len(nrow(strata_df)), .solve_one)
  }

  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0L) {
    stop("No profiles estimated. Check marginals, panel_map, and class filters.", call. = FALSE)
  }

  final_tbl <- dplyr::bind_rows(results_list)
  message(sprintf(
    "\n[estimate_profiles_convex] Done. %d row(s) | %d pathogen(s) | %d stratum combination(s).",
    nrow(final_tbl),
    dplyr::n_distinct(final_tbl$pathogen),
    if (length(strat_cols) > 0L) dplyr::n_distinct(final_tbl[, strat_cols, drop = FALSE]) else 1L
  ))
  final_tbl
}


# ===========================================================================
# Pathway 1 -- Isolate-level engine
#
# Computes resistance profiles directly from line-level isolate data in a
# three-step pipeline. This engine serves the DALY YLL and YLD calculations
# and accepts both pre-processed facility data and survey datasets.
#
#   Step 1  compute_marginal_resistance()
#     Collapses drug-level records to antibiotic-class level per isolate
#     (class R if the isolate is resistant to ANY drug in that class), then
#     computes marginal resistance rates per pathogen x class. Classes at or
#     below zero_threshold are flagged in $near_zero for downstream review.
#
#   Step 2  compute_pairwise_coresistance()
#     Builds pairwise co-resistance matrices from the collapsed class-level
#     data produced in Step 1:
#       T[k, i, j]    number of isolates tested for both class i and class j
#       R[k, i, j]    number resistant to both
#       Prev[k, i, j] = R / T  (NA where T < min_co_tested)
#
#   Step 3  compute_resistance_profiles()
#     Enumerates all 2^n binary resistance profiles and recovers profile
#     probabilities via a simplex-constrained weighted least-squares QP
#     (GBD methodology, eq. 7.5.1.3). Marginal and pairwise co-resistance
#     rates serve as constraints; the solver falls back to a uniform
#     distribution if the QP fails.
#
# Unit of analysis: isolate (one unique isolate_col value per pathogen).
# Class-level resistance rule: R_{e,k,c} = 1 if resistant to ANY drug in c.
#
# Reference: Antimicrobial Resistance Collaborators. Lancet. 2022.
# ===========================================================================


# -- Internal helper -----------------------------------------------------------

#' Validate that required columns exist in a data frame
#' @keywords internal
.check_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }
}


# -- Step 1 --------------------------------------------------------------------

#' Compute Marginal Resistance per Pathogen and Antibiotic Class
#'
#' Collapses drug-level susceptibility data to antibiotic-class level per
#' isolate (\eqn{R_{e,k,c} = 1} if resistant to \strong{any} drug in class
#' \eqn{c}), then computes marginal resistance for every pathogen x class
#' combination found in the data.
#'
#' Classes whose marginal resistance is at or below \code{zero_threshold} are
#' listed in \code{$near_zero} as a flag for downstream use -- they are
#' \strong{not} removed here.
#'
#' @param data Data frame. Pre-processed AMR data at isolate x antibiotic
#'   level (one row per isolate-antibiotic combination). Results must already
#'   be binary (\code{"S"} / \code{"R"}); no reclassification is applied.
#' @param pathogen_col Character. Column with pathogen names.
#'   Default \code{"organism_name"}.
#' @param org_group_col Character. Column with organism group labels.
#'   Default \code{"org_group"}.
#' @param isolate_col Character. Column uniquely identifying each isolate.
#'   Default \code{"isolate_id"}.
#' @param antibiotic_class_col Character. Column with the antibiotic class
#'   for each drug. Default \code{"antibiotic_class"}.
#' @param antibiotic_value_col Character. Column with susceptibility result
#'   (\code{"S"} or \code{"R"}). Default \code{"antibiotic_value"}.
#' @param zero_threshold Numeric. Classes with
#'   \code{marginal_resistance <= zero_threshold} are listed in
#'   \code{$near_zero}. Default \code{0}.
#' @param min_n_tested Integer or \code{NULL}. Minimum number of isolates that
#'   must have been tested for a pathogen-class combination to be retained.
#'   Combinations with \code{n_tested < min_n_tested} are dropped from
#'   \code{$marginal} \strong{and} from \code{$class_long}, so the exclusion
#'   propagates automatically into \code{compute_pairwise_coresistance()} and
#'   \code{compute_resistance_profiles()}. Set to \code{NULL} or \code{0} to
#'   disable the filter. Default \code{30}.
#' @param facility_col Character or \code{NULL}. Name of the column identifying
#'   the facility/site. When provided together with \code{facility_name}, data
#'   are filtered to the specified facility \strong{before} any computation.
#'   Both \code{facility_col} and \code{facility_name} must be supplied together
#'   or both left \code{NULL}. When provided, \code{facility_col} is also
#'   retained in \code{$class_long} and \code{$marginal} so that downstream
#'   steps can apply the same filter. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. The facility value to retain
#'   (matched via \code{==} against \code{facility_col}). Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Name of the column containing
#'   patient outcomes (e.g. \code{"final_outcome"}). When provided together with
#'   \code{outcome_value}, data are filtered to isolates with the specified
#'   outcome \strong{before} any computation. Both \code{outcome_col} and
#'   \code{outcome_value} must be supplied together or both left \code{NULL}.
#'   When provided, \code{outcome_col} is retained in \code{$class_long} and
#'   \code{$marginal} so that downstream steps can apply the same filter.
#'   Default \code{NULL}.
#' @param outcome_value Character or \code{NULL}. The outcome value to retain
#'   (e.g. \code{"discharged"}, \code{"dead"}; matched via \code{==} against
#'   \code{outcome_col}). Default \code{NULL}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{marginal}}{Data frame with columns: \code{pathogen_col},
#'     \code{org_group_col}, \code{antibiotic_class_col}, \code{n_tested},
#'     \code{n_resistant}, \code{marginal_resistance}. Sorted descending by
#'     \code{marginal_resistance} within each pathogen.}
#'   \item{\code{near_zero}}{Subset of \code{marginal} where
#'     \code{marginal_resistance <= zero_threshold}. These classes are
#'     candidates for exclusion in downstream profiling.}
#'   \item{\code{class_long}}{Collapsed isolate x pathogen x class data frame
#'     (columns: \code{isolate_col}, \code{pathogen_col}, \code{org_group_col},
#'     \code{antibiotic_class_col}, \code{class_result}). Pass this directly
#'     to \code{compute_pairwise_coresistance()}.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(
#'   data                 = amr_clean,
#'   pathogen_col         = "organism_name",
#'   org_group_col        = "org_group",
#'   isolate_col          = "isolate_id",
#'   antibiotic_class_col = "antibiotic_class",
#'   antibiotic_value_col = "antibiotic_value"
#' )
#'
#' marg$marginal # full marginal resistance table
#' marg$near_zero # classes flagged as near-zero
#' }
compute_marginal_resistance <- function(
  data,
  pathogen_col = "organism_name",
  org_group_col = "org_group",
  isolate_col = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_value_col = "antibiotic_value",
  zero_threshold = 0,
  min_n_tested = 30,
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL
) {
  .check_cols(data, c(
    pathogen_col, org_group_col, isolate_col,
    antibiotic_class_col, antibiotic_value_col
  ))

  # -- Optional facility filter ----------------------------------------------

  .check_paired_args(facility_col, facility_name, "facility_col", "facility_name")
  if (!is.null(facility_col)) {
    .check_cols(data, facility_col)
    n_before <- nrow(data)
    data <- data[data[[facility_col]] == facility_name, , drop = FALSE]
    if (nrow(data) == 0L) {
      stop(sprintf("No rows found for %s = '%s'.", facility_col, facility_name))
    }
    message(sprintf(
      "Facility filter applied: %s = '%s' (%d -> %d rows).",
      facility_col, facility_name, n_before, nrow(data)
    ))
  }

  # -- Optional outcome filter -----------------------------------------------

  .check_paired_args(outcome_col, outcome_value, "outcome_col", "outcome_value")
  if (!is.null(outcome_col)) {
    .check_cols(data, outcome_col)
    n_before <- nrow(data)
    data <- data[data[[outcome_col]] == outcome_value, , drop = FALSE]
    if (nrow(data) == 0L) {
      stop(sprintf("No rows found for %s = '%s'.", outcome_col, outcome_value))
    }
    message(sprintf(
      "Outcome filter applied: %s = '%s' (%d -> %d rows).",
      outcome_col, outcome_value, n_before, nrow(data)
    ))
  }

  # -- Collapse to isolate x pathogen x class -------------------------------
  #
  # R_{e,k,c} = 1 if resistant to ANY drug in class c.

  iso_grp <- c(isolate_col, pathogen_col, org_group_col, antibiotic_class_col)
  if (!is.null(facility_col)) iso_grp <- c(facility_col, iso_grp)
  if (!is.null(outcome_col)) iso_grp <- c(outcome_col, iso_grp)

  class_long <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(iso_grp))) %>%
    dplyr::summarise(
      class_result = dplyr::if_else(
        any(!!rlang::sym(antibiotic_value_col) == "R"), "R", "S"
      ),
      .groups = "drop"
    )

  message(sprintf(
    "Collapsed to class level: %d isolates | %d pathogens | %d classes.",
    dplyr::n_distinct(class_long[[isolate_col]]),
    dplyr::n_distinct(class_long[[pathogen_col]]),
    dplyr::n_distinct(class_long[[antibiotic_class_col]])
  ))

  # -- Marginal resistance per pathogen x class -----------------------------

  marg_grp <- c(pathogen_col, org_group_col, antibiotic_class_col)
  if (!is.null(facility_col)) marg_grp <- c(facility_col, marg_grp)
  if (!is.null(outcome_col)) marg_grp <- c(outcome_col, marg_grp)

  marginal <- class_long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(marg_grp))) %>%
    dplyr::summarise(
      n_tested    = dplyr::n(),
      n_resistant = sum(class_result == "R", na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      marginal_resistance = n_resistant / n_tested
    ) %>%
    dplyr::arrange(
      !!rlang::sym(pathogen_col),
      dplyr::desc(marginal_resistance)
    )

  # -- Flag near-zero classes ------------------------------------------------

  near_zero <- marginal %>%
    dplyr::filter(marginal_resistance <= zero_threshold)

  if (nrow(near_zero) > 0) {
    message(sprintf(
      "%d pathogen-class combination(s) have marginal_resistance <= %g (listed in $near_zero):",
      nrow(near_zero), zero_threshold
    ))
    print(
      near_zero[, c(
        pathogen_col, antibiotic_class_col,
        "n_tested", "n_resistant", "marginal_resistance"
      )],
      row.names = FALSE
    )
  } else {
    message(sprintf(
      "No classes with marginal_resistance <= %g.", zero_threshold
    ))
  }

  # -- Minimum-tests threshold (optional) ------------------------------------
  # Drop pathogen-class combinations with too few isolates tested.
  # class_long is filtered to match so the exclusion carries through to
  # compute_pairwise_coresistance() and compute_resistance_profiles().

  if (!is.null(min_n_tested) && min_n_tested > 0) {
    n_before <- nrow(marginal)
    excluded <- marginal[marginal$n_tested < min_n_tested, , drop = FALSE]
    marginal <- marginal[marginal$n_tested >= min_n_tested, , drop = FALSE]

    if (nrow(excluded) > 0L) {
      message(sprintf(
        "min_n_tested = %d: %d of %d pathogen-class combination(s) excluded (n_tested < %d):",
        min_n_tested, nrow(excluded), n_before, min_n_tested
      ))
      print(
        excluded[, c(pathogen_col, antibiotic_class_col, "n_tested"),
          drop = FALSE
        ],
        row.names = FALSE
      )
    } else {
      message(sprintf(
        "min_n_tested = %d: all %d pathogen-class combination(s) passed the threshold.",
        min_n_tested, n_before
      ))
    }

    # Re-derive near_zero from the already-filtered marginal
    near_zero <- marginal[marginal$marginal_resistance <= zero_threshold, ,
      drop = FALSE
    ]

    # Align class_long: remove isolate rows for excluded pathogen-class combos.
    # Join keys: pathogen x class (facility/outcome are already constant if filtered).
    join_keys <- c(pathogen_col, antibiotic_class_col)
    if (!is.null(facility_col)) join_keys <- c(facility_col, join_keys)
    if (!is.null(outcome_col)) join_keys <- c(outcome_col, join_keys)
    class_long <- dplyr::semi_join(
      class_long,
      marginal[, join_keys, drop = FALSE],
      by = join_keys
    )

    message(sprintf(
      "class_long after min_n_tested filter: %d isolate-class row(s) retained.",
      nrow(class_long)
    ))
  }

  return(list(
    marginal   = marginal,
    near_zero  = near_zero,
    class_long = class_long
  ))
}


# -- Step 2 --------------------------------------------------------------------

#' Compute Pairwise Co-resistance Matrices per Pathogen
#'
#' For every pathogen and every pair of antibiotic classes \eqn{(c_i, c_j)}
#' that were both tested:
#' \deqn{
#'   T_{k,i,j} = \sum_e \mathbf{1}(c_i \text{ tested} \land c_j \text{ tested})
#' }
#' \deqn{
#'   R_{k,i,j} = \sum_e \mathbf{1}(R_{e,k,c_i}=1 \land R_{e,k,c_j}=1)
#' }
#' \deqn{
#'   \text{Prev}_{k,i,j} = R_{k,i,j} \;/\; T_{k,i,j}
#' }
#'
#' Matrices are computed for \strong{all} tested classes -- no filtering by
#' marginal resistance or GBD core list is applied here.
#'
#' @param marginal_output The list returned by
#'   \code{compute_marginal_resistance()}. The \code{$class_long} element is
#'   used as input.
#' @param pathogen_col Character. Must match the column name used in Step 1.
#'   Default \code{"organism_name"}.
#' @param org_group_col Character. Default \code{"org_group"}.
#' @param isolate_col Character. Default \code{"isolate_id"}.
#' @param antibiotic_class_col Character. Default \code{"antibiotic_class"}.
#' @param min_co_tested Integer. Pairwise cells with fewer than
#'   \code{min_co_tested} co-tested isolates are set to \code{NA} in the
#'   prevalence matrix. Default \code{10}.
#' @param facility_col Character or \code{NULL}. Column identifying the
#'   facility/site. When provided together with \code{facility_name}, filters
#'   \code{class_long} to the specified facility before building matrices.
#'   For this to work, \code{compute_marginal_resistance()} must have been
#'   called with the same \code{facility_col} argument so that the column is
#'   present in \code{$class_long}. Both must be supplied together or both
#'   left \code{NULL}. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. Facility value to retain.
#'   Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Column containing patient
#'   outcomes. When provided together with \code{outcome_value}, filters
#'   \code{class_long} to isolates with the specified outcome before building
#'   matrices. \code{compute_marginal_resistance()} must have been called with
#'   the same \code{outcome_col} so that the column is present in
#'   \code{$class_long}. Both must be supplied together or both left
#'   \code{NULL}. Default \code{NULL}.
#' @param outcome_value Character or \code{NULL}. Outcome value to retain
#'   (e.g. \code{"discharged"}, \code{"dead"}). Default \code{NULL}.
#'
#' @return Named list. One entry per pathogen (keyed by pathogen name), each
#'   containing:
#' \describe{
#'   \item{\code{prevalence}}{Square symmetric matrix of pairwise co-resistance
#'     rates. Diagonal and cells with \code{n < min_co_tested} are \code{NA}.}
#'   \item{\code{T_matrix}}{Integer matrix of co-tested isolate counts.}
#'   \item{\code{R_matrix}}{Integer matrix of co-resistant isolate counts.}
#'   \item{\code{classes}}{Character vector of class names (row/column order).}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(amr_clean, ...)
#' co_res <- compute_pairwise_coresistance(marg)
#'
#' # Prevalence matrix for K. pneumoniae
#' co_res[["Klebsiella pneumoniae"]]$prevalence
#'
#' # Co-tested counts
#' co_res[["Klebsiella pneumoniae"]]$T_matrix
#' }
compute_pairwise_coresistance <- function(
  marginal_output,
  pathogen_col = "organism_name",
  org_group_col = "org_group",
  isolate_col = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  min_co_tested = 10,
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL
) {
  if (!is.list(marginal_output) || !"class_long" %in% names(marginal_output)) {
    stop("marginal_output must be the list returned by compute_marginal_resistance().")
  }

  class_long <- marginal_output$class_long

  .check_cols(class_long, c(
    isolate_col, pathogen_col,
    org_group_col, antibiotic_class_col
  ))

  # -- Optional facility filter ----------------------------------------------

  .check_paired_args(facility_col, facility_name, "facility_col", "facility_name")
  if (!is.null(facility_col)) {
    if (!facility_col %in% names(class_long)) {
      stop(sprintf(
        paste0(
          "facility_col '%s' not found in class_long. ",
          "Re-run compute_marginal_resistance() with the same ",
          "facility_col and facility_name to filter upstream."
        ),
        facility_col
      ))
    }
    n_before <- nrow(class_long)
    class_long <- class_long[class_long[[facility_col]] == facility_name, , drop = FALSE]
    if (nrow(class_long) == 0L) {
      stop(sprintf("No rows in class_long for %s = '%s'.", facility_col, facility_name))
    }
    message(sprintf(
      "Facility filter applied: %s = '%s' (%d -> %d rows in class_long).",
      facility_col, facility_name, n_before, nrow(class_long)
    ))
  }

  # -- Optional outcome filter -----------------------------------------------

  .check_paired_args(outcome_col, outcome_value, "outcome_col", "outcome_value")
  if (!is.null(outcome_col)) {
    if (!outcome_col %in% names(class_long)) {
      stop(sprintf(
        paste0(
          "outcome_col '%s' not found in class_long. ",
          "Re-run compute_marginal_resistance() with the same ",
          "outcome_col and outcome_value to filter upstream."
        ),
        outcome_col
      ))
    }
    n_before <- nrow(class_long)
    class_long <- class_long[class_long[[outcome_col]] == outcome_value, , drop = FALSE]
    if (nrow(class_long) == 0L) {
      stop(sprintf("No rows in class_long for %s = '%s'.", outcome_col, outcome_value))
    }
    message(sprintf(
      "Outcome filter applied: %s = '%s' (%d -> %d rows in class_long).",
      outcome_col, outcome_value, n_before, nrow(class_long)
    ))
  }

  pathogens <- sort(unique(class_long[[pathogen_col]]))
  out <- list()

  for (path in pathogens) {
    org_data <- class_long[class_long[[pathogen_col]] == path, ]
    classes <- sort(unique(org_data[[antibiotic_class_col]]))
    n_c <- length(classes)

    if (n_c < 2) {
      message(sprintf("'%s': fewer than 2 classes tested, skipping.", path))
      next
    }

    # Wide: rows = isolates, cols = classes, values = "R" / "S" / NA (not tested)
    iso_wide <- org_data %>%
      tidyr::pivot_wider(
        id_cols     = !!rlang::sym(isolate_col),
        names_from  = !!rlang::sym(antibiotic_class_col),
        values_from = class_result
        # isolates not tested for a class -> NA
      )

    # Binary matrix: 1 = R, 0 = S, NA = not tested
    bin_mat <- iso_wide[, classes, drop = FALSE] %>%
      dplyr::mutate(dplyr::across(
        dplyr::everything(),
        ~ dplyr::case_when(.x == "R" ~ 1L, .x == "S" ~ 0L, TRUE ~ NA_integer_)
      )) %>%
      as.matrix()

    # Pairwise T and R
    pairwise_T <- matrix(0L, n_c, n_c, dimnames = list(classes, classes))
    pairwise_R <- matrix(0L, n_c, n_c, dimnames = list(classes, classes))

    for (i in seq_len(n_c)) {
      for (j in i:n_c) {
        both_tested <- !is.na(bin_mat[, i]) & !is.na(bin_mat[, j])
        pairwise_T[i, j] <- pairwise_T[j, i] <- sum(both_tested)
        pairwise_R[i, j] <- pairwise_R[j, i] <- sum(
          both_tested & bin_mat[, i] == 1L & bin_mat[, j] == 1L
        )
      }
    }

    # Prevalence matrix: NA where co-tested < min_co_tested or on diagonal
    prev_mat <- ifelse(
      pairwise_T < min_co_tested,
      NA_real_,
      pairwise_R / pairwise_T
    )
    dimnames(prev_mat) <- list(classes, classes)
    diag(prev_mat) <- NA_real_

    out[[path]] <- list(
      prevalence = prev_mat,
      T_matrix   = pairwise_T,
      R_matrix   = pairwise_R,
      classes    = classes
    )

    message(sprintf(
      "'%s': %d classes, %d isolates -- co-resistance matrix built.",
      path, n_c, nrow(bin_mat)
    ))
  }

  return(out)
}


# -- Step 3 --------------------------------------------------------------------

#' Compute Resistance Profile Probabilities per Pathogen
#'
#' For each pathogen \eqn{k} with \eqn{n_k} antibiotic classes, enumerates all
#' \eqn{2^{n_k}} binary resistance profiles \eqn{\delta \in \{0,1\}^{n_k}} and
#' estimates their probabilities by solving a simplex-constrained weighted
#' least-squares Quadratic Programme (GBD equation 7.5.1.3):
#'
#' \deqn{
#'   \hat{p} = \arg\min_{p \in \Delta_{2^n}}
#'             \sum_{i=1}^{m} \frac{(m_i^\top p - v_i)^2}{\sigma_i^2}
#' }
#'
#' where \eqn{m = n(n+1)/2} data-derived linear constraints encode
#' \strong{n marginal} resistance rates and \strong{n(n-1)/2 pairwise}
#' co-resistance rates, and \eqn{\Delta} is the standard probability simplex.
#'
#' \subsection{Constraint rows in M}{
#'   \describe{
#'     \item{Marginal (rows 1 to n)}{
#'       Row \eqn{d}: \eqn{M_{d,\delta} = 1} iff \eqn{\delta_d = 1}
#'       (i.e., class \eqn{d} is resistant in profile \eqn{\delta}).
#'       Constraint: \eqn{\sum_\delta M_{d,\delta}\,p_\delta = \hat{r}_{kd}}.
#'     }
#'     \item{Pairwise (rows n+1 to m)}{
#'       Row \eqn{(d_1,d_2)}: \eqn{M_{d_1 d_2,\delta} = 1} iff
#'       \eqn{\delta_{d_1} = 1 \land \delta_{d_2} = 1}.
#'       Constraint: \eqn{\sum_\delta M_{d_1 d_2,\delta}\,p_\delta =
#'       \hat{r}_{k,d_1 d_2}}.
#'       When a pairwise estimate is unavailable (too few co-tested isolates),
#'       the product of marginals (independence assumption) is used as fallback.
#'     }
#'   }
#' }
#'
#' The QP is solved via \code{quadprog::solve.QP}. A small ridge term
#' (\code{ridge}) is added to the Hessian to guarantee strict
#' positive-definiteness. On solver failure the pathogen gets a uniform
#' distribution over all profiles.
#'
#' \subsection{Performance Notes}{
#'   This function has been optimized for speed with vectorized profile generation
#'   and label creation (10-100x faster than previous versions). However,
#'   computational complexity is still exponential in the number of classes:
#'   \itemize{
#'     \item \strong{n <= 14}: Fast (seconds to minutes)
#'     \item \strong{n = 15}: Moderate (minutes)
#'     \item \strong{n >= 16}: Slow and memory-intensive (use \code{top_n_classes})
#'   }
#'   For large datasets with many pathogens, consider using \code{top_n_classes}
#'   to limit each pathogen to its most-tested classes (e.g., \code{top_n_classes = 12}).
#' }
#'
#' @param marginal_output    List returned by \code{compute_marginal_resistance()}.
#' @param coresistance_output List returned by
#'   \code{compute_pairwise_coresistance()}.
#' @param pathogens          Character vector. Pathogen(s) to process.
#'   \code{NULL} (default) processes every pathogen in \code{marginal_output}.
#' @param top_n_pathogens    Integer or \code{NULL} (default). When set, only
#'   the top \code{top_n_pathogens} pathogens ranked by total isolates tested
#'   (sum of \code{n_tested} across all antibiotic classes, descending) are
#'   processed. Applied after the \code{pathogens} argument filter. Useful for
#'   focusing on the most data-rich pathogens -- e.g. \code{top_n_pathogens = 5}
#'   runs profiles for only the 5 most-tested pathogens.
#' @param exclude_near_zero  Logical. If \code{TRUE} (default), antibiotic
#'   classes that appear in \code{marginal_output$near_zero} for a given
#'   pathogen are excluded from profile enumeration.
#' @param top_n_classes      Integer or \code{NULL} (default). When set, only
#'   the top \code{top_n_classes} antibiotic classes ranked by \code{n_tested}
#'   (descending) are kept per pathogen before profile enumeration. Useful for
#'   capping the combinatorial explosion (2^n profiles) for pathogens tested
#'   against many drug classes -- e.g. \code{top_n_classes = 5} gives at most
#'   32 profiles. Applied after \code{exclude_near_zero}.
#' @param sigma_sq           Positive numeric. Assumed variance for each
#'   constraint (uniform). Default \code{1}.
#' @param ridge              Positive numeric. Ridge term added to the QP
#'   Hessian for numerical stability. Default \code{1e-8}.
#' @param pathogen_col       Character. Column name for pathogens. Must match
#'   the column used in Steps 1-2. Default \code{"organism_name"}.
#' @param antibiotic_class_col Character. Column name for antibiotic classes.
#'   Default \code{"antibiotic_class"}.
#' @param facility_col Character or \code{NULL}. Column identifying the
#'   facility/site. When provided together with \code{facility_name}, filters
#'   \code{marginal_output$marginal} and \code{marginal_output$near_zero} to
#'   the specified facility before profile enumeration. For this to work,
#'   \code{compute_marginal_resistance()} must have been called with the same
#'   \code{facility_col} argument. Both must be supplied together or both left
#'   \code{NULL}. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. Facility value to retain.
#'   Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Column containing patient
#'   outcomes. When provided together with \code{outcome_value}, filters
#'   \code{marginal_output$marginal} and \code{marginal_output$near_zero} to
#'   the specified outcome before profile enumeration.
#'   \code{compute_marginal_resistance()} must have been called with the same
#'   \code{outcome_col} argument. Both must be supplied together or both left
#'   \code{NULL}. Default \code{NULL}.
#' @param outcome_value Character or \code{NULL}. Outcome value to retain
#'   (e.g. \code{"discharged"}, \code{"dead"}). Default \code{NULL}.
#' @param n_cores Integer. Number of CPU cores for parallel computation.
#'   Default \code{1L} (sequential).
#'
#' @return Named list, one entry per pathogen, each a list with:
#' \describe{
#'   \item{\code{profiles}}{Data frame: \code{profile} (character label,
#'     e.g.\ \code{"RSR"}), \code{probability} (\eqn{\hat{p}_\delta}), and
#'     one binary (0/1) integer column per antibiotic class indicating whether
#'     that class is resistant (\code{1}) or susceptible (\code{0}) in each
#'     profile.}
#'   \item{\code{classes}}{Character vector of antibiotic class names used
#'     (alphabetical; bit 0 = classes[1]).}
#'   \item{\code{n_classes}}{Integer. Number of classes used.}
#'   \item{\code{constraint_residuals}}{Named numeric vector of
#'     \eqn{m_i^\top \hat{p} - v_i} for each constraint. Small absolute
#'     values indicate good constraint satisfaction.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(amr_clean)
#' co_res <- compute_pairwise_coresistance(marg)
#' rp <- compute_resistance_profiles(marg, co_res)
#'
#' # All profiles and probabilities for K. pneumoniae
#' rp[["Klebsiella pneumoniae"]]$profiles
#'
#' # Constraint residuals (quality check)
#' rp[["Klebsiella pneumoniae"]]$constraint_residuals
#'
#' # Single pathogen
#' rp_kp <- compute_resistance_profiles(
#'   marg, co_res,
#'   pathogens = "Klebsiella pneumoniae"
#' )
#' }
compute_resistance_profiles <- function(
  marginal_output,
  coresistance_output,
  pathogens = NULL,
  top_n_pathogens = NULL,
  exclude_near_zero = TRUE,
  top_n_classes = NULL,
  sigma_sq = 1,
  ridge = 1e-8,
  pathogen_col = "organism_name",
  antibiotic_class_col = "antibiotic_class",
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL,
  n_cores = 1L
) {
  # -- Input validation -------------------------------------------------------
  if (!is.list(marginal_output) ||
    !all(c("marginal", "near_zero", "class_long") %in% names(marginal_output))) {
    stop("marginal_output must be the list returned by compute_marginal_resistance().")
  }
  if (!is.list(coresistance_output)) {
    stop("coresistance_output must be the list returned by compute_pairwise_coresistance().")
  }
  has_osqp <- requireNamespace("osqp", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)
  has_quadprog <- requireNamespace("quadprog", quietly = TRUE)
  if (!has_osqp && !has_quadprog) {
    stop("Either 'osqp' + 'Matrix' (recommended) or 'quadprog' must be installed.")
  }

  if (!is.null(n_cores)) {
    n_cores <- as.integer(n_cores)
    if (is.na(n_cores) || n_cores < 1L) stop("n_cores must be a positive integer.")
  } else {
    n_cores <- 1L
  }

  if (!is.null(top_n_pathogens)) {
    if (!is.numeric(top_n_pathogens) || length(top_n_pathogens) != 1 ||
      top_n_pathogens < 1 || top_n_pathogens != round(top_n_pathogens)) {
      stop("top_n_pathogens must be a single positive integer.")
    }
    top_n_pathogens <- as.integer(top_n_pathogens)
  }

  if (!is.null(top_n_classes)) {
    if (!is.numeric(top_n_classes) || length(top_n_classes) != 1 ||
      top_n_classes < 1 || top_n_classes != round(top_n_classes)) {
      stop("top_n_classes must be a single positive integer.")
    }
    top_n_classes <- as.integer(top_n_classes)
  }

  # -- Optional facility filter ----------------------------------------------

  .check_paired_args(facility_col, facility_name, "facility_col", "facility_name")
  if (!is.null(facility_col)) {
    if (!facility_col %in% names(marginal_output$marginal)) {
      stop(sprintf(
        paste0(
          "facility_col '%s' not found in marginal_output$marginal. ",
          "Re-run compute_marginal_resistance() with the same ",
          "facility_col and facility_name to filter upstream."
        ),
        facility_col
      ))
    }
    marginal_output$marginal <- marginal_output$marginal[
      marginal_output$marginal[[facility_col]] == facility_name, ,
      drop = FALSE
    ]
    marginal_output$near_zero <- marginal_output$near_zero[
      marginal_output$near_zero[[facility_col]] == facility_name, ,
      drop = FALSE
    ]
    if (nrow(marginal_output$marginal) == 0L) {
      stop(sprintf("No rows in marginal for %s = '%s'.", facility_col, facility_name))
    }
    message(sprintf(
      "Facility filter applied: %s = '%s'.", facility_col, facility_name
    ))
  }

  # -- Optional outcome filter -----------------------------------------------

  .check_paired_args(outcome_col, outcome_value, "outcome_col", "outcome_value")
  if (!is.null(outcome_col)) {
    if (!outcome_col %in% names(marginal_output$marginal)) {
      stop(sprintf(
        paste0(
          "outcome_col '%s' not found in marginal_output$marginal. ",
          "Re-run compute_marginal_resistance() with the same ",
          "outcome_col and outcome_value to filter upstream."
        ),
        outcome_col
      ))
    }
    marginal_output$marginal <- marginal_output$marginal[
      marginal_output$marginal[[outcome_col]] == outcome_value, ,
      drop = FALSE
    ]
    marginal_output$near_zero <- marginal_output$near_zero[
      marginal_output$near_zero[[outcome_col]] == outcome_value, ,
      drop = FALSE
    ]
    if (nrow(marginal_output$marginal) == 0L) {
      stop(sprintf("No rows in marginal for %s = '%s'.", outcome_col, outcome_value))
    }
    message(sprintf(
      "Outcome filter applied: %s = '%s'.", outcome_col, outcome_value
    ))
  }

  all_pathogens <- sort(unique(marginal_output$marginal[[pathogen_col]]))

  if (!is.null(pathogens)) {
    missing_p <- setdiff(pathogens, all_pathogens)
    if (length(missing_p) > 0) {
      stop(sprintf(
        "Pathogen(s) not found in marginal_output: %s",
        paste(missing_p, collapse = ", ")
      ))
    }
    all_pathogens <- pathogens
  }

  # -- Restrict to top N most-tested pathogens ---------------------------------
  if (!is.null(top_n_pathogens)) {
    path_tested <- marginal_output$marginal %>%
      dplyr::filter(!!rlang::sym(pathogen_col) %in% all_pathogens) %>%
      dplyr::group_by(!!rlang::sym(pathogen_col)) %>%
      dplyr::summarise(total_tested = sum(n_tested), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total_tested))

    top_paths <- path_tested[[pathogen_col]][
      seq_len(min(top_n_pathogens, nrow(path_tested)))
    ]
    excluded <- setdiff(all_pathogens, top_paths)

    message(sprintf(
      "top_n_pathogens = %d: keeping %d pathogen(s), excluding %d.",
      top_n_pathogens, length(top_paths), length(excluded)
    ))
    if (length(excluded) > 0) {
      message("  Excluded: ", paste(excluded, collapse = ", "))
    }
    message("  Selected (ranked by total isolates tested):")
    print(path_tested[path_tested[[pathogen_col]] %in% top_paths, ],
      row.names = FALSE
    )

    all_pathogens <- top_paths # preserves descending rank order
  }

  out <- list()
  skipped_single <- list()

  n_pathogens_total <- length(all_pathogens)
  if (n_pathogens_total > 1) {
    message(sprintf(
      "\nProcessing %d pathogen(s) with %d core(s)...",
      n_pathogens_total, n_cores
    ))
  }

  # -- Per-pathogen worker (closure over outer env) --------------------------
  .process_one <- function(path_idx) {
    path <- all_pathogens[path_idx]

    if (n_pathogens_total > 1) {
      message(sprintf(
        "\n[%d/%d] Processing '%s'...",
        path_idx, n_pathogens_total, path
      ))
    }

    # -- Determine antibiotic classes ---------------------------------------
    marg_k <- marginal_output$marginal[
      marginal_output$marginal[[pathogen_col]] == path,
    ]

    if (exclude_near_zero) {
      nz_classes <- marginal_output$near_zero[
        marginal_output$near_zero[[pathogen_col]] == path,
        antibiotic_class_col,
        drop = TRUE
      ]
      marg_k <- marg_k[!marg_k[[antibiotic_class_col]] %in% nz_classes, ]
    }

    if (!is.null(top_n_classes) && nrow(marg_k) > top_n_classes) {
      marg_k <- marg_k[order(marg_k$n_tested, decreasing = TRUE), ][
        seq_len(top_n_classes),
      ]
      message(sprintf(
        "'%s': restricting to top %d classes by n_tested.", path, top_n_classes
      ))
    }

    classes <- sort(marg_k[[antibiotic_class_col]])
    n <- length(classes)

    if (n == 0) {
      message(sprintf("'%s': no classes remain after filtering, skipping.", path))
      return(list(type = "empty", path = path))
    }

    if (n < 2) {
      return(list(
        type = "skipped",
        path = path,
        data = data.frame(
          pathogen            = path,
          antibiotic_class    = marg_k[[antibiotic_class_col]],
          n_tested            = marg_k$n_tested,
          marginal_resistance = marg_k$marginal_resistance,
          stringsAsFactors    = FALSE
        )
      ))
    }

    n_profiles <- 2L^n

    if (n > 18L) {
      warning(sprintf(
        "'%s': %d classes -> 2^%d = %d profiles. This may be very slow and memory-intensive.",
        path, n, n, n_profiles
      ))
    }

    # -- Marginal resistance vector r_kd ------------------------------------
    r_marg <- setNames(
      marg_k$marginal_resistance,
      marg_k[[antibiotic_class_col]]
    )[classes]

    # -- Pairwise co-resistance matrix (may be NULL) ------------------------
    co_mat <- if (path %in% names(coresistance_output)) {
      coresistance_output[[path]]$prevalence
    } else {
      NULL
    }

    # -- Enumerate 2^n resistance profiles ---------------------------------
    profiles_enum <- enumerate_binary_profiles(classes)
    n_profiles <- nrow(profiles_enum)
    profile_labels <- profiles_enum$profile_delta
    profiles_mat <- as.matrix(profiles_enum[, classes, drop = FALSE])

    # -- Constraint matrix M (m x 2^n) and target vector v -----------------
    cm <- build_constraint_matrix(profiles_enum, r_marg, co_mat)
    M <- cm$M
    v <- cm$v

    if (length(cm$capped_pairs) > 0L) {
      message(sprintf(
        "'%s': %d pairwise value(s) capped to min(marginal): %s",
        path, length(cm$capped_pairs),
        paste(names(cm$capped_pairs), cm$capped_pairs, sep = " ", collapse = "; ")
      ))
    }
    if (length(cm$fallback_pairs) > 0L) {
      message(sprintf(
        "'%s': %d pair(s) used independence fallback: %s",
        path, length(cm$fallback_pairs),
        paste(cm$fallback_pairs, collapse = "; ")
      ))
    }

    # -- QP: simplex-constrained weighted least-squares ---------------------
    #
    # Minimise  (1/2) p^T H p  -  d^T p
    #   H = (2/sigma_sq) * M^T M  +  ridge * I   (guaranteed pos-def)
    #   d = (2/sigma_sq) * M^T v
    # Subject to: sum(p) = 1, p >= 0
    #
    # OSQP path: constraint matrix A = rbind([1...1], I_{2^n}) stored sparse
    #   -> O(2^n) memory vs O(4^n) for the dense diag() used by quadprog.
    # quadprog path: retained as fallback when osqp/Matrix are unavailable.
    coef <- 2.0 / sigma_sq
    H_mat <- coef * crossprod(M)
    diag(H_mat) <- diag(H_mat) + ridge
    d_qp <- coef * drop(crossprod(M, v))

    p_hat <- tryCatch(
      {
        if (has_osqp) {
          A_sp <- rbind(
            Matrix::Matrix(1.0, nrow = 1L, ncol = n_profiles, sparse = TRUE),
            Matrix::Diagonal(n_profiles)
          )
          prob <- osqp::osqp(
            P = Matrix::forceSymmetric(Matrix::Matrix(H_mat)),
            q = -d_qp,
            A = A_sp,
            l = c(1.0, rep(0.0, n_profiles)),
            u = c(1.0, rep(Inf, n_profiles)),
            pars = osqp::osqpSettings(
              verbose  = FALSE,
              eps_abs  = 1e-8,
              eps_rel  = 1e-8,
              max_iter = 10000L,
              polish   = TRUE
            )
          )
          res <- prob$solve()
          if (!(res$info$status %in% c("solved", "solved_inaccurate"))) {
            stop(paste("OSQP status:", res$info$status))
          }
          pmax(res$x, 0.0)
        } else {
          # quadprog fallback -- dense Amat, slow for n > 12
          Amat <- cbind(rep(1.0, n_profiles), diag(n_profiles))
          bvec <- c(1.0, rep(0.0, n_profiles))
          sol <- quadprog::solve.QP(
            Dmat = H_mat, dvec = d_qp,
            Amat = Amat,  bvec = bvec, meq = 1L
          )
          pmax(sol$solution, 0.0)
        }
      },
      error = function(e) {
        warning(sprintf(
          "'%s': QP solver failed (%s). Returning uniform distribution.",
          path, conditionMessage(e)
        ))
        rep(1.0 / n_profiles, n_profiles)
      }
    )

    p_hat <- p_hat / sum(p_hat)

    # -- Constraint residuals -----------------------------------------------
    residuals <- drop(M %*% p_hat) - v
    names(residuals) <- cm$constraint_names

    profiles_df <- data.frame(
      profile = profile_labels,
      probability = p_hat,
      stringsAsFactors = FALSE
    )
    profiles_df <- cbind(profiles_df, as.data.frame(profiles_mat, check.names = FALSE))

    message(sprintf(
      "'%s': n=%d classes -> %d profiles. Max |residual| = %.5f.",
      path, n, n_profiles, max(abs(residuals))
    ))

    list(
      type = "success",
      path = path,
      result = list(
        profiles              = profiles_df,
        classes               = classes,
        n_classes             = n,
        constraint_residuals  = residuals,
        constraint_targets    = setNames(cm$v, cm$constraint_names),
        constraint_names      = cm$constraint_names
      )
    )
  }

  # -- Execute: parallel or sequential ---------------------------------------
  if (n_cores > 1L) {
    all_results <- parallel::mclapply(
      seq_along(all_pathogens), .process_one,
      mc.cores = n_cores,
      mc.preschedule = FALSE # better load balancing across heterogeneous pathogens
    )
  } else {
    all_results <- lapply(seq_along(all_pathogens), .process_one)
  }

  # -- Collect results --------------------------------------------------------
  for (res in all_results) {
    if (is.null(res) || inherits(res, "try-error")) next
    if (res$type == "success") {
      out[[res$path]] <- res$result
    } else if (res$type == "skipped") {
      skipped_single[[res$path]] <- res$data
    }
  }

  # -- Summary: pathogens skipped because only 1 class remained --------------
  if (length(skipped_single) > 0) {
    skipped_tbl <- do.call(rbind, skipped_single)
    message(sprintf(
      "\n%d pathogen(s) skipped -- only 1 antibiotic class remained after filtering.",
      length(skipped_single)
    ))
    message(paste0(
      "For these pathogens the resistance profile IS the marginal resistance.\n",
      "  Tip: use top_n_classes, reduce zero_threshold, or set exclude_near_zero = FALSE",
      " to include them.\n"
    ))
    print(skipped_tbl, row.names = FALSE)
  }

  return(out)
}


# Normalize String for Joining
#
# Lowercases, trims whitespace, and removes punctuation for fuzzy joins.


# Resistance class selection (moved from prep_ast_and_syndrome.R)
# These are analysis functions that inform DALY burden attribution - they
# do not belong in the preprocessing layer.

#' Select Resistance Class for Burden Attribution
#'
#' Selects a single resistance class per event using beta-lactam hierarchy
#' and relative risk (RR) values. Prevents double-counting in DALY burden
#' estimation by choosing the most clinically relevant resistant class.
#'
#' Selection order:
#' 1. Beta-lactam hierarchy rank (Carbapenems > 4GC > 3GC > ...)
#' 2. RR value (higher relative risk = higher priority)
#' 3. Alphabetical tie-breaker
#'
#' @param data Data frame with class-level resistance and RR columns.
#' @param event_col Character. Event ID column. Default "event_id".
#' @param class_col Character. Antibiotic class column. Default "antibiotic_class".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param rr_col Character. RR value column. Default "rr_value".
#'   If missing, only hierarchy is used.
#' @param hierarchy Named numeric vector. Custom hierarchy (class name -> rank).
#'   If NULL, uses \code{get_beta_lactam_hierarchy()}.
#' @param filter_resistant Logical. If TRUE, only consider resistant (R) classes.
#'   Default TRUE.
#'
#' @return Data frame filtered to one resistance class per event.
#' @export
select_resistance_class <- function(data,
                                    event_col = "event_id",
                                    class_col = "antibiotic_class",
                                    susceptibility_col = "antibiotic_value",
                                    rr_col = "rr_value",
                                    hierarchy = NULL,
                                    filter_resistant = TRUE) {
  required_cols <- c(event_col, class_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (is.null(hierarchy)) hierarchy <- get_beta_lactam_hierarchy()

  use_rr <- rr_col %in% names(data)
  if (!use_rr) {
    message(sprintf("RR column '%s' not found. Using hierarchy only for selection.", rr_col))
  }

  n_before <- nrow(data)
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  message(sprintf(
    "Selecting resistance classes using hierarchy%s...",
    ifelse(use_rr, " + RR", "")
  ))

  if (filter_resistant) {
    data <- data %>% dplyr::filter(!!rlang::sym(susceptibility_col) == "R")
    message(sprintf("Filtered to resistant classes: %d rows", nrow(data)))
  }

  selected <- prioritize_resistance(
    data      = data,
    event_col = event_col,
    class_col = class_col,
    rr_col    = if (use_rr) rr_col else NULL,
    hierarchy = hierarchy
  )

  n_events_after <- dplyr::n_distinct(selected[[event_col]])
  message(sprintf(
    "Selected: %d rows from %d events (avg %.2f classes/event before -> 1.0 after)",
    nrow(selected), n_events_after, n_before / n_events_before
  ))

  selection_summary <- selected %>%
    dplyr::count(!!rlang::sym(class_col), name = "n_events") %>%
    dplyr::arrange(dplyr::desc(n_events)) %>%
    utils::head(10)
  message("\nTop 10 selected classes:")
  print(selection_summary)

  return(selected)
}


#' Prioritize Resistance Class
#'
#' Internal helper for \code{select_resistance_class()}. Applies beta-lactam
#' hierarchy + RR ranking to select one row per event.
#'
#' @param data Data frame.
#' @param event_col Character. Event ID column.
#' @param class_col Character. Class column.
#' @param rr_col Character or NULL. RR column; if NULL uses hierarchy only.
#' @param hierarchy Named numeric vector mapping class name to rank.
#'
#' @return Data frame with one row per event (highest priority class selected).
#' @keywords internal
prioritize_resistance <- function(data,
                                  event_col,
                                  class_col,
                                  rr_col = NULL,
                                  hierarchy) {
  # Accept either a named vector (name = class, value = rank) or an unnamed
  # character vector (position = rank) -- get_beta_lactam_hierarchy() returns
  # the latter.
  hier_classes <- if (!is.null(names(hierarchy))) names(hierarchy) else as.character(hierarchy)
  hierarchy_df <- data.frame(
    class = hier_classes,
    hierarchy_rank = seq_along(hierarchy),
    stringsAsFactors = FALSE
  )
  names(hierarchy_df)[1] <- class_col

  data <- data %>% dplyr::left_join(hierarchy_df, by = class_col)
  max_rank <- max(hierarchy_df$hierarchy_rank, na.rm = TRUE)
  data <- data %>%
    dplyr::mutate(hierarchy_rank = dplyr::coalesce(hierarchy_rank, max_rank + 1L))

  if (!is.null(rr_col) && rr_col %in% names(data)) {
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(
        hierarchy_rank, dplyr::desc(!!rlang::sym(rr_col)),
        !!rlang::sym(class_col)
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(selection_method = "hierarchy_rr", selection_confidence = "high")
  } else {
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(hierarchy_rank, !!rlang::sym(class_col)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(selection_method = "hierarchy_only", selection_confidence = "medium")
  }

  selected <- selected %>% dplyr::select(-hierarchy_rank)

  multi_class_events <- data %>%
    dplyr::group_by(!!rlang::sym(event_col)) %>%
    dplyr::summarise(n_classes = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_classes > 1)

  if (nrow(multi_class_events) > 0) {
    message(sprintf(
      "Applied selection to %d events with multiple resistant classes.",
      nrow(multi_class_events)
    ))
  }

  return(selected)
}


# ===========================================================================
# Pathway 2 -- Bayesian hierarchical multivariate probit (revised)
#
# Key design decisions in this implementation:
#
# 1. Genuine multivariate-probit likelihood: the residual correlation matrix
#    Omega is identified from data via per-event latent-noise augmentation
#    (eta). Each event contributes a D-dim latent noise vector; the binary
#    outcomes constrain those latents and thereby inform L_Omega.
#    Previously L_Omega appeared only in a prior and was never updated.
#
# 2. Correlated random effects: hospital (and lower-level) effects use
#    diag(tau) * L_corr * z_raw reparameterisation so the D-dim random
#    effect vectors are correlated across antibiotic classes.
#
# 3. Separate-pathogen design: pass `pathogen` to restrict the fit to a
#    single pathogen. Orchestrate multi-pathogen fits in the analysis repo.
#
# 4. Configurable hierarchy: 1, 2, or 3 random-effect levels supported
#    (hospital; hospital + patient; hospital + patient + admission).
#    Nested RE groups receive globally unique composite keys built here.
#
# 5. Explicit estimand: "observed_stewardship_event_mix" -- the profile
#    distribution over the actual event case-mix in the data. Standardised
#    estimands require a reference covariate distribution and are for later.
#
# 6. Complete profile enumeration: all 2^D profiles appear in output. Profiles
#    not sampled in finite Monte Carlo receive an estimated probability of 0,
#    which may reflect simulation error for rare profiles rather than true
#    structural zero probability.
#
# 7. Strengthened diagnostics: Rhat < 1.01; tail ESS; tree-depth saturation;
#    E-BFMI; data-volume summary per relevant parameter block.
#
# Missing-AST handling: Observed AST cells impose sign constraints on the latent
# variables z_aug; untested cells (NA) impose no sign constraint and are
# represented as unconstrained latent values. This is NOT equivalent to dropping
# NA cells: the full N_events x D latent matrix is always estimated. The model
# assumes testing is conditionally ignorable given included covariates and random
# effects and does NOT correct for selective testing bias or model the AST
# cascade-testing process.
#
# Estimand: "observed_stewardship_event_mix" -- the joint resistance-profile
# distribution characterising the events recorded in the monitored beds.
# This is NOT the whole-hospital profile distribution unless the stewardship
# beds are a census of all infection events. This distinction is essential for
# the AMR burden pipeline.
# ===========================================================================


# Internal helpers

#' Check pairwise co-testing overlap per hospital x class-pair
#'
#' Reports the full 2x2 co-testing breakdown (n_RR/n_RS/n_SR/n_SS), not just
#' the total co-tested count. A pair can clear \code{min_pairwise_cotested}
#' while showing no variation at all (e.g. 100 co-tested events that are all
#' resistant-resistant) -- such a pair carries no information about residual
#' cross-class correlation, so \code{sufficient} additionally requires all
#' four cells to have at least \code{min_pairwise_cell} observations.
#' @keywords internal
.check_pairwise_cotesting <- function(event_class_data, class_cols, upper_re_col,
                                      min_pairwise_cotested = 30L,
                                      min_pairwise_cell = 1L) {
  rows <- list()
  hospitals <- sort(unique(event_class_data[[upper_re_col]]))
  for (h in hospitals) {
    sub <- event_class_data[event_class_data[[upper_re_col]] == h, class_cols,
      drop = FALSE
    ]
    pairs <- utils::combn(class_cols, 2L, simplify = FALSE)
    for (pr in pairs) {
      c1 <- pr[1L]
      c2 <- pr[2L]
      v1 <- sub[[c1]]
      v2 <- sub[[c2]]
      both_tested <- !is.na(v1) & !is.na(v2)
      n_cotested <- sum(both_tested)
      n_RR <- sum(both_tested & v1 == 1L & v2 == 1L)
      n_RS <- sum(both_tested & v1 == 1L & v2 == 0L)
      n_SR <- sum(both_tested & v1 == 0L & v2 == 1L)
      n_SS <- sum(both_tested & v1 == 0L & v2 == 0L)
      has_full_cell_variation <-
        n_RR >= min_pairwise_cell && n_RS >= min_pairwise_cell &&
          n_SR >= min_pairwise_cell && n_SS >= min_pairwise_cell
      rows[[length(rows) + 1L]] <- tibble::tibble(
        !!upper_re_col := h,
        class_1 = c1,
        class_2 = c2,
        n_cotested = n_cotested,
        n_RR = n_RR,
        n_RS = n_RS,
        n_SR = n_SR,
        n_SS = n_SS,
        has_full_cell_variation = has_full_cell_variation,
        sufficient = (n_cotested >= min_pairwise_cotested) && has_full_cell_variation
      )
    }
  }
  dplyr::bind_rows(rows)
}


#' Check panel eligibility thresholds per hospital x class cell
#' @keywords internal
.validate_panel_eligibility <- function(
  event_class_data,
  class_cols,
  upper_re_col,
  pathogen_col,
  min_tested = 30L,
  min_resistant = 5L,
  min_susceptible = 5L
) {
  rows <- list()
  pathogens <- sort(unique(event_class_data[[pathogen_col]]))
  hospitals <- sort(unique(event_class_data[[upper_re_col]]))

  for (h in hospitals) {
    for (g in pathogens) {
      sub <- event_class_data[
        event_class_data[[upper_re_col]] == h &
          event_class_data[[pathogen_col]] == g, ,
        drop = FALSE
      ]
      if (nrow(sub) == 0L) next
      for (cc in class_cols) {
        vals <- sub[[cc]]
        n_tested <- sum(!is.na(vals))
        n_res <- sum(!is.na(vals) & vals == 1L)
        n_sus <- sum(!is.na(vals) & vals == 0L)
        eligible <- (n_tested >= min_tested) &
          (n_res >= min_resistant) &
          (n_sus >= min_susceptible)
        rows[[length(rows) + 1L]] <- tibble::tibble(
          !!upper_re_col := h,
          !!pathogen_col := g,
          antibiotic_class = cc,
          n_tested = n_tested,
          n_resistant = n_res,
          n_susceptible = n_sus,
          eligible = eligible
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}


#' Resolve the antibiotic-class panel used for profile enumeration
#'
#' Defines the set of classes entering the \eqn{2^D} profile enumeration for
#' one hospital x pathogen pair from the fit-time eligibility report
#' (\code{fitted_model$eligibility_report}), rather than "tested at least
#' once". A class is included only if its marginal
#' n_tested/n_resistant/n_susceptible thresholds were met at that hospital for
#' that pathogen -- this applies to \strong{every} residual structure,
#' including identity, so an identity-residual panel is never narrowed for
#' reasons that only matter for estimating a correlated residual matrix. For
#' correlated-residual fits specifically, classes that lack sufficient
#' pairwise co-testing with every other candidate class at that hospital are
#' additionally dropped iteratively (the class involved in the most
#' insufficient pairs is dropped first, repeated until all remaining pairs
#' clear the threshold), since \eqn{\Omega} for an under-co-tested pair would
#' otherwise be prior-dominated.
#'
#' The eligibility report passed in is computed by
#' \code{fit_bayesian_multivariate_probit()} on \code{event_data} \emph{after}
#' all experiment-specific filtering (pathogen filter, \code{eligible_pairs}
#' semi-join, all-NA-event drop) -- i.e. on the exact row population that was
#' fitted, not on an earlier unfiltered wide table. See the eligibility-report
#' construction in \code{fit_bayesian_multivariate_probit()}.
#'
#' @return List with:
#'   \item{classes}{Character vector of \code{class_cols} entering the panel,
#'     in \code{class_cols} order.}
#'   \item{excluded}{Character vector of \code{class_cols} NOT entering the
#'     panel (empty if none).}
#'   \item{method}{Character; \code{"marginal_only"} (identity residual, or no
#'     pairwise report available), \code{"marginal_plus_pairwise"} (correlated
#'     residual), or \code{"no_eligibility_report_available"} (defensive
#'     fallback for fit objects built some other way).}
#'   \item{reason}{Character or \code{NA}; \code{NA} when nothing was
#'     excluded, otherwise a human-readable breakdown of which classes were
#'     excluded and why (insufficient marginal support vs. insufficient
#'     pairwise co-testing).}
#' @keywords internal
.resolve_profile_class_panel <- function(
  class_cols,
  hospital,
  pathogen,
  eligibility_report,
  upper_re_col,
  pathogen_col,
  residual_structure = "identity"
) {
  marginal <- eligibility_report$marginal
  pairwise <- eligibility_report$pairwise

  if (is.null(marginal) || nrow(marginal) == 0L) {
    # No eligibility report available (fit_bayesian_multivariate_probit() always
    # computes one; this only fires for fit objects built some other way).
    # Fall back to the full class set rather than silently treating untested
    # classes as eligible.
    return(list(
      classes = class_cols, excluded = character(0L),
      method = "no_eligibility_report_available", reason = NA_character_
    ))
  }

  sub <- marginal[marginal[[upper_re_col]] == hospital &
    marginal[[pathogen_col]] == pathogen, , drop = FALSE]
  elig_classes <- sub$antibiotic_class[sub$eligible]
  elig_classes <- intersect(class_cols, elig_classes) # preserve class_cols order
  marginal_excluded <- setdiff(class_cols, elig_classes)

  use_pairwise <- identical(residual_structure, "correlated") &&
    !is.null(pairwise) && nrow(pairwise) > 0L
  method <- if (use_pairwise) "marginal_plus_pairwise" else "marginal_only"

  if (length(elig_classes) == 0L) {
    reason <- sprintf("insufficient_marginal_support: %s", paste(marginal_excluded, collapse = ", "))
    return(list(classes = character(0L), excluded = marginal_excluded, method = method, reason = reason))
  }

  if (!use_pairwise) {
    reason <- if (length(marginal_excluded) == 0L) {
      NA_character_
    } else {
      sprintf("insufficient_marginal_support: %s", paste(marginal_excluded, collapse = ", "))
    }
    return(list(classes = elig_classes, excluded = marginal_excluded, method = method, reason = reason))
  }

  # Correlated residual only: iteratively drop the class most often paired
  # with an insufficiently co-tested partner, until every remaining pair
  # clears the pairwise co-testing threshold (or fewer than 2 classes remain).
  psub <- pairwise[pairwise[[upper_re_col]] == hospital, , drop = FALSE]
  current <- elig_classes
  repeat {
    if (length(current) < 2L) break
    pr <- psub[psub$class_1 %in% current & psub$class_2 %in% current, , drop = FALSE]
    bad <- pr[!pr$sufficient, , drop = FALSE]
    if (nrow(bad) == 0L) break
    offender_counts <- table(c(bad$class_1, bad$class_2))
    worst <- names(offender_counts)[which.max(offender_counts)]
    current <- setdiff(current, worst)
  }
  pairwise_excluded <- setdiff(elig_classes, current)
  all_excluded <- union(marginal_excluded, pairwise_excluded)

  reason_parts <- c(
    if (length(marginal_excluded) > 0L) {
      sprintf("insufficient_marginal_support: %s", paste(marginal_excluded, collapse = ", "))
    },
    if (length(pairwise_excluded) > 0L) {
      sprintf("insufficient_pairwise_cotesting: %s", paste(pairwise_excluded, collapse = ", "))
    }
  )
  reason <- if (length(reason_parts) == 0L) NA_character_ else paste(reason_parts, collapse = "; ")

  list(classes = current, excluded = all_excluded, method = method, reason = reason)
}


#' Summarize a posterior D x D correlation matrix from stored fit draws
#'
#' Extracts posterior mean/median/CI for each off-diagonal cell of a named
#' correlation-matrix generated quantity (e.g. "Omega", "R_hospital",
#' "R_patient", "R_admission"). These are three distinct quantities and must
#' not be combined into one table: Omega is the event-level residual
#' cross-class correlation (only estimated when residual_structure ==
#' "correlated"); R_hospital/R_patient/R_admission are the correlation, across
#' antibiotic classes, of that random effect's own tendencies (estimated
#' regardless of residual_structure, since the hospital/patient/admission
#' random effects always use a diag_pre_multiply(tau, L_corr) parameterisation).
#' Returns NULL if the requested matrix_var is not present in fit$draws (e.g.
#' Omega on an identity-residual fit, or R_admission on a 1- or 2-RE fit).
#'
#' @param fit A fitted_model object (or its lightweight saved form) with a
#'   \code{$draws} posterior::draws_array retaining the requested matrix_var.
#' @param matrix_var Character. One of "Omega", "R_hospital", "R_patient",
#'   "R_admission".
#' @param class_cols Character vector of antibiotic class column names, in
#'   the same order used to fit the model (dimension order of matrix_var).
#' @param ci_level Numeric credible-interval level. Default 0.95.
#' @return A tibble with one row per off-diagonal class pair, or NULL if
#'   matrix_var was not estimated for this fit.
#' @export
#' @param block_index Integer or NULL. When the fit used the generic
#'   random-effect architecture (Stage 1), per-block correlation matrices are
#'   emitted as \code{R_block[r,i,j]} (block index r as the FIRST subscript),
#'   not \code{R_hospital[i,j]}/\code{R_patient[i,j]}/etc. Pass
#'   \code{matrix_var = "R_block"} and the 1-based block index here to look up
#'   a specific block's correlation matrix; leave NULL for 2-index variables
#'   (\code{Omega}, or legacy \code{R_hospital} etc. from a fit made with the
#'   old hardcoded 1re/2re/3re Stan models).
#' @export
summarize_fit_correlation_matrix <- function(fit, matrix_var, class_cols, ci_level = 0.95,
                                             block_index = NULL) {
  draws_arr <- fit$draws
  if (is.null(draws_arr)) {
    return(NULL)
  }
  var_names <- posterior::variables(draws_arr)
  D <- length(class_cols)
  alpha <- (1 - ci_level) / 2
  rows <- list()
  for (i in seq_len(D)) {
    for (j in seq_len(D)) {
      if (i >= j) next # off-diagonal, upper triangle only (matrix is symmetric)
      vname <- if (is.null(block_index)) {
        sprintf("%s[%d,%d]", matrix_var, i, j)
      } else {
        sprintf("%s[%d,%d,%d]", matrix_var, block_index, i, j)
      }
      if (!vname %in% var_names) next
      draws_ij_sub <- posterior::subset_draws(draws_arr, variable = vname)
      draws_ij <- as.numeric(posterior::extract_variable(draws_arr, vname))
      rows[[length(rows) + 1L]] <- tibble::tibble(
        class_1             = class_cols[i],
        class_2             = class_cols[j],
        correlation_mean    = mean(draws_ij),
        correlation_median  = stats::median(draws_ij),
        correlation_lower   = stats::quantile(draws_ij, probs = alpha, names = FALSE),
        correlation_upper   = stats::quantile(draws_ij, probs = 1 - alpha, names = FALSE),
        # Per-pair sampler diagnostics -- a visually attractive correlation
        # estimate is meaningless if the chains never agreed on it. Computed
        # from the SAME chain-preserving subset (not the flattened draws_ij
        # above) so rhat/ess_bulk use the actual multi-chain structure.
        rhat                = tryCatch(posterior::rhat(draws_ij_sub), error = function(e) NA_real_),
        ess_bulk            = tryCatch(posterior::ess_bulk(draws_ij_sub), error = function(e) NA_real_),
        ess_tail            = tryCatch(posterior::ess_tail(draws_ij_sub), error = function(e) NA_real_)
      )
    }
  }
  if (length(rows) == 0L) {
    return(NULL)
  }
  dplyr::bind_rows(rows)
}


# Generic Stan models -- arbitrary number (R) of random-intercept blocks.
# Replaces the six hardcoded 1re/2re/3re x identity/correlated variants with
# two variants (identity/correlated residual only), looping over blocks via
# the flattened total_re_levels/n_levels/level_start/re_idx representation
# built by prepare_random_effects().
#
# Random-effect parameterisation per block r, exactly matching the legacy
# per-level diag_pre_multiply(tau, L_corr) * z convention:
#   re_effect[, level_start[r]:(level_start[r]+n_levels[r]-1)] =
#     diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, <same slice>]
# Each event's total random-effect contribution is the SUM across all R
# blocks of that block's effect at the event's level in that block:
#   mu_mat[e] = X_event[e] * beta + sum_r re_effect[, re_idx[e, r]]'
# re_idx[e, r] is already the FLATTENED (global) level index for block r,
# so R-side code only needs one generic column lookup per block, not a
# separate hardcoded h_ev/p_ev/a_ev array per level.

.amr_probit_stan_generic_correlated <- function() {
  r"(
// Multivariate probit with data augmentation -- generic number (R) of
// random-intercept blocks, correlated residual (Omega estimated).
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> K;
  int<lower=1> R;
  int<lower=1> total_re_levels;
  array[R] int<lower=1> n_levels;
  array[R] int<lower=1> level_start;
  array[N_events, R] int<lower=1, upper=total_re_levels> re_idx;

  matrix[N_events, K] X_event;

  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  array[N] real<lower=-1,upper=1> sign_obs; // +1 resistant, -1 susceptible
  array[N] int<lower=1,upper=N_events*D> obs_idx; // matches to_vector()'s column-major order

  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta_residual;
  real<lower=1> lkj_eta_random_effect;
}
parameters {
  matrix[K, D] beta;
  matrix[D, total_re_levels] z_re;
  array[R] vector<lower=0>[D] tau_re;
  array[R] cholesky_factor_corr[D] L_corr_re;
  cholesky_factor_corr[D] L_Omega;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, total_re_levels] re_effect;
  for (r in 1:R) {
    int lo = level_start[r];
    int hi = level_start[r] + n_levels[r] - 1;
    re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi];
  }
}
model {
  to_vector(beta) ~ normal(0, prior_beta_sd);
  for (r in 1:R) {
    int lo = level_start[r];
    int hi = level_start[r] + n_levels[r] - 1;
    tau_re[r]                    ~ normal(0, prior_tau_sd);
    to_vector(z_re[, lo:hi])     ~ std_normal();
    L_corr_re[r]                 ~ lkj_corr_cholesky(lkj_eta_random_effect);
  }
  L_Omega ~ lkj_corr_cholesky(lkj_eta_residual);

  {
    // exp() only at the N observed cells -- masking after exp() risks 0 * Inf = NaN.
    matrix[N_events, D] z_aug = z_free;
    for (n in 1:N) {
      z_aug[ev_idx[n], d_idx[n]] = sign_obs[n] * exp(z_free[ev_idx[n], d_idx[n]]);
    }

    matrix[N_events, D] mu_mat = X_event * beta;
    for (e in 1:N_events) {
      for (r in 1:R) {
        mu_mat[e] += re_effect[, re_idx[e, r]]';
      }
    }

    // Batch events under the shared L_Omega. Keep these local to exclude them from CmdStan output.
    array[N_events] vector[D] z_aug_arr;
    array[N_events] vector[D] mu_arr;
    for (e in 1:N_events) {
      z_aug_arr[e] = z_aug[e]';
      mu_arr[e] = mu_mat[e]';
    }
    target += multi_normal_cholesky_lupdf(z_aug_arr | mu_arr, L_Omega);
  }
  target += sum(to_vector(z_free)[obs_idx]);
}
generated quantities {
  corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L_Omega);
  array[R] corr_matrix[D] R_block;
  for (r in 1:R) R_block[r] = multiply_lower_tri_self_transpose(L_corr_re[r]);
}
)"
}

.amr_probit_stan_generic_identity <- function() {
  r"(
// Multivariate probit -- generic number (R) of random-intercept blocks,
// identity residual (classes conditionally independent; L_Omega not
// estimated, residual covariance = I_D).
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> K;
  int<lower=1> R;
  int<lower=1> total_re_levels;
  array[R] int<lower=1> n_levels;
  array[R] int<lower=1> level_start;
  array[N_events, R] int<lower=1, upper=total_re_levels> re_idx;

  matrix[N_events, K] X_event;

  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  array[N] real<lower=-1,upper=1> sign_obs; // +1 resistant, -1 susceptible
  array[N] int<lower=1,upper=N_events*D> obs_idx; // matches to_vector()'s column-major order

  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta_random_effect;
}

parameters {
  matrix[K, D] beta;
  matrix[D, total_re_levels] z_re;
  array[R] vector<lower=0>[D] tau_re;
  array[R] cholesky_factor_corr[D] L_corr_re;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, total_re_levels] re_effect;
  for (r in 1:R) {
    int lo = level_start[r];
    int hi = level_start[r] + n_levels[r] - 1;
    re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi];
  }
}
model {
  to_vector(beta) ~ normal(0, prior_beta_sd);
  for (r in 1:R) {
    int lo = level_start[r];
    int hi = level_start[r] + n_levels[r] - 1;
    tau_re[r]                    ~ normal(0, prior_tau_sd);
    to_vector(z_re[, lo:hi])     ~ std_normal();
    L_corr_re[r]                 ~ lkj_corr_cholesky(lkj_eta_random_effect);
  }

  {
    // exp() only at the N observed cells -- masking after exp() risks 0 * Inf = NaN.
    matrix[N_events, D] z_aug = z_free;
    for (n in 1:N) {
      z_aug[ev_idx[n], d_idx[n]] = sign_obs[n] * exp(z_free[ev_idx[n], d_idx[n]]);
    }

    matrix[N_events, D] mu_mat = X_event * beta;
    for (e in 1:N_events) {
      for (r in 1:R) {
        mu_mat[e] += re_effect[, re_idx[e, r]]';
      }
    }

    // Keep these local to exclude them from CmdStan output.
    target += normal_lupdf(to_vector(z_aug) | to_vector(mu_mat), 1.0);
  }
  target += sum(to_vector(z_free)[obs_idx]);
}
generated quantities {
  array[R] corr_matrix[D] R_block;
  for (r in 1:R) R_block[r] = multiply_lower_tri_self_transpose(L_corr_re[r]);
}
)"
}

# Fixed-effects-only variants deliberately omit every random-effect parameter.
# Keeping these separate from the generic R >= 1 models avoids relying on
# zero-sized Stan arrays/matrices and leaves existing model mathematics intact.
.amr_probit_stan_fixed_correlated <- function() {
  r"(
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> K;
  matrix[N_events, K] X_event;
  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  array[N] real<lower=-1,upper=1> sign_obs;
  array[N] int<lower=1,upper=N_events*D> obs_idx;
  real<lower=0> prior_beta_sd;
  real<lower=1> lkj_eta_residual;
}
parameters {
  matrix[K, D] beta;
  cholesky_factor_corr[D] L_Omega;
  matrix[N_events, D] z_free;
}
model {
  to_vector(beta) ~ normal(0, prior_beta_sd);
  L_Omega ~ lkj_corr_cholesky(lkj_eta_residual);
  {
    matrix[N_events, D] z_aug = z_free;
    for (n in 1:N)
      z_aug[ev_idx[n], d_idx[n]] = sign_obs[n] * exp(z_free[ev_idx[n], d_idx[n]]);
    array[N_events] vector[D] z_aug_arr;
    array[N_events] vector[D] mu_arr;
    matrix[N_events, D] mu_mat = X_event * beta;
    for (e in 1:N_events) {
      z_aug_arr[e] = z_aug[e]';
      mu_arr[e] = mu_mat[e]';
    }
    target += multi_normal_cholesky_lupdf(z_aug_arr | mu_arr, L_Omega);
  }
  target += sum(to_vector(z_free)[obs_idx]);
}
generated quantities {
  corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
)"
}

.amr_probit_stan_fixed_identity <- function() {
  r"(
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> K;
  matrix[N_events, K] X_event;
  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  array[N] real<lower=-1,upper=1> sign_obs;
  array[N] int<lower=1,upper=N_events*D> obs_idx;
  real<lower=0> prior_beta_sd;
}
parameters {
  matrix[K, D] beta;
  matrix[N_events, D] z_free;
}
model {
  to_vector(beta) ~ normal(0, prior_beta_sd);
  {
    matrix[N_events, D] z_aug = z_free;
    for (n in 1:N)
      z_aug[ev_idx[n], d_idx[n]] = sign_obs[n] * exp(z_free[ev_idx[n], d_idx[n]]);
    matrix[N_events, D] mu_mat = X_event * beta;
    target += normal_lupdf(to_vector(z_aug) | to_vector(mu_mat), 1.0);
  }
  target += sum(to_vector(z_free)[obs_idx]);
}
)"
}


# fit_bayesian_multivariate_probit()

#' Fit Bayesian Hierarchical Multivariate Probit Model for Resistance Profiles
#'
#' Accepts wide-format event-level AST data (one row per organism-event, one
#' column per antibiotic class with 0 / 1 / \code{NA} values) and fits a
#' hierarchical Bayesian multivariate probit model via \pkg{cmdstanr}.
#'
#' \strong{Model (per pathogen):}
#' \deqn{Y_{ed} = \mathbf{1}(Z_{ed} > 0)}
#' \deqn{Z_{ed} = \mathbf{x}_e^\top\beta_d
#'   + \text{hospital\_effect}_{d,h(e)}
#'   \;[+\; \text{patient\_effect}_{d,p(e)}]
#'   \;[+\; \text{admission\_effect}_{d,a(e)}]
#'   + (L_\Omega\,\eta_e)_d}
#' \deqn{\eta_e \sim N(0, I_D),\quad L_\Omega \sim \text{LKJCholesky}(\eta)}
#'
#' The per-event latent noise \eqn{\eta_e} makes \eqn{L_\Omega} (and therefore
#' \eqn{\Omega = L_\Omega L_\Omega^\top}) identifiable from the observed binary
#' outcomes. Hospital (and lower-level) effects use a
#' \eqn{\text{diag}(\tau) L_{\text{corr}} z_{\text{raw}}} parameterisation so
#' the D-dim random effect vectors are correlated across antibiotic classes.
#'
#' \strong{Single-pathogen design:} pass \code{pathogen} to restrict the fit.
#' Run once per pathogen and orchestrate in the analysis repository.
#'
#' \strong{Random effects} (\code{random_effects}): an optional character
#' vector or list of grouping blocks defining the clustering hierarchy. Use
#' \code{list()} for a fixed-effects-only model. When blocks are present, the
#' first element is the upper-most level;
#' subsequent elements are nested within it. Nested levels receive globally unique
#' composite keys built internally. Any hierarchical grouping variable can occupy
#' any slot -- the labels hospital/patient/admission are semantic conventions, not
#' constraints. Isolate- or sample-event-level effects can be passed as the second
#' or third element, but note the returned object uses generic level names
#' (\code{upper_re_col}, \code{middle_re_col}, \code{lower_re_col}).
#'
#' \strong{Missing-AST handling:} Observed AST cells impose sign constraints on
#' the latent variables; untested cells impose no sign constraint and are
#' represented as unconstrained latent values. The model assumes testing is
#' conditionally ignorable given included covariates and random effects. It does
#' NOT correct for selective testing bias. See \code{residual_structure} for
#' control over residual covariance complexity.
#'
#' \strong{Fixed-effect missingness:} the function warns but does NOT silently
#' impute covariates. Imputation is the caller's responsibility; columns with
#' any remaining NA values after the call will cause Stan to fail.
#'
#' @param event_class_data Data frame. One row per organism-event. Antibiotic
#'   class columns hold 0 (susceptible), 1 (resistant), or \code{NA} (not
#'   tested). Metadata columns hold covariates and grouping variables.
#' @param class_cols Character vector. Names of the antibiotic class columns.
#'   Required.
#' @param fixed_effects Character vector. Event-level covariate column names.
#'   Required.
#' @param random_effects Character vector or named list of random-intercept
#'   blocks. Use \code{list()} for no random-effect blocks. When non-empty,
#'   the legacy character-vector form names grouping columns; the list form
#'   uses \code{name}, \code{group_col}, and optional \code{terms = "intercept"}.
#' @param profile_group_col Character scalar or \code{NULL}. Column used for
#'   downstream profile aggregation, eligibility summaries, and validation.
#'   It is independent of the random-effect specification. Defaults to the
#'   first random-effect grouping column when random effects exist; it is
#'   required for fixed-effects-only fits.
#' @param pathogen Character or \code{NULL}. When supplied, filters
#'   \code{event_class_data} to rows where \code{pathogen_col} equals this
#'   value before fitting. Recommended: fit one pathogen at a time.
#' @param pathogen_col Character. Column identifying the pathogen.
#'   Default \code{"pathogen"}.
#' @param event_id_col Character. Column in \code{event_class_data} holding
#'   unique event identifiers. Default \code{"event_id"}.
#' @param eligible_pairs Tibble or \code{NULL}. Hospital x pathogen pairs to
#'   include. \code{NULL} uses all pairs present in the data.
#' @param outcome_col Character or \code{NULL}. Patient outcome column. Only
#'   used downstream to split R_ALL and R_NF cohorts -- does not enter the
#'   probit likelihood. Default \code{NULL}.
#' @param reserve_drug_cols Character vector or \code{NULL}. Class columns to
#'   exclude from the main model.
#' @param panel_eligibility Named list. Eligibility thresholds per
#'   hospital x class cell: \code{min_tested} (default 30),
#'   \code{min_resistant} (default 5), \code{min_susceptible} (default 5).
#'   Cells not meeting thresholds are reported but fitting still proceeds.
#' @param residual_structure Character. Controls the residual covariance structure.
#'   \code{"identity"} (default): classes are conditionally independent given
#'   fixed and random effects -- residual covariance = I_D. \code{"correlated"}:
#'   estimates a full residual correlation matrix \eqn{\Omega} via LKJCholesky
#'   prior. Use \code{"correlated"} only when panel co-testing overlap is adequate
#'   (check \code{$eligibility_report$pairwise}); otherwise \eqn{\Omega} is driven
#'   mainly by the LKJ prior and adds identifiability risk.
#' @param estimand Character. Identifies the target quantity. Only
#'   \code{"observed_stewardship_event_mix"} is currently supported.
#' @param prior_config Named list. Any subset of \code{beta_sd} (default 1.5),
#'   \code{tau_sd} (default 1.0), legacy \code{lkj_eta} (default 2.0), or
#'   separate \code{lkj_eta_residual} and \code{lkj_eta_random_effect}.
#'   Legacy \code{lkj_eta} resolves to both values for exact backward compatibility.
#' @param sampler_config Named list. Sampler settings forwarded to
#'   \code{cmdstanr::sample()}: \code{chains} (4), \code{iter_warmup} (1000),
#'   \code{iter_sampling} (1000), \code{adapt_delta} (NULL, uses Stan default),
#'   \code{max_treedepth} (NULL), \code{seed} (123), \code{parallel_chains}
#'   (NULL), \code{max_param_count} (NULL -- set to a positive integer to stop
#'   if approximate parameter count exceeds the threshold),
#'   \code{save_warmup} (NULL, i.e. \code{cmdstanr}'s own default of
#'   \code{FALSE} -- set \code{TRUE} to retain warmup-phase draws in the raw
#'   CmdStan CSV for later adaptation-trajectory forensics; note this recovers
#'   per-iteration diagnostics and parameter draws during warmup, NOT the
#'   full evolving HMC mass matrix at every adaptation window -- only the
#'   final adapted metric is ever written to the CSV, regardless of this
#'   setting), \code{metric} (NULL, i.e. \code{cmdstanr}'s default
#'   \code{"diag_e"} -- set \code{"dense_e"} to adapt a full mass matrix
#'   instead of a diagonal one; substantially more expensive to adapt at
#'   high parameter counts, e.g. models with a large \code{z_free} block),
#'   \code{init} (NULL, i.e. CmdStan's default Uniform(-2,2) unconstrained
#'   init -- set to a value accepted by \code{cmdstanr::sample()}'s own
#'   \code{init} argument, e.g. a list of per-chain initial-value lists, to
#'   run an initialization-sensitivity check). Any additional entries are
#'   forwarded via \code{...}.
#' @param compute Named list controlling the Stan execution backend. Supported
#'   fields are \code{backend} (\code{"cpu"} default or \code{"opencl"}),
#'   \code{opencl_platform_id}, \code{opencl_device_id}, and
#'   \code{allow_cpu_fallback} (\code{FALSE} by default). OpenCL changes the
#'   compilation/sampling backend only; it does not change the statistical
#'   model, priors, diagnostics, or downstream estimands.
#' @param show_messages Logical. Print sampling progress. Default \code{TRUE}.
#' @param save_full_latent_diagnostics Logical. When \code{FALSE} (default),
#'   \code{diagnostics_detail$all_parameters} (per-parameter Rhat/ESS
#'   including the \code{N_events * D} latent \code{z_free} nuisance
#'   parameters -- tens of thousands of rows for realistic data) is returned
#'   with 0 rows to keep the fitted object small; the full-scope summary
#'   scalars (\code{max_rhat_full}, \code{min_ess_bulk_full},
#'   \code{min_ess_tail_full}, \code{n_params_full},
#'   \code{pct_rhat_full_above_1_01}, \code{pct_ess_bulk_full_below_100}) are
#'   always populated in \code{diagnostics} regardless. Set \code{TRUE} to
#'   retain the full per-parameter table for debugging.
#' @param ... Additional arguments forwarded to \code{cmdstanr::sample()}.
#'
#' @return Named list with elements: \code{draws}, \code{diagnostics},
#'   \code{diagnostics_detail}, \code{fit}, \code{data_long}, \code{index_maps},
#'   \code{X_event}, \code{event_re_idx}, \code{class_cols},
#'   \code{event_metadata}, \code{n_re_levels}, \code{upper_re_col},
#'   \code{middle_re_col}, \code{lower_re_col}, \code{patient_key_col},
#'   \code{admission_key_col}, \code{pathogen_col}, \code{pathogen_fitted},
#'   \code{residual_structure}, \code{estimand}, \code{prior_config_used},
#'   \code{sampler_config_used}, \code{compute_config_used}, and
#'   \code{eligibility_report}.
#'
#'   \code{diagnostics} is a one-row monitored summary. The main diagnostic
#'   fields \code{max_rhat}, \code{min_ess_bulk}, and \code{min_ess_tail}
#'   are computed over the monitored parameters retained in \code{draws}.
#'   These exclude the latent \code{z_free} nuisance parameters used for
#'   probit data augmentation. \code{converged_structural}/\code{converged_full}
#'   incorporate tree-depth saturation as well as R-hat/ESS/divergences/E-BFMI
#'   (a run that repeatedly saturates \code{max_treedepth} is not reported as
#'   converged). \code{diagnostic_status} is a canonical multi-level status --
#'   one of \code{"pass"}, \code{"warning_rhat"}, \code{"warning_low_bulk_ess"},
#'   \code{"warning_low_tail_ess"}, \code{"warning_treedepth"},
#'   \code{"fail_rhat"}, \code{"fail_energy"}, or \code{"fail_divergent"} (most
#'   severe first) -- intended as the single authoritative status field for
#'   callers, replacing ad hoc re-derivations downstream.
#'
#'   Full-fit diagnostic fields, including \code{max_rhat_all},
#'   \code{min_ess_bulk_all}, and \code{min_ess_tail_all}, are reported
#'   separately for visibility and include all Stan parameters, including
#'   \code{z_free}.
#'
#'   \code{diagnostics_detail} contains \code{monitored_parameters},
#'   \code{all_parameters} (per-parameter Rhat/ESS with \code{parameter_group}),
#'   \code{grouped} (per-\code{parameter_group} Rhat/ESS quantiles and
#'   percentage-over-threshold, not only the single worst parameter), and
#'   \code{chains} (per-chain \code{n_sampling}, \code{n_divergent},
#'   \code{n_treedepth_sat}, \code{ebfmi}, \code{mean_accept_stat},
#'   \code{mean_stepsize}, \code{mean_treedepth}, \code{max_treedepth}) -- the
#'   canonical diagnostic tables; callers should not recompute these
#'   independently from \code{fit$fit$sampler_diagnostics()}.
#' @export
fit_bayesian_multivariate_probit <- function(
  event_class_data,
  class_cols,
  fixed_effects,
  random_effects = list(),
  profile_group_col = NULL,
  pathogen = NULL,
  pathogen_col = "pathogen",
  event_id_col = "event_id",
  eligible_pairs = NULL,
  outcome_col = NULL,
  reserve_drug_cols = NULL,
  panel_eligibility = list(),
  residual_structure = c("identity", "correlated"),
  estimand = "observed_stewardship_event_mix",
  prior_config = list(),
  sampler_config = list(),
  compute = list(),
  show_messages = TRUE,
  save_full_latent_diagnostics = FALSE,
  ...
) {
  residual_structure <- match.arg(residual_structure)

  # -- Package checks ---------------------------------------------------------
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      paste0(
        "Package 'cmdstanr' is required for Pathway 2.\n",
        "  install.packages('cmdstanr', ",
        "repos = c('https://mc-stan.org/r-packages/', getOption('repos')))"
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required (installed automatically with cmdstanr).",
      call. = FALSE
    )
  }

  compute_cfg <- validate_compute_backend(compute)

  # -- Data checks ------------------------------------------------------------
  if (!is.data.frame(event_class_data) || nrow(event_class_data) == 0L) {
    stop("`event_class_data` must be a non-empty data frame.", call. = FALSE)
  }

  # -- Validate class_cols ----------------------------------------------------
  if (missing(class_cols) || length(class_cols) == 0L) {
    stop("`class_cols` is required.", call. = FALSE)
  }
  miss_cls <- setdiff(class_cols, names(event_class_data))
  if (length(miss_cls) > 0L) {
    stop(sprintf(
      "class_cols not found in event_class_data: %s",
      paste(miss_cls, collapse = ", ")
    ), call. = FALSE)
  }

  # -- Validate fixed_effects -------------------------------------------------
  if (missing(fixed_effects) || is.null(fixed_effects) || length(fixed_effects) == 0L) {
    stop("`fixed_effects` is required (no default).", call. = FALSE)
  }
  miss_fe <- setdiff(fixed_effects, names(event_class_data))
  if (length(miss_fe) > 0L) {
    stop(sprintf(
      "fixed_effects column(s) not found: %s",
      paste(miss_fe, collapse = ", ")
    ), call. = FALSE)
  }

  # -- Validate random_effects and profile grouping ----------------------------
  # Accepts EITHER the legacy character vector (c("center_name", "readmission_id"))
  # or the generic list-of-blocks spec (list(list(name=, group_col=, terms=), ...));
  # see .normalize_random_effects_spec(). Arbitrarily many blocks are allowed now
  # (Stage 1 generic random-effect architecture) -- the old "1, 2, or 3 levels"
  # ceiling was an artefact of the six hardcoded Stan variants, not a real
  # statistical limit. Composite-key safety (e.g. ensuring a patient/admission ID
  # is scoped within the right parent level) is the CALLER's responsibility --
  # prepare_random_effects() below detects nested-vs-crossed from the actual data
  # rather than assuming it, so an accidentally-non-unique group_col shows up as
  # "crossed_with_previous" rather than silently merging distinct groups.
  if (is.null(random_effects)) random_effects <- list()
  .re_blocks_early <- .normalize_random_effects_spec(random_effects)
  .re_group_cols_early <- if (length(.re_blocks_early) > 0L) {
    vapply(.re_blocks_early, function(b) b$group_col, character(1L))
  } else character(0)
  miss_re <- setdiff(.re_group_cols_early, names(event_class_data))
  if (length(miss_re) > 0L) {
    stop(sprintf(
      "random_effects group_col(s) not found: %s",
      paste(miss_re, collapse = ", ")
    ), call. = FALSE)
  }

  upper_re_col <- profile_group_col %||%
    if (length(.re_group_cols_early) > 0L) .re_group_cols_early[1L] else NULL
  if (is.null(upper_re_col) || !is.character(upper_re_col) || length(upper_re_col) != 1L ||
      !nzchar(upper_re_col)) {
    stop("`profile_group_col` is required when `random_effects` is empty.", call. = FALSE)
  }
  if (!upper_re_col %in% names(event_class_data)) {
    stop(sprintf("profile_group_col '%s' not found in event_class_data.", upper_re_col), call. = FALSE)
  }

  # -- Validate pathogen_col and optionally filter to one pathogen -----------
  if (!pathogen_col %in% names(event_class_data)) {
    stop(sprintf("pathogen_col '%s' not found in event_class_data.", pathogen_col),
      call. = FALSE
    )
  }

  event_data <- event_class_data
  pathogen_fitted <- NULL
  if (!is.null(pathogen)) {
    if (!pathogen %in% event_data[[pathogen_col]]) {
      stop(sprintf("pathogen '%s' not found in column '%s'.", pathogen, pathogen_col),
        call. = FALSE
      )
    }
    event_data <- event_data[event_data[[pathogen_col]] == pathogen, , drop = FALSE]
    pathogen_fitted <- pathogen
    message(sprintf(
      "[fit_bayesian_multivariate_probit] Filtered to pathogen '%s': %d events.",
      pathogen, nrow(event_data)
    ))
  }

  # -- Enforce single-pathogen fit (mandatory) --------------------------------
  n_pathogens <- dplyr::n_distinct(event_data[[pathogen_col]])
  if (n_pathogens != 1L) {
    stop(sprintf(
      paste0(
        "Bayesian multivariate-probit fitting requires exactly one pathogen ",
        "(found %d). Supply `pathogen = '<name>'` or pre-filter event_class_data."
      ),
      n_pathogens
    ), call. = FALSE)
  }

  # -- Resolve prior_config ---------------------------------------------------
  pc <- list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0)
  for (nm in names(prior_config)) pc[[nm]] <- prior_config[[nm]]
  # Legacy lkj_eta deliberately resolves to both controls, reproducing the
  # historical prior exactly. Explicit controls may subsequently differ.
  pc$lkj_eta_residual <- .null_default(pc$lkj_eta_residual, pc$lkj_eta)
  pc$lkj_eta_random_effect <- .null_default(pc$lkj_eta_random_effect, pc$lkj_eta)
  if (!is.numeric(pc$beta_sd) || pc$beta_sd <= 0) {
    stop("`prior_config$beta_sd` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(pc$tau_sd) || pc$tau_sd <= 0) {
    stop("`prior_config$tau_sd` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(pc$lkj_eta_residual) || pc$lkj_eta_residual < 1 ||
      !is.numeric(pc$lkj_eta_random_effect) || pc$lkj_eta_random_effect < 1) {
    stop("Resolved residual and random-effect LKJ eta values must be >= 1.", call. = FALSE)
  }

  # -- Resolve sampler_config -------------------------------------------------
  sc_defaults <- list(
    chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
    seed = 123L, adapt_delta = NULL, max_treedepth = NULL,
    parallel_chains = NULL, max_param_count = NULL,
    save_warmup = NULL, metric = NULL, init = NULL
  )
  for (nm in names(sampler_config)) sc_defaults[[nm]] <- sampler_config[[nm]]
  sc <- sc_defaults

  message(sprintf(
    "[fit_bayesian_multivariate_probit] Priors: beta~N(0,%.2g) | tau~HN(0,%.2g) | LKJ residual(%.2g) | LKJ RE(%.2g)",
    pc$beta_sd, pc$tau_sd, pc$lkj_eta_residual, pc$lkj_eta_random_effect
  ))

  # -- Remove reserve drug columns --------------------------------------------
  if (!is.null(reserve_drug_cols)) {
    class_cols <- setdiff(class_cols, reserve_drug_cols)
    message(sprintf(
      "[fit_bayesian_multivariate_probit] %d reserve drug column(s) excluded.",
      length(reserve_drug_cols)
    ))
  }
  if (length(class_cols) < 2L) {
    stop("At least 2 class columns are required for the multivariate probit.", call. = FALSE)
  }

  # -- Filter to eligible hospital-pathogen pairs -----------------------------
  if (!is.null(eligible_pairs)) {
    ep_req <- c(upper_re_col, pathogen_col)
    miss_ep <- setdiff(ep_req, names(eligible_pairs))
    if (length(miss_ep) > 0L) {
      stop(sprintf("`eligible_pairs` must have columns: %s", paste(ep_req, collapse = ", ")),
        call. = FALSE
      )
    }
    event_data <- dplyr::semi_join(event_data, eligible_pairs, by = ep_req)
    if (nrow(event_data) == 0L) {
      stop("No events remain after filtering to eligible_pairs.", call. = FALSE)
    }
    message(sprintf(
      "[fit_bayesian_multivariate_probit] %d events after eligible_pairs filter.",
      nrow(event_data)
    ))
  }

  # -- Drop all-NA event rows (events with no observed AST in any class_col) --
  has_any_obs <- rowSums(!is.na(event_data[, class_cols, drop = FALSE])) > 0L
  n_all_na <- sum(!has_any_obs)
  if (n_all_na > 0L) {
    warning(sprintf(
      paste0(
        "%d event(s) have no observed AST values across all selected class ",
        "columns and will be removed before model construction. Ensure this ",
        "is expected (e.g. events with only reserve-drug results)."
      ),
      n_all_na
    ), call. = FALSE)
    event_data <- event_data[has_any_obs, , drop = FALSE]
    if (nrow(event_data) == 0L) {
      stop("No events remain after removing all-NA event rows.", call. = FALSE)
    }
    message(sprintf(
      "[fit_bayesian_multivariate_probit] %d events remain after all-NA filter.",
      nrow(event_data)
    ))
  }

  # -- Panel eligibility report (warn, do not filter) -------------------------
  pe <- list(
    min_tested = 30L, min_resistant = 5L, min_susceptible = 5L,
    min_pairwise_cotested = 30L, min_pairwise_cell = 1L
  )
  for (nm in names(panel_eligibility)) pe[[nm]] <- panel_eligibility[[nm]]

  eligibility_report <- .validate_panel_eligibility(
    event_data, class_cols, upper_re_col, pathogen_col,
    min_tested = pe$min_tested, min_resistant = pe$min_resistant,
    min_susceptible = pe$min_susceptible
  )

  # Pairwise co-testing overlap check. min_pairwise_cell requires nonzero
  # (or user-raised) observations in every one of the 2x2 RR/RS/SR/SS cells,
  # not just an adequate total n_cotested -- a pair can clear the total
  # threshold while showing no variation at all, which carries no
  # information about residual cross-class correlation.
  elig_action <- .null_default(panel_eligibility$eligibility_action, "warn")
  co_test_report <- .check_pairwise_cotesting(
    event_data, class_cols, upper_re_col,
    min_pairwise_cotested = pe$min_pairwise_cotested,
    min_pairwise_cell = pe$min_pairwise_cell
  )
  n_low_cotested <- sum(!co_test_report$sufficient, na.rm = TRUE)
  if (n_low_cotested > 0L) {
    if (identical(residual_structure, "correlated")) {
      msg <- sprintf(
        paste0(
          "%d hospital x class-pair combination(s) have fewer than %d co-tested events. ",
          "Omega_{jk} for these pairs may be informed mainly by the LKJ prior. ",
          "See $eligibility_report$pairwise. Consider identity residuals or dropping classes with low overlap."
        ),
        n_low_cotested, pe$min_pairwise_cotested
      )
    } else {
      msg <- sprintf(
        paste0(
          "%d hospital x class-pair combination(s) have fewer than %d co-tested events. ",
          "Residual Omega is not estimated because residual_structure='identity'. ",
          "Joint profile estimates for weakly co-tested class pairs may be less data-supported. ",
          "See $eligibility_report$pairwise."
        ),
        n_low_cotested, pe$min_pairwise_cotested
      )
    }
    if (elig_action == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  n_ineligible <- sum(!eligibility_report$eligible)
  if (n_ineligible > 0L) {
    warning(
      sprintf(
        paste0(
          "%d hospital x class cell(s) below marginal eligibility thresholds ",
          "(min_tested=%d, min_resistant=%d, min_susceptible=%d). ",
          "Results for these cells may be prior-dominated. See $eligibility_report."
        ),
        n_ineligible, pe$min_tested, pe$min_resistant, pe$min_susceptible
      ),
      call. = FALSE
    )
  }

  # -- Prepare the generic random-effect representation -----------------------
  # Builds flattened level indices for an arbitrary number of blocks (Stage 1:
  # random-intercept-only). Composite-key scoping (e.g. ensuring a patient ID
  # doesn't collide across hospitals) is the CALLER's responsibility -- passing
  # an already-scoped column (or verifying global uniqueness beforehand, as
  # anumaan-analysis does) -- prepare_random_effects() detects the resulting
  # nested/crossed relationship from the data rather than assuming it.
  re_prep <- prepare_random_effects(event_data, random_effects)
  n_re_levels <- re_prep$R # kept as a field name for backward-compat call sites

  # -- Validate event identifier (real uniqueness check, not seq_len trick) ---
  if (!event_id_col %in% names(event_data)) {
    # No external event ID: assign a warning-level fallback.
    warning(sprintf(
      paste0(
        "event_id_col '%s' not found. Using row position as event key. ",
        "Ensure event_class_data has one row per organism-event -- no duplicates."
      ),
      event_id_col
    ), call. = FALSE)
    event_data[[event_id_col]] <- seq_len(nrow(event_data))
  } else {
    if (anyDuplicated(event_data[[event_id_col]])) {
      stop(sprintf(
        paste0(
          "event_id_col '%s' contains duplicate values. ",
          "Each row must correspond to exactly one organism-event."
        ),
        event_id_col
      ), call. = FALSE)
    }
  }

  # -- Canonical event index (assigned BEFORE pivoting to long) ---------------
  event_data$.event_idx <- seq_len(nrow(event_data))
  N_events <- nrow(event_data)

  if (N_events > 1000L) {
    warning(sprintf(
      paste0(
        "N_events = %d with D = %d classes adds %d latent parameters (z_free). ",
        "Sampling may be slow. Use parallel_chains in sampler_config."
      ),
      N_events, length(class_cols), N_events * length(class_cols)
    ), call. = FALSE)
  }

  # -- Columns to carry into long format ---------------------------------------
  meta_carry <- unique(c(
    ".event_idx", re_prep$group_cols, pathogen_col,
    if (!is.null(outcome_col) && outcome_col %in% names(event_data)) outcome_col,
    fixed_effects
  ))
  meta_carry <- meta_carry[!is.null(meta_carry) & meta_carry %in% names(event_data)]

  # -- Pivot to long format ---------------------------------------------------
  data_long <- event_data %>%
    dplyr::select(dplyr::all_of(c(meta_carry, class_cols))) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(class_cols),
      names_to  = "antibiotic_class",
      values_to = "resistance_binary"
    ) %>%
    dplyr::filter(!is.na(.data$resistance_binary)) %>%
    dplyr::mutate(resistance_binary = as.integer(.data$resistance_binary))

  if (nrow(data_long) == 0L) {
    stop("No observed (event x class) pairs after removing NAs.", call. = FALSE)
  }

  # -- Fail on fixed-effect missingness (do NOT impute; stop with clear msg) --
  # Checking after the pivot would count one missing value per observed class, not per event.
  fe_df_check <- event_data[, fixed_effects, drop = FALSE]
  na_cols <- vapply(fixed_effects, function(cc) sum(is.na(fe_df_check[[cc]])), integer(1L))
  if (any(na_cols > 0L)) {
    bad <- paste(sprintf("'%s' (%d NA)", names(na_cols)[na_cols > 0L], na_cols[na_cols > 0L]),
      collapse = ", "
    )
    stop(sprintf(
      paste0(
        "Fixed-effect column(s) contain NA values: %s. ",
        "Resolve missingness before calling fit_bayesian_multivariate_probit(). ",
        "Options: complete-case filter, median + indicator, explicit 'Unknown' level, ",
        "or multiple imputation -- but the decision belongs in the analysis config."
      ),
      bad
    ), call. = FALSE)
  }
  for (cc in fixed_effects) {
    if (is.character(fe_df_check[[cc]])) fe_df_check[[cc]] <- factor(fe_df_check[[cc]])
  }

  # -- Build integer index maps -----------------------------------------------
  # Per-block level indices were already computed once, generically, by
  # prepare_random_effects() above (re_prep$group_index / $flat_group_index,
  # one row per event in event_data's original row order == .event_idx). No
  # need to re-derive per-block indices through the long-format pivot.
  class_levels <- class_cols

  for (r in seq_len(re_prep$R)) {
    n_params_block <- re_prep$n_levels[r] * length(class_levels)
    if (n_params_block > 10000L) {
      warning(
        sprintf(
          "Large random-effect block '%s': %d levels x %d classes = %d parameters.",
          re_prep$block_names[r], re_prep$n_levels[r], length(class_levels), n_params_block
        ),
        call. = FALSE
      )
    }
  }

  data_long <- data_long %>%
    dplyr::mutate(
      d_idx  = match(.data$antibiotic_class, class_levels),
      ev_idx = as.integer(.data$.event_idx)
    )

  # -- Validate ev_idx (must span 1..N_events without gaps) -------------------
  stopifnot(all(data_long$ev_idx >= 1L & data_long$ev_idx <= N_events))
  stopifnot(!anyDuplicated(event_data$.event_idx))

  # event_data is already in .event_idx order; the long table would duplicate design rows by class.
  X_event_mat <- stats::model.matrix(~., data = fe_df_check)

  D <- length(class_levels)
  K <- ncol(X_event_mat)
  re_idx_mat <- re_prep$flat_group_index # N_events x R, already event-ordered

  # Sign (+1 resistant/-1 susceptible) and linear index per observed cell, for Stan's z_aug/Jacobian.
  sign_obs <- as.integer(2L * data_long$resistance_binary - 1L)
  obs_idx <- data_long$ev_idx + (data_long$d_idx - 1L) * N_events

  re_label <- if (re_prep$R == 0L) "none (fixed-only)" else {
    paste(sprintf("%s(%d)", re_prep$block_names, re_prep$n_levels), collapse = "+")
  }
  message(sprintf(
    "[fit_bayesian_multivariate_probit] %d obs | %d events | D=%d | RE blocks=%s",
    nrow(data_long), N_events, D, re_label
  ))

  # -- Build Stan data list -----------------------------------------------------
  stan_data <- list(
    N               = nrow(data_long), # observation pairs (for Jacobian)
    N_events        = N_events,
    D               = D,
    K               = K,
    X_event         = unname(X_event_mat), # N_events x K
    ev_idx          = as.integer(data_long$ev_idx),
    d_idx           = as.integer(data_long$d_idx),
    sign_obs        = as.integer(sign_obs),
    obs_idx         = as.integer(obs_idx),
    prior_beta_sd   = as.numeric(pc$beta_sd)
  )
  if (re_prep$R > 0L) {
    stan_data <- c(stan_data, list(
      R               = re_prep$R,
      total_re_levels = re_prep$total_re_levels,
      n_levels        = as.array(as.integer(re_prep$n_levels)),
      level_start     = as.array(as.integer(re_prep$level_start)),
      re_idx          = matrix(as.integer(re_idx_mat), nrow = N_events),
      prior_tau_sd    = as.numeric(pc$tau_sd),
      lkj_eta_random_effect = as.numeric(pc$lkj_eta_random_effect)
    ))
  }
  if (identical(residual_structure, "correlated")) {
    stan_data$lkj_eta_residual <- as.numeric(pc$lkj_eta_residual)
  }

  # -- Parameter-count preflight -----------------------------------------------
  n_z_free <- N_events * D
  n_re_raw <- re_prep$total_re_levels * D
  n_beta <- K * D
  n_tau <- D * re_prep$R
  n_L_corr <- (D * (D - 1L) %/% 2L) * re_prep$R
  n_L_omega <- if (residual_structure == "correlated") D * (D - 1L) %/% 2L else 0L
  n_param_approx <- n_z_free + n_re_raw + n_beta + n_tau + n_L_corr + n_L_omega
  message(sprintf(
    paste0(
      "[fit_bayesian_multivariate_probit] Approximate parameter count:\n",
      "  N_events=%d x D=%d z_free                    = %d\n",
      "  %s\n",
      "  total_re_levels=%d x D=%d random-effect raw   = %d\n",
      "  beta(%dx%d) + tau(%d) + L_corr(%d) + L_Omega(%d)\n",
      "  Total approx: %d"
    ),
    N_events, D, n_z_free,
    paste(sprintf(
      "  block '%s': %d levels x %d classes = %d params",
      re_prep$block_names, re_prep$n_levels, D, re_prep$n_levels * D
    ), collapse = "\n"),
    re_prep$total_re_levels, D, n_re_raw,
    K, D, n_tau, n_L_corr, n_L_omega,
    n_param_approx
  ))
  if (!is.null(sc$max_param_count) && n_param_approx > as.integer(sc$max_param_count)) {
    stop(sprintf(
      paste0(
        "Approximate parameter count (%d) exceeds sampler_config$max_param_count (%d). ",
        "Reduce N_events, D, or number of RE groups, or increase the threshold."
      ),
      n_param_approx, as.integer(sc$max_param_count)
    ), call. = FALSE)
  }

  stan_code <- if (re_prep$R == 0L && identical(residual_structure, "correlated")) {
    .amr_probit_stan_fixed_correlated()
  } else if (re_prep$R == 0L) {
    .amr_probit_stan_fixed_identity()
  } else if (identical(residual_structure, "correlated")) {
    .amr_probit_stan_generic_correlated()
  } else {
    .amr_probit_stan_generic_identity()
  }
  message(sprintf(
    "[fit_bayesian_multivariate_probit] Compiling %s Stan model (backend=%s)...",
    if (re_prep$R == 0L) "fixed-only" else sprintf("generic %d-block", re_prep$R), compute_cfg$backend
  ))
  compile_result <- tryCatch(
    .amr_compile_stan_backend(
      stan_code = stan_code,
      residual_structure = residual_structure,
      compute_config = compute_cfg,
      compile = TRUE,
      quiet = !isTRUE(show_messages)
    ),
    error = function(e) {
      if (identical(compute_cfg$backend, "opencl") && isTRUE(compute_cfg$allow_cpu_fallback)) {
        warning(sprintf("OpenCL compile failed; retrying on CPU. Reason: %s", conditionMessage(e)),
          call. = FALSE)
        .amr_compile_stan_backend(
          stan_code = stan_code,
          residual_structure = residual_structure,
          compute_config = utils::modifyList(compute_cfg, list(
            backend = "cpu",
            opencl_platform_id = NULL,
            opencl_device_id = NULL
          )),
          compile = TRUE,
          quiet = !isTRUE(show_messages)
        )
      } else {
        stop(conditionMessage(e), call. = FALSE)
      }
    }
  )
  mod <- compile_result$model
  actual_compute_cfg <- compile_result$compute_config
  backend_fallback <- FALSE
  backend_fallback_reason <- NULL
  if (!identical(compute_cfg$backend, actual_compute_cfg$backend)) {
    backend_fallback <- TRUE
    backend_fallback_reason <- sprintf(
      "Requested %s backend but compiled %s backend.",
      compute_cfg$backend, actual_compute_cfg$backend
    )
  }

  # -- Build sample() call args -----------------------------------------------
  sample_args <- list(
    data          = stan_data,
    seed          = as.integer(sc$seed),
    chains        = as.integer(sc$chains),
    iter_warmup   = as.integer(sc$iter_warmup),
    iter_sampling = as.integer(sc$iter_sampling),
    refresh       = if (show_messages) 200L else 0L
  )
  if (!is.null(sc$adapt_delta)) sample_args$adapt_delta <- sc$adapt_delta
  if (!is.null(sc$max_treedepth)) sample_args$max_treedepth <- sc$max_treedepth
  if (!is.null(sc$parallel_chains)) sample_args$parallel_chains <- as.integer(sc$parallel_chains)
  # save_warmup/metric/init are forwarded as-is to cmdstanr::sample() -- see
  # the sampler_config roxygen entry above for what save_warmup=TRUE does and
  # does not recover (per-iteration draws/diagnostics during warmup; NOT the
  # full evolving mass matrix at every adaptation window).
  if (!is.null(sc$save_warmup)) sample_args$save_warmup <- sc$save_warmup
  if (!is.null(sc$metric)) sample_args$metric <- sc$metric
  if (!is.null(sc$init)) sample_args$init <- sc$init
  extra_args <- list(...)
  sample_args <- c(sample_args, extra_args)

  # -- Sample -----------------------------------------------------------------
  message(sprintf(
    "[fit_bayesian_multivariate_probit] Sampling: %d chains x (%d warmup + %d sampling)...",
    sc$chains, sc$iter_warmup, sc$iter_sampling
  ))
  fit <- tryCatch(
    .amr_sample_with_backend(mod, sample_args, actual_compute_cfg),
    error = function(e) {
      if (identical(actual_compute_cfg$backend, "opencl") && isTRUE(actual_compute_cfg$allow_cpu_fallback)) {
        warning(sprintf("OpenCL sampling failed; retrying on CPU. Reason: %s", conditionMessage(e)),
          call. = FALSE)
        backend_fallback <<- TRUE
        backend_fallback_reason <<- sprintf("OpenCL sampling failed: %s", conditionMessage(e))
        actual_compute_cfg <<- utils::modifyList(actual_compute_cfg, list(
          backend = "cpu",
          opencl_platform_id = NULL,
          opencl_device_id = NULL
        ))
        cpu_compile <- .amr_compile_stan_backend(
          stan_code = stan_code,
          residual_structure = residual_structure,
          compute_config = actual_compute_cfg,
          compile = TRUE,
          quiet = !isTRUE(show_messages)
        )
        mod <<- cpu_compile$model
        actual_compute_cfg <<- cpu_compile$compute_config
        .amr_sample_with_backend(mod, sample_args, actual_compute_cfg)
      } else {
        stop(conditionMessage(e), call. = FALSE)
      }
    }
  )

  # -- Extract draws (exclude z_free; it is large and not needed downstream) --
  keep_vars <- c(
    "beta", if (re_prep$R > 0L) c("re_effect", "tau_re", "R_block"),
    if (residual_structure == "correlated") c("L_Omega", "Omega"),
    "lp__"
  )
  keep_vars <- unique(keep_vars[!is.null(keep_vars)])
  draws_arr <- tryCatch(
    fit$draws(variables = keep_vars, format = "draws_array"),
    error = function(e) fit$draws(format = "draws_array")
  )

  # Canonical parameter grouping, generic over an arbitrary number of
  # declared RE blocks: re_effect[d, level] and R_block[r][i,j] are grouped
  # by the DECLARED BLOCK NAME (via re_prep$level_start/level_end for
  # re_effect's flattened level index, and the literal block index r for
  # R_block/tau_re) rather than a hardcoded hospital/patient/admission regex.
  .probit_parameter_group <- function(parameter) {
    out <- rep("other", length(parameter))
    out[grepl("^beta(\\[|$)", parameter)] <- "beta"
    out[grepl("^(Omega|L_Omega)(\\[|$)", parameter)] <- "Omega"
    out[grepl("^z_free(\\[|$)", parameter)] <- "z_free"
    out[parameter == "lp__"] <- "lp__"

    re_effect_hits <- grep("^re_effect\\[", parameter)
    if (length(re_effect_hits) > 0L) {
      level_num <- suppressWarnings(as.integer(
        sub("^re_effect\\[\\d+,(\\d+)\\]$", "\\1", parameter[re_effect_hits])
      ))
      block_of_level <- vapply(level_num, function(lv) {
        if (is.na(lv)) {
          return(NA_integer_)
        }
        hit <- which(lv >= re_prep$level_start & lv <= re_prep$level_end)
        if (length(hit) == 0L) NA_integer_ else hit[1L]
      }, integer(1L))
      out[re_effect_hits] <- ifelse(is.na(block_of_level), "re_effect",
        re_prep$block_names[block_of_level]
      )
    }

    for (r in seq_len(re_prep$R)) {
      pat <- sprintf("^(tau_re|R_block)\\[%d(,|\\])", r)
      out[grepl(pat, parameter)] <- re_prep$block_names[r]
    }
    out
  }
  .combine_diag_tables <- function(rhat_tbl, ess_tbl) {
    out <- merge(rhat_tbl, ess_tbl, by = "variable", all = TRUE)
    out$parameter_group <- .probit_parameter_group(out$variable)
    out[, c(
      "variable", "parameter_group",
      setdiff(names(out), c("variable", "parameter_group"))
    ),
    drop = FALSE
    ]
  }
  .diag_extreme <- function(tbl, value_col, direction = c("max", "min")) {
    direction <- match.arg(direction)
    if (is.null(tbl) || !value_col %in% names(tbl)) {
      return(list(value = NA_real_, parameter = NA_character_, group = NA_character_))
    }
    vals <- suppressWarnings(as.numeric(tbl[[value_col]]))
    ok <- !is.na(vals)
    if (!any(ok)) {
      return(list(value = NA_real_, parameter = NA_character_, group = NA_character_))
    }
    idx <- base::which(ok)[if (identical(direction, "max")) which.max(vals[ok]) else which.min(vals[ok])]
    list(
      value = vals[[idx]],
      parameter = tbl$variable[[idx]],
      group = tbl$parameter_group[[idx]]
    )
  }

  # -- Strengthened diagnostics -----------------------------------------------
  # Primary diagnostics use the same monitored parameter set saved in draws_arr.
  # Full diagnostics include z_free and are reported separately for visibility.
  rhat_tbl_monitored <- fit$summary(variables = keep_vars, "rhat")
  ess_tbl_monitored <- fit$summary(variables = keep_vars, "ess_bulk", "ess_tail")
  rhat_tbl_all <- fit$summary(variables = NULL, "rhat")
  ess_tbl_all <- fit$summary(variables = NULL, "ess_bulk", "ess_tail")
  monitored_diag_tbl <- .combine_diag_tables(rhat_tbl_monitored, ess_tbl_monitored)
  all_diag_tbl <- .combine_diag_tables(rhat_tbl_all, ess_tbl_all)
  samp_diag <- fit$sampler_diagnostics(format = "matrix")

  n_divergent <- sum(samp_diag[, "divergent__"], na.rm = TRUE)
  effective_max_td <- as.integer(.null_default(sc$max_treedepth, 10L))
  n_treedepth <- if ("treedepth__" %in% colnames(samp_diag)) {
    sum(samp_diag[, "treedepth__"] >= effective_max_td,
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }
  max_rhat_monitored <- .diag_extreme(monitored_diag_tbl, "rhat", "max")
  min_bulk_monitored <- .diag_extreme(monitored_diag_tbl, "ess_bulk", "min")
  min_tail_monitored <- .diag_extreme(monitored_diag_tbl, "ess_tail", "min")
  max_rhat_all <- .diag_extreme(all_diag_tbl, "rhat", "max")
  min_bulk_all <- .diag_extreme(all_diag_tbl, "ess_bulk", "min")
  min_tail_all <- .diag_extreme(all_diag_tbl, "ess_tail", "min")

  max_rhat_structural <- max_rhat_monitored$value
  min_ess_bulk_structural <- min_bulk_monitored$value
  min_ess_tail_structural <- min_tail_monitored$value

  max_rhat_full <- max_rhat_all$value
  min_ess_bulk_full <- min_bulk_all$value
  min_ess_tail_full <- min_tail_all$value

  # -- Canonical per-chain diagnostics table -----------------------------------
  # This is the single source of chain-level diagnostics (chain id, sampling
  # count, divergences, tree-depth saturation, E-BFMI, and HMC tuning summary).
  # anumaan-analysis previously recomputed an equivalent table independently
  # from `fit$fit$sampler_diagnostics()` (including re-deriving E-BFMI with the
  # same formula); that duplicate should be replaced with this table.
  n_chains_used <- as.integer(sc$chains)
  chain_id <- if ("chain__" %in% colnames(samp_diag)) {
    samp_diag[, "chain__"]
  } else {
    n_per_chain <- nrow(samp_diag) %/% n_chains_used
    rep(seq_len(n_chains_used), each = n_per_chain)[seq_len(nrow(samp_diag))]
  }
  chain_ids <- sort(unique(chain_id))

  .safe_col <- function(name) if (name %in% colnames(samp_diag)) samp_diag[, name] else rep(NA_real_, nrow(samp_diag))
  energy_col <- .safe_col("energy__")
  divergent_col <- .safe_col("divergent__")
  treedepth_col <- .safe_col("treedepth__")
  accept_col <- .safe_col("accept_stat__")
  stepsize_col <- .safe_col("stepsize__")

  chain_diag_tbl <- dplyr::bind_rows(lapply(chain_ids, function(ch) {
    idx <- chain_id == ch
    e_ch <- suppressWarnings(as.numeric(energy_col[idx]))
    ebfmi_ch <- tryCatch(
      {
        if (sum(!is.na(e_ch)) < 2L) {
          NA_real_
        } else {
          var_e <- stats::var(e_ch, na.rm = TRUE)
          if (is.na(var_e) || var_e < .Machine$double.eps) {
            NA_real_
          } else {
            mean(diff(e_ch)^2, na.rm = TRUE) / var_e
          }
        }
      },
      error = function(e) NA_real_
    )
    tibble::tibble(
      chain               = ch,
      n_sampling          = sum(idx),
      n_divergent         = sum(suppressWarnings(as.numeric(divergent_col[idx])) > 0, na.rm = TRUE),
      n_treedepth_sat     = sum(suppressWarnings(as.numeric(treedepth_col[idx])) >= effective_max_td, na.rm = TRUE),
      ebfmi               = ebfmi_ch,
      mean_accept_stat    = mean(suppressWarnings(as.numeric(accept_col[idx])), na.rm = TRUE),
      mean_stepsize       = mean(suppressWarnings(as.numeric(stepsize_col[idx])), na.rm = TRUE),
      mean_treedepth      = mean(suppressWarnings(as.numeric(treedepth_col[idx])), na.rm = TRUE),
      max_treedepth       = suppressWarnings(max(as.numeric(treedepth_col[idx]), na.rm = TRUE))
    )
  }))

  ebfmi_vals <- chain_diag_tbl$ebfmi
  ebfmi_min <- min(ebfmi_vals, na.rm = TRUE)
  ebfmi_mean <- mean(ebfmi_vals, na.rm = TRUE)
  if (is.infinite(ebfmi_min)) ebfmi_min <- NA_real_
  if (is.nan(ebfmi_mean)) ebfmi_mean <- NA_real_
  ebfmi_by_chain <- if (all(is.na(ebfmi_vals))) {
    NA_character_
  } else {
    paste(sprintf("chain%s=%.3f", chain_diag_tbl$chain, ebfmi_vals), collapse = "; ")
  }

  # Tree-depth saturation is included here (it was previously computed and
  # warned about below but silently excluded from the pass/fail flags, which
  # meant a run that repeatedly saturated max_treedepth could still report
  # converged_structural = TRUE).
  treedepth_ok <- isTRUE(is.na(n_treedepth) || n_treedepth == 0L)

  converged_structural <-
    isTRUE(max_rhat_structural < 1.01) &&
      isTRUE(min_ess_bulk_structural >= 100) &&
      isTRUE(is.na(min_ess_tail_structural) || min_ess_tail_structural >= 100) &&
      isTRUE(n_divergent == 0L) &&
      treedepth_ok &&
      isTRUE(is.na(ebfmi_min) || ebfmi_min >= 0.3)
  # Full-scope analogue of converged_structural (includes z_free). Reported
  # for visibility -- with tens of thousands of z_free parameters this will
  # often be FALSE even when converged_structural is TRUE, which is expected,
  # not a failure. converged_structural remains the recommended pass/fail
  # flag for the resistance-profile model.
  converged_full <-
    isTRUE(max_rhat_full < 1.01) &&
      isTRUE(min_ess_bulk_full >= 100) &&
      isTRUE(is.na(min_ess_tail_full) || min_ess_tail_full >= 100) &&
      isTRUE(n_divergent == 0L) &&
      treedepth_ok &&
      isTRUE(is.na(ebfmi_min) || ebfmi_min >= 0.3)
  # Narrower, more diagnostic signal than "converged_full is FALSE": TRUE only
  # when the structural parameters have converged cleanly AND the full scope
  # (i.e. the z_free latent block specifically) has not -- isolating "the
  # augmented-data nuisance variables mix less well" from "the model itself
  # has not converged" (the latter already shows up as converged_structural
  # = FALSE and needs no separate flag).
  latent_diagnostic_warning <-
    isTRUE(converged_structural) &&
      (isTRUE(max_rhat_full > 1.01) || isTRUE(min_ess_bulk_full < 100))

  # Canonical multi-level status (supersedes the separate, duplicated
  # `.probit_diagnostic_status()` previously computed in anumaan-analysis from
  # a re-derived chain-diagnostics table). Severity order: divergences and
  # poor energy exploration are the least recoverable (fail); tree-depth
  # saturation and elevated Rhat/low ESS are recoverable with more
  # warmup/sampling/adapt_delta (warning). This is computed on the STRUCTURAL
  # scope, consistent with converged_structural.
  diagnostic_status <- dplyr::case_when(
    isTRUE(n_divergent > 0L) ~ "fail_divergent",
    isTRUE(!is.na(ebfmi_min) && ebfmi_min < 0.3) ~ "fail_energy",
    isTRUE(!is.na(max_rhat_structural) && max_rhat_structural > 1.05) ~ "fail_rhat",
    isTRUE(!is.na(n_treedepth) && n_treedepth > 0L) ~ "warning_treedepth",
    isTRUE(!is.na(min_ess_bulk_structural) && min_ess_bulk_structural < 100) ~ "warning_low_bulk_ess",
    isTRUE(!is.na(min_ess_tail_structural) && min_ess_tail_structural < 100) ~ "warning_low_tail_ess",
    isTRUE(!is.na(max_rhat_structural) && max_rhat_structural > 1.01) ~ "warning_rhat",
    TRUE ~ "pass"
  )

  # -- Warnings are based on the STRUCTURAL scope; z_free stragglers are ------
  # -- reported separately below so they don't masquerade as model failure. --
  if (max_rhat_structural > 1.01) {
    warning(sprintf(
      "Convergence concern: structural max Rhat = %.3f (> 1.01). %s",
      max_rhat_structural,
      "Increase iter_warmup or adapt_delta."
    ), call. = FALSE)
  }
  if (min_ess_bulk_structural < 100L) {
    warning(sprintf(
      "Low structural bulk ESS: %.0f (< 100). Profile probabilities may be unreliable.",
      min_ess_bulk_structural
    ), call. = FALSE)
  }
  if (!is.na(min_ess_tail_structural) && min_ess_tail_structural < 100L) {
    warning(sprintf("Low structural tail ESS: %.0f (< 100).", min_ess_tail_structural), call. = FALSE)
  }
  if (!is.na(n_treedepth) && n_treedepth > 0L) {
    warning(sprintf(
      "%d iteration(s) saturated max treedepth. Increase max_treedepth.",
      n_treedepth
    ), call. = FALSE)
  }
  if (!is.na(ebfmi_min) && ebfmi_min < 0.3) {
    warning(sprintf(
      "Low minimum chain E-BFMI = %.3f (< 0.3). Possible poor geometry.",
      ebfmi_min
    ), call. = FALSE)
  }
  if (n_divergent > 0L) {
    warning(sprintf("%d divergent transition(s). Increase adapt_delta.", n_divergent),
      call. = FALSE
    )
  }
  if (latent_diagnostic_warning) {
    message(sprintf(
      paste0(
        "[fit_bayesian_multivariate_probit] NOTE: latent z_free block (%d parameters) has ",
        "diagnostic stragglers (full-scope max Rhat = %.3f, min ESS bulk = %.0f) even though ",
        "structural parameters converged cleanly. This is informational only and does NOT ",
        "affect converged_structural -- see $diagnostics$latent_diagnostic_warning."
      ),
      n_z_free, max_rhat_full, min_ess_bulk_full
    ))
  }

  diag_tbl <- tibble::tibble(
    n_chains = as.integer(sc$chains),
    iter_warmup = as.integer(sc$iter_warmup),
    iter_sampling = as.integer(sc$iter_sampling),
    n_re_levels = as.integer(re_prep$R),
    n_observed_pairs = nrow(data_long),
    n_events = N_events,
    n_classes = D,
    # n_upper_groups: backward-compatible alias for "the first declared RE
    # block's level count" (by convention the hospital-equivalent level).
    n_upper_groups = as.integer(re_prep$n_levels[1L]),
    re_block_levels = paste(sprintf("%s=%d", re_prep$block_names, re_prep$n_levels), collapse = "; "),
    n_divergent = as.integer(n_divergent),
    n_treedepth_sat = as.integer(n_treedepth),
    ebfmi_min = round(ebfmi_min, 4L),
    ebfmi_mean = round(ebfmi_mean, 4L),
    ebfmi_by_chain = ebfmi_by_chain,
    max_rhat_structural = round(max_rhat_structural, 4L),
    min_ess_bulk_structural = round(min_ess_bulk_structural, 1L),
    min_ess_tail_structural = round(min_ess_tail_structural, 1L),
    max_rhat_full = round(max_rhat_full, 4L),
    min_ess_bulk_full = round(min_ess_bulk_full, 1L),
    min_ess_tail_full = round(min_ess_tail_full, 1L),
    # Full scope (incl. z_free) is summarised by these scalars regardless of
    # save_full_latent_diagnostics; the per-parameter table itself (tens of
    # thousands of rows for realistic N_events x D) is only ever returned in
    # diagnostics_detail$all_parameters when that flag is TRUE.
    n_params_full = nrow(all_diag_tbl),
    pct_rhat_full_above_1_01 = round(100 * mean(suppressWarnings(as.numeric(all_diag_tbl$rhat)) > 1.01, na.rm = TRUE), 2L),
    pct_ess_bulk_full_below_100 = round(100 * mean(suppressWarnings(as.numeric(all_diag_tbl$ess_bulk)) < 100, na.rm = TRUE), 2L),
    parameter_with_max_rhat = max_rhat_monitored$parameter,
    parameter_group_with_max_rhat = max_rhat_monitored$group,
    parameter_with_min_ess_bulk = min_bulk_monitored$parameter,
    parameter_group_with_min_ess_bulk = min_bulk_monitored$group,
    parameter_with_min_ess_tail = min_tail_monitored$parameter,
    parameter_group_with_min_ess_tail = min_tail_monitored$group,
    converged_structural = converged_structural,
    converged_full = converged_full,
    latent_diagnostic_warning = latent_diagnostic_warning,
    diagnostic_status = diagnostic_status,
    # backward-compatible aliases (structural diagnostics)
    max_rhat = round(max_rhat_structural, 4L),
    min_ess_bulk = round(min_ess_bulk_structural, 1L),
    min_ess_tail = round(min_ess_tail_structural, 1L)
  )

  # -- Grouped (per parameter_group) diagnostic summary ------------------------
  # Reports quantiles/percentages, not only the single worst parameter, so a
  # single outlier does not stand in for the whole group's convergence.
  # Consolidates what anumaan-analysis previously recomputed independently in
  # `.probit_group_diagnostics()`.
  grouped_diag_tbl <- if (nrow(monitored_diag_tbl) == 0L) {
    tibble::tibble()
  } else {
    monitored_diag_tbl %>%
      dplyr::mutate(
        rhat = suppressWarnings(as.numeric(.data$rhat)),
        ess_bulk = suppressWarnings(as.numeric(.data$ess_bulk)),
        ess_tail = suppressWarnings(as.numeric(.data$ess_tail))
      ) %>%
      dplyr::group_by(.data$parameter_group) %>%
      dplyr::summarise(
        n_parameters = dplyr::n(),
        median_rhat = round(stats::median(.data$rhat, na.rm = TRUE), 4L),
        p95_rhat = round(stats::quantile(.data$rhat, 0.95, na.rm = TRUE), 4L),
        max_rhat = round(max(.data$rhat, na.rm = TRUE), 4L),
        pct_rhat_above_1_01 = round(100 * mean(.data$rhat > 1.01, na.rm = TRUE), 1L),
        median_ess_bulk = round(stats::median(.data$ess_bulk, na.rm = TRUE), 1L),
        p5_ess_bulk = round(stats::quantile(.data$ess_bulk, 0.05, na.rm = TRUE), 1L),
        min_ess_bulk = round(min(.data$ess_bulk, na.rm = TRUE), 1L),
        pct_ess_bulk_below_100 = round(100 * mean(.data$ess_bulk < 100, na.rm = TRUE), 1L),
        min_ess_tail = round(min(.data$ess_tail, na.rm = TRUE), 1L),
        .groups = "drop"
      )
  }

  message(sprintf(
    paste0(
      "[fit_bayesian_multivariate_probit] Done. structural max_Rhat=%.3f (full incl. z_free=%.3f) | ",
      "structural min_ESS_bulk=%.0f (full=%.0f) | divergent=%d | converged_structural=%s"
    ),
    max_rhat_structural, max_rhat_full,
    min_ess_bulk_structural, min_ess_bulk_full,
    n_divergent, converged_structural
  ))

  list(
    draws = draws_arr,
    diagnostics = diag_tbl,
    diagnostics_detail = list(
      monitored_parameters = tibble::as_tibble(monitored_diag_tbl),
      # Full per-parameter table (includes z_free) is large -- tens of
      # thousands of rows for realistic N_events x D -- and is only retained
      # when explicitly requested. diag_tbl$max_rhat_full/min_ess_bulk_full/
      # min_ess_tail_full/n_params_full/pct_rhat_full_above_1_01/
      # pct_ess_bulk_full_below_100 are always populated regardless.
      all_parameters = if (isTRUE(save_full_latent_diagnostics)) {
        tibble::as_tibble(all_diag_tbl)
      } else {
        tibble::as_tibble(all_diag_tbl)[0L, ]
      },
      grouped = grouped_diag_tbl,
      chains = chain_diag_tbl
    ),
    fit = fit,
    data_long = tibble::as_tibble(data_long),
    index_maps = list(
      class_levels      = class_levels,
      # backward-compatible alias: first declared block's level set
      upper_levels      = if (re_prep$R > 0L) re_prep$level_maps[[1L]] else NULL
    ),
    X_event = X_event_mat,
    # Generic random-effect representation (Stage 1): the self-contained
    # object from prepare_random_effects(), including flat_group_index
    # (N_events x R, one flattened level index per event per declared
    # block) -- the SINGLE generic replacement for the old event_re_idx
    # list of hardcoded h_ev/p_ev/a_ev arrays. Every downstream
    # mu-reconstruction site should sum over blocks via re_contribution()
    # against this object rather than hand-summing named arrays.
    random_effects_prep = re_prep,
    event_re_idx = re_prep$flat_group_index, # N_events x R (generic; see random_effects_prep)
    class_cols = class_cols,
    event_metadata = tibble::as_tibble(event_data),
    n_re_levels = re_prep$R,
    upper_re_col = upper_re_col,
    profile_group_col = upper_re_col,
    pathogen_col = pathogen_col,
    pathogen_fitted = pathogen_fitted,
    residual_structure = residual_structure,
    estimand = estimand,
    prior_config_used = pc,
    panel_eligibility_used = pe,
    sampler_config_used = sc,
    compute_config_used = list(
      requested_backend = compute_cfg$backend,
      actual_backend = actual_compute_cfg$backend,
      stan_opencl_enabled = identical(actual_compute_cfg$backend, "opencl"),
      opencl_platform_id = actual_compute_cfg$opencl_platform_id,
      opencl_device_id = actual_compute_cfg$opencl_device_id,
      allow_cpu_fallback = compute_cfg$allow_cpu_fallback,
      backend_fallback = backend_fallback,
      backend_fallback_reason = backend_fallback_reason,
      cmdstan_version = tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NA_character_),
      cmdstanr_version = tryCatch(as.character(utils::packageVersion("cmdstanr")), error = function(e) NA_character_),
      compile_cache_key = compile_result$cache_key,
      compiled_basename = compile_result$basename
    ),
    eligibility_report = list(
      marginal  = eligibility_report,
      pairwise  = co_test_report
    )
  )
}


#' Conditional profile probabilities under a correlated residual, via Gibbs
#'
#' For one hospital-pathogen panel at one posterior draw, samples the missing
#' latent class dimensions conditional on the observed resistance SIGNS of
#' the tested classes (we only ever observe \eqn{\mathrm{sign}(Z_{ed})}, never
#' the exact latent value), via a Gibbs sampler on the truncated multivariate
#' normal \eqn{Z_e \sim N(\mu_e, \Omega)}, vectorized across all events in the
#' panel simultaneously (the per-dimension conditional-normal regression
#' coefficients depend only on \eqn{\Omega} -- fixed across events at a given
#' draw -- so one matrix operation updates every event's Gibbs state at once).
#'
#' At every iteration, dimensions with an observed value are drawn from a
#' truncated normal restricted to the half-line matching the observed sign
#' (inverse-CDF sampling); dimensions without an observed value are drawn
#' from the unrestricted conditional normal. Because the observed dimensions
#' are truncation-constrained at every iteration, every kept iteration's
#' pattern is automatically consistent with the observed cells -- profiles
#' inconsistent with observed data are never visited, rather than assigned
#' zero probability after the fact. Counting each kept iteration's full
#' binary pattern against the enumerated profile list gives an empirical
#' conditional profile-probability vector per event that sums to exactly 1
#' by construction (every iteration falls into exactly one profile).
#'
#' \code{Dp == 1} (a single-class panel) has no joint structure to speak of
#' and is handled directly via \eqn{\Phi(\mu_{ed})}, bypassing Gibbs entirely.
#'
#' @param mu_hp Numeric matrix, n_hp x Dp. Linear predictor for this draw.
#' @param Omega_hp Numeric matrix, Dp x Dp. Residual correlation matrix for
#'   this draw (already symmetrized/regularised by the caller).
#' @param obs_panel Numeric matrix, n_hp x Dp, values in \{0, 1, NA\}.
#' @param profile_bin Integer matrix, n_profiles x Dp, 0/1 (from
#'   \code{enumerate_binary_profiles()}); row \eqn{j} must correspond to
#'   integer code \eqn{j - 1} in binary (column \eqn{d} = bit \eqn{d - 1}),
#'   which is how \code{enumerate_binary_profiles()} orders rows.
#' @param n_burnin,n_kept Integer. Gibbs iterations discarded as burn-in and
#'   iterations kept for the empirical profile-probability count,
#'   respectively, at this single posterior draw. Total Gibbs cost scales
#'   with \code{(n_burnin + n_kept) * n_posterior_draws_for_profiles}, so
#'   correlated-residual profile generation should typically use fewer
#'   posterior draws than an identity-residual fit of the same size.
#'
#' @return Numeric matrix, n_hp x n_profiles: empirical conditional profile
#'   probability per event, rows summing to 1.
#' @keywords internal
.gibbs_conditional_profile_probs <- function(mu_hp, Omega_hp, obs_panel, profile_bin,
                                             n_burnin = 10L, n_kept = 20L) {
  n_hp <- nrow(mu_hp)
  Dp <- ncol(mu_hp)
  n_profiles <- nrow(profile_bin)
  pow2 <- 2L^(seq_len(Dp) - 1L)

  if (Dp == 1L) {
    # No joint structure with a single class; identical to the analytic
    # identity-residual formula. profile_bin row 1 = code 0 = "S" (bit
    # unset), row 2 = code 1 = "R" (bit set) -- see enumerate_binary_profiles().
    p <- stats::pnorm(mu_hp[, 1])
    p_eff <- obs_panel[, 1]
    miss <- is.na(p_eff)
    p_eff[miss] <- p[miss]
    prob_mat <- matrix(0, n_hp, n_profiles)
    prob_mat[, 1] <- 1 - p_eff
    prob_mat[, 2] <- p_eff
    return(prob_mat)
  }

  known_mask <- !is.na(obs_panel)

  # Per-dimension conditional-normal regression coefficients depend only on
  # Omega (fixed across events at this draw) -- precomputed once. Near-
  # singular Omega_oo (very high |correlation| among the conditioning
  # dimensions) is handled with escalating controlled jitter; if that still
  # doesn't produce a valid (finite, positive) conditional variance, this
  # stops with an explicit, diagnosable reason rather than silently
  # returning NaN/negative-variance draws.
  beta_list <- vector("list", Dp)
  cond_sd <- numeric(Dp)
  for (d in seq_len(Dp)) {
    other <- setdiff(seq_len(Dp), d)
    Omega_oo <- Omega_hp[other, other, drop = FALSE]

    solved <- FALSE
    for (jitter in c(1e-8, 1e-6, 1e-4, 1e-2)) {
      Omega_oo_inv <- tryCatch(solve(Omega_oo + diag(jitter, Dp - 1L)), error = function(e) NULL)
      if (is.null(Omega_oo_inv)) next
      beta_d <- as.vector(Omega_hp[d, other] %*% Omega_oo_inv)
      cond_var <- Omega_hp[d, d] - sum(Omega_hp[d, other] * beta_d)
      if (is.finite(cond_var) && cond_var > 1e-9) {
        solved <- TRUE
        break
      }
    }
    if (!solved) {
      stop(sprintf(
        paste0(
          "[.gibbs_conditional_profile_probs] Omega is too near-singular to condition on ",
          "dimension %d even after jitter up to 1e-2 (conditioning set size %d). This usually ",
          "means the fitted residual correlation matrix has classes that are (near-)perfectly ",
          "correlated -- inspect fit$diagnostics_detail for this draw/hospital-pathogen pair, ",
          "or drop one of the collinear classes from the panel."
        ),
        d, Dp - 1L
      ), call. = FALSE)
    }

    beta_list[[d]] <- beta_d
    cond_sd[d] <- sqrt(cond_var)
  }

  # Initialise Z: observed dimensions nudged onto the correct side of 0 if
  # mu itself has the wrong sign; missing dimensions start at mu.
  Z <- mu_hp
  for (d in seq_len(Dp)) {
    k <- known_mask[, d]
    if (!any(k)) next
    want_pos <- obs_panel[k, d] == 1
    cur <- Z[k, d]
    wrong <- (want_pos & cur <= 0) | (!want_pos & cur >= 0)
    if (any(wrong)) {
      cur[wrong] <- ifelse(want_pos[wrong], abs(cur[wrong]) + 0.1, -abs(cur[wrong]) - 0.1)
      Z[k, d] <- cur
    }
  }

  prob_mat <- matrix(0, n_hp, n_profiles)
  n_total_iter <- n_burnin + n_kept

  for (g in seq_len(n_total_iter)) {
    for (d in seq_len(Dp)) {
      other <- setdiff(seq_len(Dp), d)
      resid <- Z[, other, drop = FALSE] - mu_hp[, other, drop = FALSE]
      mu_cond <- mu_hp[, d] + as.vector(resid %*% beta_list[[d]])
      sd_cond <- cond_sd[d]

      k <- known_mask[, d]
      z_new <- numeric(n_hp)
      if (any(!k)) z_new[!k] <- stats::rnorm(sum(!k), mu_cond[!k], sd_cond)
      if (any(k)) {
        want_pos <- obs_panel[k, d] == 1
        p0 <- stats::pnorm(0, mu_cond[k], sd_cond) # P(Z < 0 | ...)
        u <- ifelse(want_pos, stats::runif(sum(k), p0, 1), stats::runif(sum(k), 0, p0))
        u <- pmin(pmax(u, 1e-10), 1 - 1e-10) # guard qnorm(0)/qnorm(1) = +-Inf
        z_new[k] <- stats::qnorm(u, mu_cond[k], sd_cond)
      }
      Z[, d] <- z_new
    }
    if (g > n_burnin) {
      signs <- (Z > 0) * 1L
      idx <- as.vector(signs %*% pow2) + 1L
      prob_mat[cbind(seq_len(n_hp), idx)] <- prob_mat[cbind(seq_len(n_hp), idx)] + 1
    }
  }

  prob_mat / n_kept
}


# compute_event_profile_probabilities()

#' Compute Observed-Plus-Imputed Resistance Profile Probabilities
#'
#' For each event, retains every observed (tested) AST cell exactly and
#' computes a probability distribution only over the genuinely untested
#' (\code{NA}) antibiotic-class cells in that event's profile panel. Any
#' enumerated profile inconsistent with the event's observed cells receives
#' probability 0 -- this function never overwrites an observed R/S result.
#'
#' \strong{Identity residual structure} (\code{fitted_model$residual_structure
#' == "identity"}): classes are conditionally independent given the linear
#' predictor, so the missing-cell probabilities are computed analytically as
#' \eqn{P(Y_{ed}=1 \mid \theta) = \Phi(\mu_{ed})} for each missing class
#' \eqn{d}, and the profile probability over the missing dimensions is the
#' exact product of per-class Bernoulli probabilities -- no latent-variable
#' simulation is used, so there is no added latent-profile-simulation noise.
#' Rows carry \code{profile_generation_method = "conditional_analytic_identity"}.
#'
#' \strong{Correlated residual structure}: the missing latent
#' dimensions are sampled conditional on the observed resistance SIGNS of the
#' tested classes via a Gibbs sampler on the truncated multivariate normal --
#' see \code{.gibbs_conditional_profile_probs()}. This replaces the earlier
#' unconditional latent-\eqn{Z} simulation, which resampled tested cells
#' along with untested ones and was therefore not observed-plus-imputed.
#' Rows carry \code{profile_generation_method = "conditional_gibbs_correlated"}
#' and ARE eligible for downstream DALY use (subject to the same panel-support
#' and sampler-acceptability gates as identity-residual profiles -- see
#' \code{aggregate_profiles_for_daly()}). Because Gibbs is only run for
#' \code{n_gibbs_burnin + n_gibbs_kept} iterations per draw rather than to
#' full convergence, correlated-residual profiles carry additional Monte
#' Carlo error beyond identity-residual profiles at the same
#' \code{n_posterior_draws_for_profiles}; increase \code{n_gibbs_kept} (and/or
#' \code{n_gibbs_burnin}) for lower-noise estimates at higher compute cost.
#'
#' Both residual structures produce, per event per posterior draw, a profile
#' probability vector over the panel's \eqn{2^{D_p}} enumerated profiles that
#' (a) sums to exactly 1 and (b) assigns exactly 0 to any profile
#' inconsistent with that event's observed cells -- enforced by construction
#' in both cases (an exact product for identity; Gibbs never visiting an
#' inconsistent pattern for correlated), and checked with a hard
#' \code{stop()} for the identity path since it has no other source of
#' Monte Carlo noise to blur a genuine bug.
#'
#' \strong{Class panels}: the antibiotic-class panel enumerated for each
#' hospital x pathogen pair is drawn from the approved eligibility rules
#' computed at fit time (\code{fitted_model$eligibility_report}: marginal
#' n_tested/n_resistant/n_susceptible thresholds, plus pairwise co-testing
#' sufficiency for correlated-residual fits) -- not from "tested at least
#' once". See \code{.resolve_profile_class_panel()}.
#'
#' \strong{Estimand:} The posterior distribution over the observed event
#' case-mix, conditional on each event's own observed AST results. This is
#' labelled \code{"observed_stewardship_event_mix"}.
#'
#' @param fitted_model List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param n_posterior_draws_for_profiles Integer. Number of posterior draws
#'   averaged over for both the event-level profile probabilities and the
#'   draw-level aggregates. Subsampled without replacement when total draws
#'   exceed this value. Default \code{2000L}.
#' @param outcome_col Character or \code{NULL}. Patient outcome column in
#'   \code{fitted_model$event_metadata}. When \code{NULL}, all events are
#'   treated as having a known outcome and R_NF is not computed separately.
#' @param nonfatal_values Character vector. Outcome values for the non-fatal
#'   cohort (R_NF). Default covers common discharge/survival labels.
#' @param seed Integer. Random seed. Default \code{123L}.
#' @param n_gibbs_burnin,n_gibbs_kept Integer. Only used when
#'   \code{fitted_model$residual_structure == "correlated"} -- see
#'   \code{.gibbs_conditional_profile_probs()}. Defaults \code{10L}/\code{20L}.
#' @param posterior_draw_indices Optional unique indices into the fitted
#'   posterior draws. When supplied, these replace random draw subsampling.
#' @param event_indices Optional canonical event indices to include. When
#'   supplied, only observed events with these indices are processed.
#'
#' @return Named list: \code{event_profiles} (event-level posterior mean
#'   observed-plus-imputed profile probabilities, with
#'   \code{n_classes_observed}/\code{n_classes_missing} per event) and
#'   \code{aggregate_draws} (per-draw \code{R_ALL} [truly all events in the
#'   hp pair], \code{R_KNOWN_OUTCOME} [known-outcome cohort], and \code{R_NF}
#'   [nonfatal cohort], used by \code{aggregate_profiles_for_daly()} for
#'   credible intervals). Both tibbles contain all \eqn{2^D} profiles per
#'   hospital-pathogen pair and carry \code{profile_generation_method} and
#'   \code{panel_reason} columns.
#' @export
compute_event_profile_probabilities <- function(
  fitted_model,
  n_posterior_draws_for_profiles = 2000L,
  outcome_col = NULL,
  nonfatal_values = c(
    "Discharged", "Survived", "Alive",
    "discharged", "survived", "alive"
  ),
  seed = 123L,
  n_gibbs_burnin = 10L,
  n_gibbs_kept = 20L,
  posterior_draw_indices = NULL,
  event_indices = NULL
) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)
  }

  set.seed(as.integer(seed))

  draws <- fitted_model$draws
  class_cols <- fitted_model$class_cols
  idx_maps <- fitted_model$index_maps
  event_meta <- fitted_model$event_metadata
  n_re_levels <- fitted_model$n_re_levels
  upper_re_col <- fitted_model$upper_re_col
  pathogen_col <- fitted_model$pathogen_col
  residual_structure <- .null_default(fitted_model$residual_structure, "identity")
  eligibility_report <- fitted_model$eligibility_report

  D <- length(idx_maps$class_levels)
  re_prep <- fitted_model$random_effects_prep

  X_event_sim <- fitted_model$X_event
  event_re_idx <- fitted_model$event_re_idx # N_events x R flattened level index

  K <- ncol(X_event_sim)

  # -- Thin draws to n_posterior_draws_for_profiles ---------------------------
  draws_mat <- posterior::as_draws_matrix(draws)
  n_total <- nrow(draws_mat)
  if (is.null(posterior_draw_indices)) {
    S <- min(as.integer(n_posterior_draws_for_profiles), n_total)
    draw_idx <- if (S < n_total) sort(sample.int(n_total, S)) else seq_len(n_total)
  } else {
    draw_idx <- as.integer(posterior_draw_indices)
    if (length(draw_idx) < 1L || anyNA(draw_idx) || any(draw_idx < 1L | draw_idx > n_total) || anyDuplicated(draw_idx)) {
      stop("`posterior_draw_indices` must be unique valid indices into fitted_model$draws.", call. = FALSE)
    }
    S <- length(draw_idx)
  }
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  # -- Helper: extract [S, d1, d2] array from draws ---------------------------
  .arr <- function(prefix, d1, d2) {
    cols <- as.vector(outer(
      seq_len(d1), seq_len(d2),
      function(a, b) sprintf("%s[%d,%d]", prefix, a, b)
    ))
    array(draws_mat[, cols, drop = FALSE], dim = c(S, d1, d2))
  }

  beta_arr <- .arr("beta", K, D) # S x K x D
  re_eff_arr <- if (re_prep$R > 0L) {
    .arr("re_effect", D, re_prep$total_re_levels)
  } else NULL
  L_omega_arr <- if (residual_structure == "correlated") {
    .arr("L_Omega", D, D)
  } # S x D x D (lower triangular)
  else {
    NULL
  }

  # -- Canonical event table --------------------------------------------------
  stopifnot(".event_idx" %in% names(event_meta))
  stopifnot(!anyDuplicated(event_meta$.event_idx))

  # Use event_meta to identify which events appear in data_long (have observations)
  has_obs <- event_meta$.event_idx %in% fitted_model$data_long$ev_idx
  event_meta_obs <- event_meta[has_obs, , drop = FALSE]
  if (!is.null(event_indices)) {
    event_indices <- as.integer(event_indices)
    event_meta_obs <- event_meta_obs[event_meta_obs$.event_idx %in% event_indices, , drop = FALSE]
    if (nrow(event_meta_obs) == 0L) {
      stop("`event_indices` selects no events with observed AST data.", call. = FALSE)
    }
  }
  N_ev <- nrow(event_meta_obs)

  # Event-level design matrix and RE indices (from stored event_re_idx, the
  # generic N_events x R flattened level-index matrix -- see random_effects_prep).
  ev_row_idx <- event_meta_obs$.event_idx # 1-based canonical event indices
  X_event <- X_event_sim[ev_row_idx, , drop = FALSE] # N_ev x K
  flat_re_idx_obs <- event_re_idx[ev_row_idx, , drop = FALSE] # N_ev x R

  # Observed AST matrix for these events, aligned to class_cols (0/1/NA). This
  # is the ground truth that must be preserved: profile generation below only
  # ever computes probabilities for the NA (untested) cells.
  obs_ast_mat <- as.matrix(event_meta_obs[, class_cols, drop = FALSE])
  storage.mode(obs_ast_mat) <- "double"

  # -- Outcome flags per event ------------------------------------------------
  oc_col_use <- if (!is.null(outcome_col) && outcome_col %in% names(event_meta_obs)) {
    outcome_col
  } else {
    NULL
  }

  if (!is.null(oc_col_use)) {
    ov <- event_meta_obs[[oc_col_use]]
    is_known_outcome <- !is.na(ov)
    is_nonfatal <- !is.na(ov) & (ov %in% nonfatal_values)
  } else {
    is_known_outcome <- rep(TRUE, N_ev)
    is_nonfatal <- rep(TRUE, N_ev)
  }

  # -- Hospital-pathogen profile class panels ---------------------------------
  # Panels come from the approved eligibility rules (fitted_model$eligibility_
  # report), not from "tested at least once" -- see .resolve_profile_class_panel().
  hp_pairs <- unique(event_meta_obs[, c(upper_re_col, pathogen_col), drop = FALSE])
  hp_keys <- paste(hp_pairs[[upper_re_col]], hp_pairs[[pathogen_col]], sep = "||")

  hp_panel <- stats::setNames(lapply(seq_len(nrow(hp_pairs)), function(r) {
    .resolve_profile_class_panel(
      class_cols          = class_cols,
      hospital            = hp_pairs[[upper_re_col]][r],
      pathogen            = hp_pairs[[pathogen_col]][r],
      eligibility_report  = eligibility_report,
      upper_re_col        = upper_re_col,
      pathogen_col        = pathogen_col,
      residual_structure  = residual_structure
    )
  }), hp_keys)

  # Map hp_keys -> row indices in event_meta_obs (canonical, not event_meta!)
  hp_ev_idx <- stats::setNames(lapply(hp_keys, function(key) {
    parts <- strsplit(key, "||", fixed = TRUE)[[1L]]
    which(event_meta_obs[[upper_re_col]] == parts[1L] &
      event_meta_obs[[pathogen_col]] == parts[2L])
  }), hp_keys)

  imputation_method <- if (identical(residual_structure, "correlated")) {
    "conditional_gibbs_correlated"
  } else {
    "conditional_analytic_identity"
  }

  # -- Per-hp-key static info: panel, observed matrix, cohort masks ------------
  key_info <- stats::setNames(lapply(hp_keys, function(key) {
    d_hp <- match(hp_panel[[key]]$classes, class_cols)
    ev_idx <- hp_ev_idx[[key]]
    n_hp <- length(ev_idx)
    if (n_hp == 0L || length(d_hp) < 1L) {
      return(list(n_hp = n_hp, Dp = length(d_hp)))
    }

    panel_classes <- class_cols[d_hp]
    enum_df <- enumerate_binary_profiles(panel_classes)
    profile_bin <- as.matrix(enum_df[, panel_classes, drop = FALSE]) # n_profiles x Dp, 0/1
    obs_panel <- obs_ast_mat[ev_idx, d_hp, drop = FALSE] # n_hp x Dp, 0/1/NA

    list(
      d_hp = d_hp,
      ev_idx = ev_idx,
      n_hp = n_hp,
      Dp = length(d_hp),
      class_set = paste(panel_classes, collapse = "|"),
      all_lbls = enum_df$profile_delta,
      n_profiles = length(enum_df$profile_delta),
      profile_bin = profile_bin,
      obs_panel = obs_panel,
      n_missing_per_event = rowSums(is.na(obs_panel)),
      n_observed_per_event = length(d_hp) - rowSums(is.na(obs_panel)),
      known_out = is_known_outcome[ev_idx],
      nonfatal = is_nonfatal[ev_idx],
      n_events_all = n_hp,
      n_events_known = sum(is_known_outcome[ev_idx]),
      n_events_nf = sum(is_nonfatal[ev_idx]),
      n_profile_classes = length(d_hp),
      panel_eligibility_method = hp_panel[[key]]$method,
      classes_excluded = paste(hp_panel[[key]]$excluded, collapse = "|"),
      classes_excluded_reason = hp_panel[[key]]$reason
    )
  }), hp_keys)

  # Accumulator: posterior-mean event-profile sum. Shared by both residual
  # structures now -- identity computes prob_mat analytically, correlated
  # computes it via Gibbs-MC counting (see .gibbs_conditional_profile_probs()
  # below), but both produce an n_hp x n_profiles matrix whose rows sum to 1
  # and are consistent with observed cells, so the same accumulation and
  # sum-to-1 check apply uniformly to both.
  event_prob_sum <- stats::setNames(lapply(hp_keys, function(key) {
    ki <- key_info[[key]]
    if (ki$n_hp == 0L || ki$Dp < 1L) {
      return(NULL)
    }
    matrix(0, ki$n_hp, ki$n_profiles)
  }), hp_keys)

  message(sprintf(
    "[compute_event_profile_probabilities] %d draws | %d events | %d hp-pairs | %d RE level(s) | method=%s",
    S, N_ev, length(hp_keys), n_re_levels, imputation_method
  ))

  # hp-keys contributing rows, precomputed once to preallocate both result lists.
  valid_keys <- hp_keys[vapply(hp_keys, function(key) {
    ki <- key_info[[key]]
    ki$n_hp > 0L && ki$Dp >= 1L
  }, logical(1))]

  event_profile_rows <- vector("list", sum(vapply(valid_keys, function(key) key_info[[key]]$n_hp, integer(1))))
  aggregate_draw_rows <- vector("list", S * length(valid_keys))
  .event_i <- 0L
  .agg_i <- 0L

  # -- Main draw loop: mu_all computed once per draw, shared across hp-keys ---
  for (s in seq_len(S)) {
    beta_s <- matrix(beta_arr[s, , ], nrow = K, ncol = D)
    re_term <- if (re_prep$R > 0L) {
      re_eff_s <- matrix(re_eff_arr[s, , ], nrow = D, ncol = re_prep$total_re_levels)
      re_contribution(re_eff_s, flat_re_idx_obs)
    } else {
      matrix(0, nrow = N_ev, ncol = D)
    }
    mu_all <- (X_event %*% beta_s) + re_term # N_ev x D

    Omega_s_full <- if (identical(residual_structure, "correlated")) {
      tcrossprod(matrix(L_omega_arr[s, , ], nrow = D, ncol = D))
    } else {
      NULL
    }

    for (key in valid_keys) {
      ki <- key_info[[key]]
      mu_hp <- mu_all[ki$ev_idx, ki$d_hp, drop = FALSE] # n_hp x Dp

      if (identical(residual_structure, "correlated")) {
        # Gibbs-samples missing latent dims conditional on observed R/S signs
        # (not exact latent values), truncated to those signs every iteration.
        Omega_hp <- Omega_s_full[ki$d_hp, ki$d_hp, drop = FALSE]
        Omega_hp <- (Omega_hp + t(Omega_hp)) / 2 + diag(1e-9, ki$Dp)
        prob_mat <- .gibbs_conditional_profile_probs(
          mu_hp = mu_hp, Omega_hp = Omega_hp, obs_panel = ki$obs_panel,
          profile_bin = ki$profile_bin, n_burnin = n_gibbs_burnin, n_kept = n_gibbs_kept
        )
      } else {
        # -- Identity residual: exact analytic imputation over missing cells only.
        # Classes are conditionally independent given mu, so no latent-variable
        # simulation is needed at all.
        p_miss <- stats::pnorm(mu_hp) # marginal P(R) per class, in (0,1)
        # Effective per-class probability: the observed 0/1 value where tested,
        # Phi(mu) where untested. For a profile consistent with the observed
        # cells every "known" factor below evaluates to exactly 1; for an
        # inconsistent profile it evaluates to exactly 0 -- the same product
        # both enforces the hard observed-cell constraint and computes the
        # missing-cell probability, with no log(0) / 0 * -Inf edge cases because
        # we never take logs.
        p_eff <- ki$obs_panel
        miss_mask <- is.na(ki$obs_panel)
        p_eff[miss_mask] <- p_miss[miss_mask]

        prob_mat <- matrix(1, ki$n_hp, ki$n_profiles)
        for (d in seq_len(ki$Dp)) {
          col_is1 <- ki$profile_bin[, d] == 1L
          f_d <- matrix(NA_real_, ki$n_hp, ki$n_profiles)
          if (any(col_is1)) f_d[, col_is1] <- p_eff[, d]
          if (any(!col_is1)) f_d[, !col_is1] <- 1 - p_eff[, d]
          prob_mat <- prob_mat * f_d
        }
      }

      # Hard check, per draw, before any averaging: every event's enumerated
      # profile probabilities must sum to 1 (observed cells fix a hard subset
      # to consistent/inconsistent, missing cells are a proper Bernoulli
      # product) -- a violation indicates a bug in the panel/observed-cell
      # logic, not a data issue, so this stops rather than silently
      # renormalising.
      row_sums <- rowSums(prob_mat)
      bad_rows <- which(abs(row_sums - 1) > 1e-6)
      if (length(bad_rows) > 0L) {
        stop(sprintf(
          paste0(
            "[compute_event_profile_probabilities] Profile probabilities do not sum to 1 ",
            "for %d event(s) in hospital-pathogen pair '%s' at posterior draw %d ",
            "(max abs deviation = %.3e). This indicates a bug in the panel/observed-cell ",
            "logic and must be fixed before results are trusted."
          ),
          length(bad_rows), key, s, max(abs(row_sums[bad_rows] - 1))
        ), call. = FALSE)
      }

      event_prob_sum[[key]] <- event_prob_sum[[key]] + prob_mat

      .agg_i <- .agg_i + 1L
      aggregate_draw_rows[[.agg_i]] <- tibble::tibble(
        !!upper_re_col := hp_pairs[[upper_re_col]][match(key, hp_keys)],
        !!pathogen_col := hp_pairs[[pathogen_col]][match(key, hp_keys)],
        profile_class_set = ki$class_set,
        profile_delta = ki$all_lbls,
        draw_s = s,
        R_ALL_s = colMeans(prob_mat),
        R_KNOWN_OUTCOME_s = if (any(ki$known_out)) {
          colMeans(prob_mat[ki$known_out, , drop = FALSE])
        } else {
          rep(NA_real_, ki$n_profiles)
        },
        R_NF_s = if (any(ki$nonfatal)) {
          colMeans(prob_mat[ki$nonfatal, , drop = FALSE])
        } else {
          rep(NA_real_, ki$n_profiles)
        },
        profile_generation_method = imputation_method,
        n_profile_classes = ki$n_profile_classes,
        panel_eligibility_method = ki$panel_eligibility_method,
        classes_excluded = ki$classes_excluded,
        classes_excluded_reason = ki$classes_excluded_reason,
        n_events_all = ki$n_events_all,
        n_events_known_outcome = ki$n_events_known,
        n_events_nonfatal = ki$n_events_nf
      )
    }
  }

  # -- Finalize identity-residual event-level posterior means ------------------
  for (key in valid_keys) {
    ki <- key_info[[key]]
    h_nm <- hp_pairs[[upper_re_col]][match(key, hp_keys)]
    k_nm <- hp_pairs[[pathogen_col]][match(key, hp_keys)]
    event_prob_mean <- event_prob_sum[[key]] / S
    for (ev_i in seq_len(ki$n_hp)) {
      .event_i <- .event_i + 1L
      event_profile_rows[[.event_i]] <- tibble::tibble(
        !!upper_re_col := h_nm,
        !!pathogen_col := k_nm,
        event_idx = event_meta_obs$.event_idx[ki$ev_idx[ev_i]],
        profile_class_set = ki$class_set,
        profile_delta = ki$all_lbls,
        profile_probability = event_prob_mean[ev_i, ],
        profile_generation_method = imputation_method,
        n_profile_classes = ki$n_profile_classes,
        panel_eligibility_method = ki$panel_eligibility_method,
        classes_excluded = ki$classes_excluded,
        classes_excluded_reason = ki$classes_excluded_reason,
        n_classes_observed = ki$n_observed_per_event[ev_i],
        n_classes_missing = ki$n_missing_per_event[ev_i],
        fully_model_imputed = ki$n_observed_per_event[ev_i] == 0L
      )
    }
  }

  message(sprintf(
    "[compute_event_profile_probabilities] Done. %d event-profile rows | %d draw-aggregate rows.",
    length(event_profile_rows), length(aggregate_draw_rows)
  ))

  list(
    event_profiles  = dplyr::bind_rows(event_profile_rows),
    aggregate_draws = dplyr::bind_rows(aggregate_draw_rows)
  )
}

#' Assess Numerical Stability of Correlated Profile Completion
#'
#' Repeats the existing conditional truncated-MVN Gibbs profile completion for
#' fixed posterior draws and a fixed, reproducibly selected subset of incomplete
#' events.  It does not refit Stan or alter the fitted model.  The comparison
#' therefore isolates finite Gibbs-chain length.  The conditional completion
#' algorithm is the correct profile-completion algorithm; this diagnostic checks
#' whether its finite Monte Carlo run length is adequate for a particular fit,
#' especially when posterior residual correlations are high.
#'
#' The default 0.01/0.02 thresholds are numerical-stability conventions for
#' this diagnostic, not universal statistical thresholds.
#'
#' @param fitted_model Object returned by [fit_bayesian_multivariate_probit()].
#' @param schedules Named list of `list(burnin=, kept=)` schedules.  It must
#'   contain `baseline`, `medium`, and `long`.
#' @param n_posterior_draws Number of fixed posterior draws to use.
#' @param max_events Maximum number of incomplete events to assess.
#' @param seed Seed used once to select draws/events and then reset before each
#'   schedule.  Thus schedules share posterior draws, events, and RNG streams.
#' @param tolerance Maximum aggregate difference for status `"pass"`; twice
#'   this value is the `"warning"` upper bound.
#' @param event_selection Currently `"stratified"` or `"all"`.
#' @return A structured list with summary, event/aggregate comparisons, selected
#'   events, schedules, posterior draw indices, and a separate stability status.
#' @export
assess_gibbs_profile_stability <- function(
  fitted_model,
  schedules = list(
    baseline = list(burnin = 30L, kept = 50L),
    medium = list(burnin = 60L, kept = 100L),
    long = list(burnin = 100L, kept = 250L)
  ),
  n_posterior_draws = 200L,
  max_events = 250L,
  seed = 123L,
  tolerance = 0.01,
  event_selection = c("stratified", "all")
) {
  event_selection <- match.arg(event_selection)
  required <- c("baseline", "medium", "long")
  if (!all(required %in% names(schedules))) {
    stop("`schedules` must contain named baseline, medium, and long schedules.", call. = FALSE)
  }
  schedules <- schedules[required]
  for (nm in required) {
    x <- schedules[[nm]]
    if (!is.list(x) || !all(c("burnin", "kept") %in% names(x)) ||
        length(x$burnin) != 1L || length(x$kept) != 1L ||
        is.na(x$burnin) || is.na(x$kept) || x$burnin < 0 || x$kept < 1) {
      stop(sprintf("Schedule `%s` must supply non-negative `burnin` and positive `kept`.", nm), call. = FALSE)
    }
    schedules[[nm]] <- list(burnin = as.integer(x$burnin), kept = as.integer(x$kept))
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be one positive number.", call. = FALSE)
  }

  residual_structure <- .null_default(fitted_model$residual_structure, "identity")
  empty <- tibble::tibble()
  schedule_summary <- tibble::tibble(
    residual_structure = residual_structure,
    n_selected_events = 0L,
    n_posterior_draws = 0L,
    baseline_burnin = schedules$baseline$burnin, baseline_kept = schedules$baseline$kept,
    medium_burnin = schedules$medium$burnin, medium_kept = schedules$medium$kept,
    long_burnin = schedules$long$burnin, long_kept = schedules$long$kept,
    max_abs_aggregate_difference_baseline_long = NA_real_,
    mean_abs_aggregate_difference_baseline_long = NA_real_,
    max_abs_aggregate_difference_medium_long = NA_real_,
    mean_abs_aggregate_difference_medium_long = NA_real_,
    p95_event_max_abs_diff_baseline_long = NA_real_,
    max_event_max_abs_diff_baseline_long = NA_real_,
    gibbs_stability_status = "not_applicable"
  )
  if (!identical(residual_structure, "correlated")) {
    return(list(summary = schedule_summary, event_comparison = empty,
      aggregate_comparison = empty, selected_events = empty,
      schedules_used = schedules, posterior_draw_indices = integer(), seed = as.integer(seed),
      tolerance = tolerance, status = "not_applicable",
      rng_scheme = "Not applicable: identity-residual completion is analytic."))
  }

  event_meta <- fitted_model$event_metadata
  class_cols <- fitted_model$class_cols
  has_obs <- event_meta$.event_idx %in% fitted_model$data_long$ev_idx
  candidate <- event_meta[has_obs, , drop = FALSE]
  n_missing <- rowSums(is.na(candidate[, class_cols, drop = FALSE]))
  candidate <- candidate[n_missing > 0L, , drop = FALSE]
  n_missing <- n_missing[n_missing > 0L]
  if (!nrow(candidate)) {
    return(list(summary = schedule_summary, event_comparison = empty,
      aggregate_comparison = empty, selected_events = empty,
      schedules_used = schedules, posterior_draw_indices = integer(), seed = as.integer(seed),
      tolerance = tolerance, status = "not_applicable",
      rng_scheme = "Not applicable: no incomplete observed events."))
  }
  observed_n <- rowSums(!is.na(candidate[, class_cols, drop = FALSE]))
  observed_r <- rowSums(candidate[, class_cols, drop = FALSE] == 1, na.rm = TRUE)
  pattern <- ifelse(observed_n == 0L, "no_observed_class",
    ifelse(observed_r / observed_n > 0.5, "mostly_R",
      ifelse(observed_r / observed_n < 0.5, "mostly_S", "mixed_RS")))
  selected <- tibble::as_tibble(candidate) %>%
    dplyr::mutate(n_classes_missing = n_missing,
      missingness_stratum = ifelse(.data$n_classes_missing == 1L, "missing_1", "missing_2_plus"),
      observed_pattern_stratum = pattern)
  if (event_selection == "stratified" && nrow(selected) > max_events) {
    grp_cols <- c(fitted_model$upper_re_col, "missingness_stratum", "observed_pattern_stratum")
    key <- do.call(paste, c(selected[grp_cols], sep = "||"))
    groups <- split(seq_len(nrow(selected)), key)
    set.seed(as.integer(seed))
    quota <- max(1L, floor(max_events / length(groups)))
    take <- unlist(lapply(groups, function(ii) sample(ii, min(length(ii), quota))), use.names = FALSE)
    if (length(take) < max_events) {
      rest <- setdiff(seq_len(nrow(selected)), take)
      take <- c(take, sample(rest, min(length(rest), max_events - length(take))))
    }
    selected <- selected[sort(take), , drop = FALSE]
  } else if (nrow(selected) > max_events) {
    set.seed(as.integer(seed)); selected <- selected[sort(sample.int(nrow(selected), max_events)), , drop = FALSE]
  }
  selected_events <- selected %>% dplyr::transmute(event_idx = .data$.event_idx,
    n_classes_missing = .data$n_classes_missing,
    missingness_stratum = .data$missingness_stratum,
    observed_pattern_stratum = .data$observed_pattern_stratum,
    !!fitted_model$upper_re_col := .data[[fitted_model$upper_re_col]],
    !!fitted_model$pathogen_col := .data[[fitted_model$pathogen_col]])

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)
  }
  total_draws <- nrow(posterior::as_draws_matrix(fitted_model$draws))
  set.seed(as.integer(seed) + 1L)
  n_use <- min(as.integer(n_posterior_draws), total_draws)
  draw_idx <- if (n_use < total_draws) sort(sample.int(total_draws, n_use)) else seq_len(total_draws)

  # Resetting to the same seed makes each schedule consume the same initial RNG
  # stream.  Together with explicit `draw_idx`/event IDs this is common-random-
  # number comparison; only burn-in/retained Gibbs iterations differ.
  outputs <- lapply(names(schedules), function(nm) {
    set.seed(as.integer(seed) + 2L)
    compute_event_profile_probabilities(fitted_model,
      seed = as.integer(seed) + 2L, n_gibbs_burnin = schedules[[nm]]$burnin,
      n_gibbs_kept = schedules[[nm]]$kept, posterior_draw_indices = draw_idx,
      event_indices = selected_events$event_idx)
  })
  names(outputs) <- names(schedules)
  compare <- function(a, b, label) {
    ae <- outputs[[a]]$event_profiles %>%
      dplyr::select(event_idx, profile_class_set, profile_delta, p_a = profile_probability) %>%
      dplyr::inner_join(outputs[[b]]$event_profiles %>%
        dplyr::select(event_idx, profile_class_set, profile_delta, p_b = profile_probability),
        by = c("event_idx", "profile_class_set", "profile_delta")) %>%
      dplyr::group_by(.data$event_idx, .data$profile_class_set) %>%
      dplyr::summarise(comparison = label, max_abs_diff = max(abs(.data$p_a - .data$p_b)),
        mean_abs_diff = mean(abs(.data$p_a - .data$p_b)), l1_distance = sum(abs(.data$p_a - .data$p_b)),
        .groups = "drop")
    aa <- outputs[[a]]$aggregate_draws %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(fitted_model$upper_re_col, fitted_model$pathogen_col,
        "profile_class_set", "profile_delta")))) %>%
      dplyr::summarise(p_a = mean(.data$R_ALL_s), .groups = "drop") %>%
      dplyr::inner_join(outputs[[b]]$aggregate_draws %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(fitted_model$upper_re_col, fitted_model$pathogen_col,
          "profile_class_set", "profile_delta")))) %>%
        dplyr::summarise(p_b = mean(.data$R_ALL_s), .groups = "drop"),
        by = c(fitted_model$upper_re_col, fitted_model$pathogen_col, "profile_class_set", "profile_delta")) %>%
      dplyr::transmute(dplyr::across(dplyr::all_of(c(fitted_model$upper_re_col, fitted_model$pathogen_col,
        "profile_class_set", "profile_delta"))), comparison = label, abs_diff = abs(.data$p_a - .data$p_b))
    list(event = ae, aggregate = aa)
  }
  bl <- compare("baseline", "long", "baseline_vs_long")
  ml <- compare("medium", "long", "medium_vs_long")
  event_comparison <- dplyr::bind_rows(bl$event, ml$event)
  aggregate_comparison <- dplyr::bind_rows(bl$aggregate, ml$aggregate)
  max_agg <- max(bl$aggregate$abs_diff, na.rm = TRUE)
  status <- if (max_agg <= tolerance) "pass" else if (max_agg <= 2 * tolerance) "warning" else "fail"
  summary <- schedule_summary %>% dplyr::mutate(
    n_selected_events = nrow(selected_events), n_posterior_draws = length(draw_idx),
    max_abs_aggregate_difference_baseline_long = max_agg,
    mean_abs_aggregate_difference_baseline_long = mean(bl$aggregate$abs_diff),
    max_abs_aggregate_difference_medium_long = max(ml$aggregate$abs_diff, na.rm = TRUE),
    mean_abs_aggregate_difference_medium_long = mean(ml$aggregate$abs_diff),
    p95_event_max_abs_diff_baseline_long = stats::quantile(bl$event$max_abs_diff, 0.95, names = FALSE),
    max_event_max_abs_diff_baseline_long = max(bl$event$max_abs_diff), gibbs_stability_status = status)
  list(summary = summary, event_comparison = event_comparison,
    aggregate_comparison = aggregate_comparison, selected_events = selected_events,
    schedules_used = schedules, posterior_draw_indices = draw_idx, seed = as.integer(seed),
    tolerance = tolerance, status = status,
    rng_scheme = "Shared explicit posterior draw indices/event IDs; RNG reset to seed + 2 before each schedule.")
}

#' Aggregate Posterior Profile Draws into R_ALL / R_KNOWN_OUTCOME / R_NF Summaries
#'
#' Summarises draw-level R_ALL, R_KNOWN_OUTCOME, and R_NF values from
#' \code{compute_event_profile_probabilities()} into posterior mean and
#' credible interval per hospital x pathogen x profile combination, and
#' derives explicit \code{eligible_for_YLL}/\code{eligible_for_YLD} flags
#' (with reasons) rather than asserting universal usability.
#'
#' \strong{Cohort definitions:} \code{R_ALL} is the mean over truly all events
#' in the hospital-pathogen-panel cohort (no outcome filter). \code{R_KNOWN_OUTCOME}
#' restricts to events with a known patient outcome (this is what earlier
#' versions of this function called \code{R_ALL} -- renamed because it is not
#' actually all events). \code{R_NF} restricts to the non-fatal subset of the
#' known-outcome cohort. All three are computed per posterior draw and then
#' summarised here; \code{n_draws_all}/\code{n_draws_known_outcome}/\code{n_draws_nf}
#' count only the non-\code{NA} draws feeding each summary (a cohort that is
#' empty for a given hospital-pathogen pair yields \code{NA} draws, which must
#' not be silently counted as valid support).
#'
#' \strong{Eligibility:} \code{eligible_for_profile_inference} is a hard
#' boolean, \code{TRUE} only for rows generated by a conditional
#' (observed-plus-imputed) method -- \code{"conditional_analytic_identity"}
#' or \code{"conditional_gibbs_correlated"}; see
#' \code{compute_event_profile_probabilities()}. \code{FALSE} for any other
#' \code{profile_generation_method} value, including the legacy
#' \code{"unconditional_simulation_correlated_not_daly_eligible"} tag from an
#' older \code{profile_output} -- callers must not rely on the descriptive
#' \code{profile_generation_method} string alone.
#' \code{eligible_for_YLL}/\code{eligible_for_YLD} additionally require:
#' \code{sampler_acceptable} (passed in by the caller, e.g. from
#' \code{fitted_model$diagnostics$converged_structural}); the relevant event
#' cohort has at least \code{min_n_events} events; and at least
#' \code{min_n_draws} valid posterior draws contributed to the relevant
#' summary. \code{exclusion_reason_YLL}/\code{exclusion_reason_YLD} record
#' why, distinguishing a cohort that is completely empty
#' (\code{"no_known_outcome_events"} / \code{"no_nonfatal_events"}) from one
#' that is merely small (\code{"too_few_known_outcome_events"} /
#' \code{"too_few_nonfatal_events"}).
#'
#' @param profile_output List returned by \code{compute_event_profile_probabilities()}.
#' @param hospital_col Character. Hospital column in \code{aggregate_draws}.
#'   Must match the upper RE column used during fitting. Default \code{"hospital"}.
#' @param pathogen_col Character. Pathogen column. Default \code{"pathogen"}.
#' @param estimand Character. Estimand label to attach to output.
#'   Default \code{"observed_stewardship_event_mix"}.
#' @param ci_level Numeric. Credible interval coverage. Default \code{0.95}.
#' @param min_n_events Integer. Minimum events in the relevant cohort
#'   (known-outcome for YLL, nonfatal for YLD) for eligibility. Default \code{10L}.
#' @param min_n_draws Integer. Minimum valid posterior draws in the relevant
#'   summary for eligibility. Default \code{100L}.
#' @param sampler_acceptable Logical, scalar or one value per row of
#'   \code{profile_output$aggregate_draws} groups. Whether the fit's sampler
#'   diagnostics are acceptable (e.g. \code{fitted_model$diagnostics$converged_structural}
#'   from \code{fit_bayesian_multivariate_probit()}); \code{FALSE} makes every
#'   row ineligible for YLL/YLD regardless of panel/count support. Default
#'   \code{TRUE} (i.e. does not gate on sampler status unless the caller
#'   supplies it -- \code{estimate_resistance_profiles()} supplies it
#'   automatically).
#'
#' @return Tibble with one row per hospital x pathogen x profile: R_ALL,
#'   R_KNOWN_OUTCOME, and R_NF posterior mean and credible interval, panel
#'   composition (\code{n_profile_classes}, \code{panel_eligibility_method},
#'   \code{classes_excluded}, \code{classes_excluded_reason}), event and draw
#'   counts, and profile/YLL/YLD eligibility flags with reasons.
#' @export
aggregate_profiles_for_daly <- function(
  profile_output,
  hospital_col = "hospital",
  pathogen_col = "pathogen",
  estimand = "observed_stewardship_event_mix",
  ci_level = 0.95,
  min_n_events = 10L,
  min_n_draws = 100L,
  sampler_acceptable = TRUE
) {
  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q

  if (!is.list(profile_output) ||
    !all(c("aggregate_draws", "event_profiles") %in% names(profile_output))) {
    stop("`profile_output` must be the list returned by compute_event_profile_probabilities().",
      call. = FALSE
    )
  }

  agg_draws <- profile_output$aggregate_draws
  if (nrow(agg_draws) == 0L) {
    warning("aggregate_draws is empty. Returning empty tibble.", call. = FALSE)
    return(tibble::tibble())
  }

  # Detect the actual hospital column name from the draws tibble
  if (!hospital_col %in% names(agg_draws)) {
    non_hospital_cols <- c(
      pathogen_col, "profile_class_set", "profile_delta", "draw_s",
      "R_ALL_s", "R_KNOWN_OUTCOME_s", "R_NF_s",
      "profile_generation_method",
      "n_profile_classes", "panel_eligibility_method",
      "classes_excluded", "classes_excluded_reason",
      "n_events_all", "n_events_known_outcome", "n_events_nonfatal"
    )
    candidate <- setdiff(names(agg_draws), non_hospital_cols)
    if (length(candidate) == 1L) {
      hospital_col <- candidate[1L]
      message(sprintf(
        "[aggregate_profiles_for_daly] Using '%s' as hospital column.",
        hospital_col
      ))
    } else {
      stop(sprintf("hospital_col '%s' not found in aggregate_draws.", hospital_col),
        call. = FALSE
      )
    }
  }

  result <- agg_draws %>%
    dplyr::group_by(
      .data[[hospital_col]],
      .data[[pathogen_col]],
      .data$profile_class_set,
      .data$profile_delta
    ) %>%
    dplyr::summarise(
      R_ALL_mean = mean(.data$R_ALL_s, na.rm = TRUE),
      R_ALL_lower = stats::quantile(.data$R_ALL_s, lo_q, na.rm = TRUE),
      R_ALL_upper = stats::quantile(.data$R_ALL_s, hi_q, na.rm = TRUE),
      R_KNOWN_OUTCOME_mean = mean(.data$R_KNOWN_OUTCOME_s, na.rm = TRUE),
      R_KNOWN_OUTCOME_lower = stats::quantile(.data$R_KNOWN_OUTCOME_s, lo_q, na.rm = TRUE),
      R_KNOWN_OUTCOME_upper = stats::quantile(.data$R_KNOWN_OUTCOME_s, hi_q, na.rm = TRUE),
      R_NF_mean = mean(.data$R_NF_s, na.rm = TRUE),
      R_NF_lower = stats::quantile(.data$R_NF_s, lo_q, na.rm = TRUE),
      R_NF_upper = stats::quantile(.data$R_NF_s, hi_q, na.rm = TRUE),
      # n_draws_* count valid POSTERIOR DRAWS contributing to this profile's
      # summary, not event-profile rows -- each retained fitted state
      # contributes exactly one row per (hospital, pathogen, profile_class_set,
      # profile_delta) group here, so this is a draw count by construction.
      # Event-level detail lives only in profile_output$event_profiles.
      n_draws_all = sum(!is.na(.data$R_ALL_s)),
      n_draws_known_outcome = sum(!is.na(.data$R_KNOWN_OUTCOME_s)),
      n_draws_nf = sum(!is.na(.data$R_NF_s)),
      n_events_all = dplyr::first(.data$n_events_all),
      n_events_known_outcome = dplyr::first(.data$n_events_known_outcome),
      n_events_nonfatal = dplyr::first(.data$n_events_nonfatal),
      profile_generation_method = dplyr::first(.data$profile_generation_method),
      n_profile_classes = dplyr::first(.data$n_profile_classes),
      panel_eligibility_method = dplyr::first(.data$panel_eligibility_method),
      classes_excluded = dplyr::first(.data$classes_excluded),
      classes_excluded_reason = dplyr::first(.data$classes_excluded_reason),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # mean(all-NA, na.rm = TRUE) returns NaN (not NA); normalise so "no data"
      # is represented consistently across mean/lower/upper.
      dplyr::across(
        c(.data$R_ALL_mean, .data$R_KNOWN_OUTCOME_mean, .data$R_NF_mean),
        ~ dplyr::if_else(is.nan(.x), NA_real_, .x)
      ),
      profile_label = .data$profile_delta,
      profile_set_type = "facility_bayesian_probit",
      estimand = estimand,
      estimator = "bayesian_multivariate_probit",
      sampler_acceptable = sampler_acceptable,
      # Only the two conditional (observed-plus-imputed) generation methods
      # are DALY-eligible; any other value, including the legacy
      # "unconditional_simulation_correlated_not_daly_eligible" tag, is not.
      eligible_for_profile_inference =
        .data$profile_generation_method %in%
          c("conditional_analytic_identity", "conditional_gibbs_correlated"),
      low_draws_all_flag = .data$n_draws_all < min_n_draws,
      low_draws_known_outcome_flag = .data$n_draws_known_outcome < min_n_draws,
      low_draws_nf_flag = .data$n_draws_nf < min_n_draws,
      low_events_flag = .data$n_events_all < min_n_events,
      eligible_for_YLL =
        .data$eligible_for_profile_inference &
          .data$sampler_acceptable &
          .data$n_events_known_outcome >= min_n_events &
          !.data$low_draws_known_outcome_flag,
      eligible_for_YLD =
        .data$eligible_for_profile_inference &
          .data$sampler_acceptable &
          .data$n_events_nonfatal >= min_n_events &
          !.data$low_draws_nf_flag,
      exclusion_reason_YLL = dplyr::case_when(
        .data$eligible_for_YLL ~ NA_character_,
        !.data$eligible_for_profile_inference ~ .data$profile_generation_method,
        !.data$sampler_acceptable ~ "sampler_not_acceptable",
        .data$n_events_known_outcome == 0L ~ "no_known_outcome_events",
        .data$n_events_known_outcome < min_n_events ~ "too_few_known_outcome_events",
        .data$low_draws_known_outcome_flag ~ "too_few_valid_draws",
        TRUE ~ "not_eligible"
      ),
      exclusion_reason_YLD = dplyr::case_when(
        .data$eligible_for_YLD ~ NA_character_,
        !.data$eligible_for_profile_inference ~ .data$profile_generation_method,
        !.data$sampler_acceptable ~ "sampler_not_acceptable",
        .data$n_events_nonfatal == 0L ~ "no_nonfatal_events",
        .data$n_events_nonfatal < min_n_events ~ "too_few_nonfatal_events",
        .data$low_draws_nf_flag ~ "too_few_valid_draws",
        TRUE ~ "not_eligible"
      ),
      # eligible_for_YLL/YLD (above) only ever required MARGINAL support --
      # profile generation succeeding, the sampler being acceptable, and
      # enough events/draws. They say nothing about whether the model's
      # *joint* (multi-class) resistance profile can be trusted, which the
      # exp_07 finding showed is a materially different question (excellent
      # marginal calibration, poor pairwise/complete-profile calibration
      # under an identity residual). eligible_for_marginal_YLL/YLD makes that
      # existing, narrower meaning explicit under its own name.
      # eligible_for_joint_profile_YLL/YLD starts identical to the marginal
      # flags here (validation hasn't run yet at profile-generation time) and
      # is subsequently AND-ed down by generate_bayesian_probit_validation()
      # once pairwise/complete-profile calibration results are available --
      # mirroring how eligible_for_YLL/YLD is itself narrowed post-validation.
      eligible_for_marginal_YLL = .data$eligible_for_YLL,
      eligible_for_marginal_YLD = .data$eligible_for_YLD,
      eligible_for_joint_profile_YLL = .data$eligible_for_YLL,
      eligible_for_joint_profile_YLD = .data$eligible_for_YLD,
      exclusion_reason_joint_profile_YLL = .data$exclusion_reason_YLL,
      exclusion_reason_joint_profile_YLD = .data$exclusion_reason_YLD,
      dplyr::across(dplyr::starts_with("R_"), ~ round(.x, 6L))
    ) %>%
    dplyr::select(
      dplyr::all_of(c(hospital_col, pathogen_col)),
      profile_set_type, estimand, profile_class_set, profile_delta, profile_label,
      R_ALL_mean, R_ALL_lower, R_ALL_upper,
      R_KNOWN_OUTCOME_mean, R_KNOWN_OUTCOME_lower, R_KNOWN_OUTCOME_upper,
      R_NF_mean, R_NF_lower, R_NF_upper,
      n_profile_classes, panel_eligibility_method, classes_excluded, classes_excluded_reason,
      n_events_all, n_events_known_outcome, n_events_nonfatal,
      n_draws_all, n_draws_known_outcome, n_draws_nf,
      profile_generation_method, sampler_acceptable,
      eligible_for_profile_inference,
      eligible_for_YLL, exclusion_reason_YLL,
      eligible_for_YLD, exclusion_reason_YLD,
      eligible_for_marginal_YLL, eligible_for_marginal_YLD,
      eligible_for_joint_profile_YLL, exclusion_reason_joint_profile_YLL,
      eligible_for_joint_profile_YLD, exclusion_reason_joint_profile_YLD,
      low_draws_all_flag, low_draws_known_outcome_flag, low_draws_nf_flag, low_events_flag,
      estimator
    )

  message(sprintf(
    "[aggregate_profiles_for_daly] %d hospital-pathogen-profile rows (%d eligible for YLL, %d eligible for YLD).",
    nrow(result), sum(result$eligible_for_YLL), sum(result$eligible_for_YLD)
  ))

  result
}


# estimate_resistance_profiles()  -- top-level dispatcher for both pathways

#' Estimate Resistance Profiles: Pathway 1 (Convex) or Pathway 2 (Bayesian)
#'
#' Top-level dispatcher. Runs either the convex optimisation pathway
#' (Pathway 1, aggregate surveillance data) or the Bayesian hierarchical
#' multivariate probit pathway (Pathway 2, facility-level AST data).
#'
#' For Pathway 2, pass \code{pathogen} to fit one pathogen at a time -- the
#' recommended workflow. Run this function once per pathogen and collect the
#' results in the analysis repository.
#'
#' @param data For Pathway 1: marginals tibble from
#'   \code{compute_marginals_from_data()}. For Pathway 2: wide-format
#'   event-level tibble (one row per event; class columns 0/1/NA).
#' @param method Character. \code{"convex"} or \code{"bayesian"}.
#' @param pairwise Tibble or \code{NULL}. Pathway 1 only.
#' @param panel_map Named list or \code{NULL}. Pathway 1 only.
#' @param class_cols Character vector. Pathway 2 only. Required.
#' @param fixed_effects Character vector. Pathway 2 only. Required.
#' @param random_effects Character vector or named list of random-intercept
#'   blocks. Pathway 2 only. Use \code{list()} for a fixed-effects-only model.
#' @param profile_group_col Character scalar or \code{NULL}. Pathway 2 only.
#'   Grouping column for profile aggregation and validation. Required when
#'   \code{random_effects} is empty.
#' @param pathogen Character or \code{NULL}. Pathway 2 only. When supplied,
#'   filters data to a single pathogen before fitting.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param eligible_pairs Tibble or \code{NULL}. Pathway 2 only.
#' @param outcome_col Character or \code{NULL}. Pathway 2 only. Patient
#'   outcome column for R_ALL / R_NF split.
#' @param nonfatal_values Character vector. Outcome values for R_NF cohort.
#' @param panel_eligibility Named list. Eligibility thresholds; passed to
#'   \code{fit_bayesian_multivariate_probit()}.
#' @param estimand Character. Estimand label. Default
#'   \code{"observed_stewardship_event_mix"}.
#' @param prior_config Named list. Pathway 2 only. Priors.
#' @param sampler_config Named list. Pathway 2 only. Structured sampler
#'   settings: \code{chains} (4), \code{iter_warmup} (1000),
#'   \code{iter_sampling} (1000), \code{seed} (123), \code{adapt_delta},
#'   \code{max_treedepth}, \code{parallel_chains}.
#' @param residual_structure Character. Pathway 2. \code{"identity"} (default)
#'   or \code{"correlated"}. See \code{fit_bayesian_multivariate_probit()}.
#' @param n_posterior_draws_for_profiles Integer. Pathway 2. Number of posterior
#'   draws used for MVN simulation of profiles. Default \code{2000L}.
#'
#' @return Named list: \code{profiles}, \code{eligibility}, \code{diagnostics},
#'   \code{fitted_models}, \code{config_used}.
#' @export
estimate_resistance_profiles <- function(
  data,
  method = c("convex", "bayesian"),
  # Pathway 1
  pairwise = NULL,
  panel_map = NULL,
  # Pathway 2 -- required with no defaults
  class_cols = NULL,
  fixed_effects = NULL,
  random_effects = list(),
  profile_group_col = NULL,
  pathogen = NULL,
  pathogen_col = "pathogen",
  eligible_pairs = NULL,
  outcome_col = NULL,
  nonfatal_values = c(
    "Discharged", "Survived", "Alive",
    "discharged", "survived", "alive"
  ),
  panel_eligibility = list(),
  residual_structure = c("identity", "correlated"),
  estimand = "observed_stewardship_event_mix",
  prior_config = list(),
  sampler_config = list(),
  n_posterior_draws_for_profiles = 2000L
) {
  method <- match.arg(method)
  residual_structure <- match.arg(residual_structure)

  config_used <- list(
    method                         = method,
    pathogen_col                   = pathogen_col,
    pathogen                       = pathogen,
    residual_structure             = residual_structure,
    random_effects                 = random_effects,
    profile_group_col              = profile_group_col,
    estimand                       = estimand,
    n_posterior_draws_for_profiles = n_posterior_draws_for_profiles,
    prior_config                   = prior_config,
    sampler_config                 = sampler_config
  )

  # ---- Pathway 1 -----------------------------------------------------------
  if (method == "convex") {
    message("[estimate_resistance_profiles] Running Pathway 1 (convex optimisation)...")

    profiles_tbl <- estimate_profiles_convex(
      marginals    = data,
      pairwise     = pairwise,
      panel_map    = panel_map,
      pathogen_col = pathogen_col
    )

    diag_tbl <- profiles_tbl %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        n_profiles          = dplyr::n(),
        pct_converged       = mean(.data$convergence_flag, na.rm = TRUE),
        max_abs_residual    = max(.data$max_abs_residual, na.rm = TRUE),
        any_underdetermined = any(.data$identifiability_flag, na.rm = TRUE),
        .groups             = "drop"
      )

    return(list(
      profiles      = profiles_tbl,
      eligibility   = NULL,
      diagnostics   = diag_tbl,
      fitted_models = list(convex_profiles = profiles_tbl),
      config_used   = config_used
    ))
  }

  # ---- Pathway 2 -----------------------------------------------------------
  message("[estimate_resistance_profiles] Running Pathway 2 (Bayesian multivariate probit)...")
  message(sprintf("  Estimand: %s", estimand))
  if (!is.null(pathogen)) {
    message(sprintf("  Pathogen: %s", pathogen))
  }

  if (is.null(class_cols)) {
    stop("`class_cols` is required for method = 'bayesian'.", call. = FALSE)
  }
  if (is.null(fixed_effects)) {
    stop("`fixed_effects` is required for method = 'bayesian'.", call. = FALSE)
  }

  fitted_mod <- fit_bayesian_multivariate_probit(
    event_class_data   = data,
    class_cols         = class_cols,
    fixed_effects      = fixed_effects,
    random_effects     = random_effects,
    profile_group_col  = profile_group_col,
    pathogen           = pathogen,
    pathogen_col       = pathogen_col,
    eligible_pairs     = eligible_pairs,
    outcome_col        = outcome_col,
    panel_eligibility  = panel_eligibility,
    residual_structure = residual_structure,
    estimand           = estimand,
    prior_config       = prior_config,
    sampler_config     = sampler_config
  )

  profile_probs <- compute_event_profile_probabilities(
    fitted_model                   = fitted_mod,
    n_posterior_draws_for_profiles = as.integer(n_posterior_draws_for_profiles),
    outcome_col                    = outcome_col,
    nonfatal_values                = nonfatal_values,
    seed                           = .null_default(sampler_config$seed, 123L)
  )

  profiles_tbl <- aggregate_profiles_for_daly(
    profile_output = profile_probs,
    hospital_col = fitted_mod$upper_re_col,
    pathogen_col = pathogen_col,
    estimand = estimand,
    sampler_acceptable = isTRUE(fitted_mod$diagnostics$converged_structural)
  )

  eligibility_tbl <- if (!is.null(eligible_pairs)) {
    eligible_pairs %>%
      dplyr::count(.data[[fitted_mod$upper_re_col]], name = "n_eligible_pathogens") %>%
      dplyr::left_join(
        profiles_tbl %>%
          dplyr::group_by(.data[[fitted_mod$upper_re_col]]) %>%
          dplyr::summarise(n_profiles_estimated = dplyr::n(), .groups = "drop"),
        by = fitted_mod$upper_re_col
      )
  } else {
    NULL
  }

  message("[estimate_resistance_profiles] Pathway 2 complete.")

  list(
    profiles      = profiles_tbl,
    eligibility   = eligibility_tbl,
    diagnostics   = fitted_mod$diagnostics,
    fitted_models = list(bayesian_probit = fitted_mod),
    config_used   = config_used
  )
}
