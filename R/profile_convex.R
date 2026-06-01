# profile_convex.R
# Pathway 1: Convex optimisation-based resistance profile estimation.
#
# Functions implemented here:
#   validate_profile_inputs()   -- column checks, format detection, stratification
#   preprocess_for_profiles()   -- full preprocessing pipeline
#
# Downstream functions (estimate_profiles_convex, bootstrap, etc.) will follow
# in this same file.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# validate_profile_inputs()
# ---------------------------------------------------------------------------

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
#'   data        = my_ast_data,
#'   col_map     = list(
#'     isolate_col   = "isolate_id",
#'     pathogen_col  = "organism",
#'     ast_col       = "result",
#'     patient_col   = "pid",
#'     date_col      = "culture_date",
#'     geography_col = "state",
#'     specimen_col  = "specimen",
#'     age_col       = "age_years",
#'     dob_col       = NULL,
#'     antibiotic_col = "drug_name",
#'     class_col     = NULL,
#'     location_col  = "ward",
#'     outcome_col   = "discharge_status"
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
  if (nrow(data) == 0L)    stop("`data` has zero rows.")

  results <- list()

  # ------------------------------------------------------------------
  # 1. Format detection
  # ------------------------------------------------------------------
  abx_col <- .null_default(col_map$antibiotic_col, "antibiotic_name")
  iso_col <- .null_default(col_map$isolate_col,    "isolate_id")

  has_abx_col   <- abx_col %in% names(data)
  has_iso_col   <- iso_col %in% names(data)

  if (has_abx_col && has_iso_col) {
    # Long format: each isolate appears in multiple rows (one per antibiotic)
    max_rows_per_iso <- max(table(data[[iso_col]]), na.rm = TRUE)
    detected_format  <- if (max_rows_per_iso > 1L) "long" else "wide"
  } else if (has_abx_col) {
    detected_format <- "long"
  } else {
    detected_format <- "wide"
  }

  results <- .add_check(
    results, "format_detection", "pass",
    sprintf("Data detected as '%s' format. %d rows x %d columns.",
            detected_format, nrow(data), ncol(data))
  )

  # ------------------------------------------------------------------
  # 2. Mandatory columns (all formats)
  # ------------------------------------------------------------------
  mandatory_always <- c(
    isolate_col   = .null_default(col_map$isolate_col,   "isolate_id"),
    pathogen_col  = .null_default(col_map$pathogen_col,  "pathogen"),
    ast_col       = .null_default(col_map$ast_col,       "ast_value"),
    patient_col   = .null_default(col_map$patient_col,   "patient_id"),
    date_col      = .null_default(col_map$date_col,      "date_of_culture"),
    geography_col = .null_default(col_map$geography_col, "state"),
    specimen_col  = .null_default(col_map$specimen_col,  "specimen_type")
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

  # ------------------------------------------------------------------
  # 3. Mandatory in long format: antibiotic_name column
  # ------------------------------------------------------------------
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

  # ------------------------------------------------------------------
  # 4. Age OR DOB -- at least one mandatory
  # ------------------------------------------------------------------
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
      c(if (has_age) sprintf("age ('%s')", age_col),
        if (has_dob) sprintf("dob ('%s')", dob_col)),
      collapse = " and "
    )
    results <- .add_check(
      results, "mandatory_age_or_dob", "pass",
      sprintf("Age/DOB check passed. Found: %s.", found_desc)
    )
  }

  # ------------------------------------------------------------------
  # 5. Optional columns -- warn but do not stop
  # ------------------------------------------------------------------
  eff_outcome_col <- .null_default(outcome_col, col_map$outcome_col)

  optional_checks <- list(
    class_col    = list(
      col = col_map$class_col,
      msg_absent = "Antibiotic class column '%s' absent. Will derive classes from antibiotic names via WHO AWaRe reference."
    ),
    location_col = list(
      col = col_map$location_col,
      msg_absent = "Location column '%s' absent. Ward/ICU/OPD subgroup outputs will not be available."
    ),
    outcome_col  = list(
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

  # ------------------------------------------------------------------
  # 6. AST value content check
  # ------------------------------------------------------------------
  ast_col_name <- .null_default(col_map$ast_col, "ast_value")
  if (ast_col_name %in% names(data)) {
    ast_raw   <- toupper(trimws(as.character(data[[ast_col_name]])))
    valid_sir <- c("S", "I", "R")
    n_valid   <- sum(ast_raw %in% valid_sir, na.rm = TRUE)
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

  # ------------------------------------------------------------------
  # 7. Isolate ID uniqueness (wide format only)
  # ------------------------------------------------------------------
  if (detected_format == "wide" && has_iso_col) {
    n_dup <- sum(duplicated(data[[iso_col]], incomparables = NA))
    if (n_dup > 0L) {
      results <- .add_check(
        results, "isolate_id_uniqueness", "warn",
        sprintf("%d duplicate isolate IDs in '%s' (wide format expects one row per isolate).",
                n_dup, iso_col),
        n_dup
      )
    } else {
      results <- .add_check(
        results, "isolate_id_uniqueness", "pass",
        sprintf("All isolate IDs in '%s' are unique.", iso_col)
      )
    }
  }

  # ------------------------------------------------------------------
  # 8. Stratification-specific checks
  # ------------------------------------------------------------------
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
      results <- .add_check(results, "outcome_col", "fail",
        sprintf("outcome_col '%s' requested but not found in data.", eff_outcome_col))
    } else {
      n_na_out <- sum(is.na(data[[eff_outcome_col]]))
      results <- .add_check(results, "outcome_col", "pass",
        sprintf("Outcome column '%s' found. %d NA value(s).", eff_outcome_col, n_na_out), n_na_out)
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
      extra    <- if (dim == "year")
        " Year will be extracted as integer from this column."
      else ""
      results <- .add_check(
        results, paste0("stratify_", dim), "pass",
        sprintf("Stratification by %s: column '%s' found with %d unique value(s).%s",
                dim, col, n_strata, extra),
        n_strata
      )
    }
  }

  # ------------------------------------------------------------------
  # 9. Collate, report, and stop/warn
  # ------------------------------------------------------------------
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
      sprintf("%d mandatory check(s) failed:\n%s",
              n_fail, paste(sprintf("  - %s", fail_msgs), collapse = "\n")),
      call. = FALSE
    )
  }

  attr(results_tbl, "detected_format") <- detected_format
  invisible(results_tbl)
}


# ---------------------------------------------------------------------------
# preprocess_for_profiles()
# ---------------------------------------------------------------------------

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
#'   data        = raw_ast,
#'   col_map     = list(
#'     isolate_col   = "isolate_id",
#'     pathogen_col  = "organism",
#'     ast_col       = "result",
#'     patient_col   = "pid",
#'     date_col      = "culture_date",
#'     geography_col = "state",
#'     specimen_col  = "specimen",
#'     age_col       = "age",
#'     antibiotic_col = "drug_name"
#'   ),
#'   panel_map   = list(
#'     "Klebsiella pneumoniae" = c("3GC", "Carbapenems", "Fluoroquinolones"),
#'     "Escherichia coli"      = c("3GC", "Fluoroquinolones", "Aminoglycosides")
#'   ),
#'   stratify_by = c("geography", "year"),
#'   outcome_col = "final_outcome"
#' )
#'
#' result$data_wide           # ready for compute_marginals_from_data()
#' result$preprocessing_log   # step-by-step row counts
#' result$panel_exclusions    # classes dropped per organism
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
    who_table   = NULL
) {

  # ------------------------------------------------------------------
  # Argument checks
  # ------------------------------------------------------------------
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
         call. = FALSE)
  }

  eff_outcome_col <- .null_default(outcome_col, col_map$outcome_col)

  # Resolve column names with defaults
  isolate_col   <- .null_default(col_map$isolate_col,    "isolate_id")
  pathogen_col  <- .null_default(col_map$pathogen_col,   "pathogen")
  ast_col       <- .null_default(col_map$ast_col,        "ast_value")
  patient_col   <- .null_default(col_map$patient_col,    "patient_id")
  date_col      <- .null_default(col_map$date_col,       "date_of_culture")
  geography_col <- .null_default(col_map$geography_col,  "state")
  specimen_col  <- .null_default(col_map$specimen_col,   "specimen_type")
  age_col       <- col_map$age_col
  dob_col       <- col_map$dob_col
  antibiotic_col <- .null_default(col_map$antibiotic_col, "antibiotic_name")
  class_col     <- col_map$class_col
  location_col  <- col_map$location_col

  plog       <- list()
  exclusions <- list()

  # ------------------------------------------------------------------
  # Step 0: Input validation
  # ------------------------------------------------------------------
  message("\n[preprocess_for_profiles] Step 0: Validating inputs...")
  validation_report <- validate_profile_inputs(
    data        = data,
    col_map     = col_map,
    stratify_by = stratify_by,
    outcome_col = eff_outcome_col
  )
  detected_format <- attr(validation_report, "detected_format")
  plog <- .log_step(plog, "0_validate", nrow(data), nrow(data),
                    sprintf("Validation passed. Format: %s.", detected_format))

  # ------------------------------------------------------------------
  # Step 1: Wide -> long conversion
  # ------------------------------------------------------------------
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
    known_meta    <- known_meta[!is.null(known_meta) & known_meta %in% names(data)]
    candidate_abx <- setdiff(names(data), known_meta)

    if (length(candidate_abx) == 0L) {
      stop("Wide format detected but no antibiotic/class columns found after excluding metadata columns.",
           call. = FALSE)
    }
    message(sprintf("  Wide -> long: pivoting %d antibiotic column(s).", length(candidate_abx)))

    data <- prep_pivot_ast_wide_to_long(
      data             = data,
      antibiotic_cols  = candidate_abx,
      id_cols          = known_meta,
      antibiotic_name_col  = antibiotic_col,
      antibiotic_value_col = ast_col,
      remove_missing   = FALSE   # keep NA rows; will be handled in harmonisation
    )
    plog <- .log_step(plog, "1_wide_to_long", n_in, nrow(data),
                      sprintf("Wide -> long: pivoted %d antibiotic column(s).", length(candidate_abx)))
  } else {
    plog <- .log_step(plog, "1_format", n_in, nrow(data), "Long format confirmed. No pivot needed.")
  }

  # ------------------------------------------------------------------
  # Step 2: AST harmonisation
  # ------------------------------------------------------------------
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
  plog <- .log_step(plog, "2_harmonise_ast", n_in, nrow(data),
                    sprintf("AST values harmonised into '%s'. I values recoded.", harmonised_ast_col))

  # ------------------------------------------------------------------
  # Step 3: Antibiotic name standardisation + class assignment
  # ------------------------------------------------------------------
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
      plog <- .log_step(plog, "3a_standardize_abx", n_in, nrow(data),
                        "Antibiotic names standardised via WHO AWaRe; antibiotic_class derived.")
    } else {
      stop(sprintf(
        "No antibiotic class column ('%s') and no antibiotic name column ('%s') found. Cannot derive classes.",
        .null_default(col_map$class_col, "antibiotic_class"), antibiotic_col
      ), call. = FALSE)
    }
  } else {
    plog <- .log_step(plog, "3a_class_col", n_in, nrow(data),
                      sprintf("Antibiotic class column '%s' already present. Standardisation skipped.", class_col))
  }

  # Ensure class_col is set
  if (!"antibiotic_class" %in% names(data) && class_col %in% names(data)) {
    data$antibiotic_class <- data[[class_col]]
    class_col <- "antibiotic_class"
  }

  plog <- .log_step(plog, "3b_class_assignment", n_in, nrow(data),
                    sprintf("Antibiotic class column resolved to '%s'.", class_col))

  # ------------------------------------------------------------------
  # Step 4: Class-level collapse
  #   Rule: any R in class => class = R   (ertapenem NA + imipenem S + meropenem R -> Carbapenems R)
  #   Counting unit: isolate_id (not patient_id)
  # ------------------------------------------------------------------
  message("\n[preprocess_for_profiles] Step 4: Collapsing to class level (any R -> class R)...")
  n_in <- nrow(data)

  # prep_collapse_class_level() hardcodes a reference to "antibiotic_normalized"
  # for the drugs_tested summary column. Ensure it exists so the call never
  # errors when the user supplied antibiotic_class directly (skipping Step 3).
  if (!"antibiotic_normalized" %in% names(data)) {
    data$antibiotic_normalized <- if (antibiotic_col %in% names(data))
      data[[antibiotic_col]]
    else
      NA_character_
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
  carry_cols    <- setdiff(carry_cols, grouping_keys)

  data_collapsed <- prep_collapse_class_level(
    data             = data,
    event_col        = isolate_col,      # isolate_id is the counting unit
    organism_col     = pathogen_col,
    class_col        = class_col,
    susceptibility_col = harmonised_ast_col,
    extra_cols       = if (length(carry_cols) > 0L) carry_cols else NULL
  )

  plog <- .log_step(plog, "4_collapse_class", n_in, nrow(data_collapsed),
                    sprintf(
                      "Class-level collapse: %d drug rows -> %d isolate-class rows. Unit: %s.",
                      n_in, nrow(data_collapsed), isolate_col
                    ))
  data <- data_collapsed

  # ------------------------------------------------------------------
  # Step 4b: Isolate deduplication
  #   Same isolate_id x antibiotic combination can appear more than once
  #   (e.g. duplicate data entry, multi-join artefact). Keep the worst
  #   phenotype per isolate x drug (R > I > S > NA).
  # ------------------------------------------------------------------
  message("\n[preprocess_for_profiles] Step 4b: Deduplicating isolate x drug rows...")
  n_in <- nrow(data)

  abx_dedup_col <- if ("antibiotic_normalized" %in% names(data)) "antibiotic_normalized"
                   else if (antibiotic_col %in% names(data))      antibiotic_col
                   else NULL

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
    plog <- .log_step(plog, "4b_dedup", n_in, nrow(data),
                      sprintf("Dedup: %d duplicate isolate-drug row(s) removed (kept worst phenotype).",
                              n_removed))
  } else {
    plog <- .log_step(plog, "4b_dedup", n_in, nrow(data),
                      "Dedup skipped: no antibiotic name column identified.", status = "warn")
  }

  # ------------------------------------------------------------------
  # Step 5: Year derivation -- always performed when date_col is present.
  #   Year is always added to the output so downstream stratification
  #   by year does not require re-running preprocessing.
  # ------------------------------------------------------------------
  message("\n[preprocess_for_profiles] Step 5: Deriving year from date column...")
  n_in <- nrow(data)

  if (date_col %in% names(data)) {
    parsed_dates <- suppressWarnings(as.Date(as.character(data[[date_col]])))
    data$year    <- as.integer(format(parsed_dates, "%Y"))
    n_missing_yr <- sum(is.na(data$year))
    plog <- .log_step(plog, "5_year_derivation", n_in, nrow(data),
                      sprintf("Year derived from '%s'. %d missing year value(s).", date_col, n_missing_yr))
  } else {
    warning(sprintf("Date column '%s' not found. 'year' column not added.", date_col),
            call. = FALSE)
    plog <- .log_step(plog, "5_year_derivation", n_in, nrow(data),
                      sprintf("Skipped: date column '%s' not present.", date_col), status = "warn")
  }

  # ------------------------------------------------------------------
  # Step 6: Panel filter
  #   For each pathogen, retain only classes listed in panel_map.
  #   Everything else is logged as excluded with reason_code.
  # ------------------------------------------------------------------
  message("\n[preprocess_for_profiles] Step 6: Applying organism-specific class panels...")
  n_in <- nrow(data)

  panel_tbl <- dplyr::bind_rows(lapply(names(panel_map), function(org) {
    tibble::tibble(
      !!pathogen_col   := org,
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
  pathogens_in_data  <- unique(data[[pathogen_col]])
  pathogens_no_panel <- setdiff(pathogens_in_data, names(panel_map))
  if (length(pathogens_no_panel) > 0L) {
    warning(sprintf(
      "%d pathogen(s) in data have no entry in panel_map and will be excluded from profiling: %s",
      length(pathogens_no_panel),
      paste(head(pathogens_no_panel, 5L), collapse = ", ")
    ), call. = FALSE)
    exclusions <- c(exclusions, list(tibble::tibble(
      !!pathogen_col  := pathogens_no_panel,
      antibiotic_class = NA_character_,
      reason_code     = "no_panel_defined",
      reason_text     = "Pathogen has no entry in panel_map; excluded from profiling."
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

  plog <- .log_step(plog, "6_panel_filter", n_in, nrow(data_filtered),
                    sprintf("Panel filter: %d -> %d rows. %d exclusion(s) logged.",
                            n_in, nrow(data_filtered), nrow(excluded_rows)))
  data <- data_filtered

  # ------------------------------------------------------------------
  # Step 7: Long -> wide with antibiotic classes as columns
  #   One row per isolate_id. Class columns: S / R / NA.
  # ------------------------------------------------------------------
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
    data             = data,
    event_col        = isolate_col,
    antibiotic_col   = "antibiotic_class",
    susceptibility_col = "class_resistance",
    prefix           = "",                    # class names become column names directly
    keep_cols        = wide_keep_cols
  )

  # Sanitise column names (prep_create_wide_ast_matrix already does this, but be explicit)
  names(data_wide) <- gsub("[^A-Za-z0-9_]", "_", names(data_wide))
  names(data_wide) <- gsub("_{2,}", "_",          names(data_wide))

  plog <- .log_step(plog, "7_wide_pivot", n_in, nrow(data_wide),
                    sprintf("Long -> wide: %d isolates x %d class column(s).",
                            nrow(data_wide),
                            ncol(data_wide) - length(wide_keep_cols) - 1L))

  # ------------------------------------------------------------------
  # Collate outputs
  # ------------------------------------------------------------------
  col_map_resolved <- col_map
  col_map_resolved$ast_col_harmonised <- harmonised_ast_col
  col_map_resolved$class_col          <- "antibiotic_class"
  col_map_resolved$outcome_col        <- eff_outcome_col

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
    data_wide        = tibble::as_tibble(data_wide),
    preprocessing_log = log_tbl,
    panel_exclusions = exclusions_tbl,
    col_map_resolved = col_map_resolved,
    detected_format  = detected_format
  )
}


# ---------------------------------------------------------------------------
# check_profile_constraints()
# ---------------------------------------------------------------------------

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
#' marg    <- compute_marginal_resistance(amr_clean)
#' co_res  <- compute_pairwise_coresistance(marg)
#' rp      <- compute_resistance_profiles(marg, co_res)
#'
#' # Using stored constraint residuals only
#' checks  <- check_profile_constraints(rp)
#'
#' # Full re-verification with original rates
#' checks  <- check_profile_constraints(
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
    marginals        = NULL,
    pairwise         = NULL,
    tolerance        = 1e-6,
    pathogen_col     = "pathogen",
    class_col        = "antibiotic_class",
    rate_col         = "marginal_resistance",
    class1_col       = "antibiotic_class_1",
    class2_col       = "antibiotic_class_2",
    pairwise_rate_col = "pairwise_resistance_prevalence"
) {
  if (!is.list(profiles_output) || length(profiles_output) == 0L) {
    stop("`profiles_output` must be the non-empty named list from compute_resistance_profiles().",
         call. = FALSE)
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

    prof_df  <- entry$profiles
    p_hat    <- prof_df$probability
    n_prof   <- length(p_hat)
    classes  <- if ("classes" %in% names(entry)) entry$classes else
      setdiff(names(prof_df), non_class_cols)

    .row <- function(type, name, target, reconstructed, source) {
      abs_res <- if (is.na(target)) NA_real_ else abs(reconstructed - target)
      pass    <- if (is.na(target)) reconstructed >= -tolerance else abs_res < tolerance
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

    # ------------------------------------------------------------------
    # Check 1: Non-negativity -- all p_hat >= 0
    # ------------------------------------------------------------------
    min_p   <- min(p_hat, na.rm = TRUE)
    all_rows <- c(all_rows, list(.row(
      "nonneg", "min_probability",
      target        = NA_real_,
      reconstructed = min_p,
      source        = "recomputed"
    )))

    # ------------------------------------------------------------------
    # Check 2: Sum-to-one
    # ------------------------------------------------------------------
    p_sum <- sum(p_hat, na.rm = TRUE)
    all_rows <- c(all_rows, list(.row(
      "sum_to_one", "sum_probability",
      target        = 1.0,
      reconstructed = p_sum,
      source        = "recomputed"
    )))

    # ------------------------------------------------------------------
    # Check 3a: Marginal constraints from stored constraint_residuals
    # ------------------------------------------------------------------
    if ("constraint_residuals" %in% names(entry)) {
      resid   <- entry$constraint_residuals
      targets <- if ("constraint_targets" %in% names(entry))
                   entry$constraint_targets
                 else
                   setNames(rep(NA_real_, length(resid)), names(resid))

      marg_nms <- names(resid)[grepl("^marg_", names(resid))]
      pair_nms <- names(resid)[grepl("^pair_", names(resid))]

      for (nm in marg_nms) {
        tgt   <- targets[[nm]]
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
        tgt   <- targets[[nm]]
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

    # ------------------------------------------------------------------
    # Check 3b: Marginal constraints recomputed from supplied marginals
    # ------------------------------------------------------------------
    if (!is.null(marginals) && pathogen_col %in% names(marginals)) {
      marg_k <- marginals[marginals[[pathogen_col]] == path, ]
      marg_k <- marg_k[marg_k[[class_col]] %in% classes, ]

      if (nrow(marg_k) > 0L) {
        bin_mat <- as.matrix(prof_df[, classes, drop = FALSE])
        storage.mode(bin_mat) <- "double"

        for (i in seq_len(nrow(marg_k))) {
          cls    <- marg_k[[class_col]][i]
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

    # ------------------------------------------------------------------
    # Check 3c: Pairwise constraints recomputed from supplied pairwise
    # ------------------------------------------------------------------
    if (!is.null(pairwise) && pathogen_col %in% names(pairwise)) {
      pair_k <- pairwise[pairwise[[pathogen_col]] == path, ]

      if (nrow(pair_k) > 0L && class1_col %in% names(pair_k) &&
          class2_col %in% names(pair_k) && pairwise_rate_col %in% names(pair_k)) {

        bin_mat <- as.matrix(prof_df[, classes, drop = FALSE])
        storage.mode(bin_mat) <- "double"

        for (i in seq_len(nrow(pair_k))) {
          c1     <- pair_k[[class1_col]][i]
          c2     <- pair_k[[class2_col]][i]
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

  # ------------------------------------------------------------------
  # Per-pathogen summary printed to console
  # ------------------------------------------------------------------
  if (nrow(result_tbl) > 0L) {
    summary_tbl <- result_tbl %>%
      dplyr::group_by(pathogen) %>%
      dplyr::summarise(
        n_checks = dplyr::n(),
        n_pass   = sum(pass, na.rm = TRUE),
        n_fail   = sum(!pass, na.rm = TRUE),
        max_abs_residual = max(abs_residual, na.rm = TRUE),
        .groups  = "drop"
      )
    message("\n[check_profile_constraints] Results (tolerance = ", tolerance, "):")
    for (i in seq_len(nrow(summary_tbl))) {
      r      <- summary_tbl[i, ]
      status <- if (r$n_fail == 0L) "[OK]" else "[!!]"
      message(sprintf(
        "  %s %-40s  %d/%d passed | max |residual| = %.2e",
        status, r$pathogen, r$n_pass, r$n_checks, r$max_abs_residual
      ))
    }
  }

  result_tbl
}


# ---------------------------------------------------------------------------
# bootstrap_profiles_convex()
# ---------------------------------------------------------------------------

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
#' marg   <- compute_marginal_resistance(amr_clean)
#' co_res <- compute_pairwise_coresistance(marg)
#'
#' boot   <- bootstrap_profiles_convex(
#'   marginals          = marg$marginal,
#'   coresistance_output = co_res,
#'   B                  = 500,
#'   seed               = 42
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
    B                   = 500L,
    seed                = 123L,
    alpha               = 0.05,
    n_cores             = 1L,
    exclude_near_zero   = TRUE,
    top_n_classes       = NULL,
    sigma_sq            = 1,
    ridge               = 1e-8,
    pathogen_col        = "organism_name",
    class_col           = "antibiotic_class",
    n_tested_col        = "n_tested",
    n_resistant_col     = "n_resistant",
    org_group_col       = "org_group"
) {
  # ------------------------------------------------------------------
  # Input validation
  # ------------------------------------------------------------------
  required_marg <- c(pathogen_col, class_col, n_tested_col, n_resistant_col)
  missing_marg  <- setdiff(required_marg, names(marginals))
  if (length(missing_marg) > 0L) {
    stop(sprintf(
      "Columns not found in `marginals`: %s",
      paste(missing_marg, collapse = ", ")
    ), call. = FALSE)
  }

  B      <- as.integer(B)
  n_cores <- as.integer(n_cores)
  if (is.na(B) || B < 1L)      stop("`B` must be a positive integer.", call. = FALSE)
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

  # ------------------------------------------------------------------
  # Point estimate on original marginals (stored as attribute)
  # ------------------------------------------------------------------
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
    marginal_output     = .build_marginal_output(marginals),
    coresistance_output = if (is.null(coresistance_output)) list() else coresistance_output,
    exclude_near_zero   = exclude_near_zero,
    top_n_classes       = top_n_classes,
    sigma_sq            = sigma_sq,
    ridge               = ridge,
    pathogen_col        = pathogen_col,
    antibiotic_class_col = class_col,
    n_cores             = 1L
  )

  if (length(point_out) == 0L) {
    stop("Point estimate returned no profiles. Check marginals and class filters.", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # Bootstrap replicates
  # ------------------------------------------------------------------
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
    n_t    <- marg_b[[n_tested_col]]
    r_obs  <- pmin(pmax(marg_b$marginal_resistance, 0), 1)  # guard [0,1]
    r_obs[is.na(r_obs)] <- 0

    n_r_b                    <- stats::rbinom(nrow(marg_b), size = n_t, prob = r_obs)
    marg_b$marginal_resistance <- n_r_b / n_t

    # Resample pairwise T and R matrices if available
    co_b <- coresistance_output
    if (!is.null(co_b) && length(co_b) > 0L) {
      co_b <- lapply(co_b, function(co_path) {
        T_mat <- co_path$T_matrix
        R_mat <- co_path$R_matrix
        prev  <- co_path$prevalence
        n_c   <- nrow(T_mat)
        R_b   <- matrix(0L, n_c, n_c, dimnames = dimnames(T_mat))
        prev_b <- matrix(NA_real_, n_c, n_c, dimnames = dimnames(T_mat))

        for (i in seq_len(n_c)) {
          for (j in i:n_c) {
            if (i == j) next
            t_ij <- T_mat[i, j]
            if (is.na(t_ij) || t_ij == 0L) next
            p_ij  <- if (!is.na(prev[i, j])) prev[i, j] else 0
            r_ij  <- stats::rbinom(1L, size = t_ij, prob = pmin(pmax(p_ij, 0), 1))
            R_b[i, j] <- R_b[j, i] <- r_ij
            pv <- r_ij / t_ij
            prev_b[i, j] <- prev_b[j, i] <- pv
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
      error = function(e) NULL   # failed replicate returns NULL
    )
  }

  if (n_cores > 1L) {
    all_reps <- parallel::mclapply(seq_len(B), .one_replicate,
                                   mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    all_reps <- lapply(seq_len(B), .one_replicate)
  }

  # Remove NULL (failed) replicates
  valid_reps <- Filter(Negate(is.null), all_reps)
  n_valid    <- length(valid_reps)

  if (n_valid == 0L) stop("All bootstrap replicates failed.", call. = FALSE)
  if (n_valid < B)   warning(sprintf("%d/%d bootstrap replicates failed and were dropped.", B - n_valid, B),
                             call. = FALSE)

  # ------------------------------------------------------------------
  # Aggregate: collect p_hat vectors per pathogen x profile
  # ------------------------------------------------------------------
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
    n_prof        <- nrow(point_prof)
    uniform_p     <- 1 / n_prof
    converged     <- apply(prob_matrix, 2, function(col) !all(abs(col - uniform_p) < 1e-10))
    n_converged   <- sum(converged, na.rm = TRUE)
    conv_matrix   <- prob_matrix[, converged, drop = FALSE]

    if (ncol(conv_matrix) == 0L) {
      warning(sprintf("'%s': no converged bootstrap replicates. Intervals will be NA.", path),
              call. = FALSE)
      lower <- rep(NA_real_, n_prof)
      upper <- rep(NA_real_, n_prof)
      pmean <- rep(NA_real_, n_prof)
      pmed  <- rep(NA_real_, n_prof)
    } else {
      lower <- apply(conv_matrix, 1, stats::quantile, probs = lo_q, na.rm = TRUE)
      upper <- apply(conv_matrix, 1, stats::quantile, probs = hi_q, na.rm = TRUE)
      pmean <- apply(conv_matrix, 1, mean, na.rm = TRUE)
      pmed  <- apply(conv_matrix, 1, stats::median, na.rm = TRUE)
    }

    result_df <- tibble::tibble(
      profile                  = point_prof$profile,
      probability_mean         = round(pmean, 6L),
      probability_median       = round(pmed,  6L),
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


# ---------------------------------------------------------------------------
# enumerate_binary_profiles()
# ---------------------------------------------------------------------------

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
  if (length(classes) > 20L)
    warning(sprintf("2^%d = %d profiles may be very slow. Consider top_n_classes.", length(classes), 2L^length(classes)), call. = FALSE)

  n         <- length(classes)
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
  label_df  <- as.data.frame(char_mat,    check.names = FALSE)
  colnames(label_df) <- classes
  binary_df <- as.data.frame(profiles_mat, check.names = FALSE)

  tibble::as_tibble(
    cbind(data.frame(profile_delta = do.call(paste0, label_df),
                     stringsAsFactors = FALSE, check.names = FALSE),
          binary_df),
    .name_repair = "minimal"
  )
}


# ---------------------------------------------------------------------------
# build_constraint_matrix()
# ---------------------------------------------------------------------------

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
  classes    <- names(r_marg)
  n          <- length(classes)

  # Binary indicator matrix (2^n x n)
  profiles_mat <- as.matrix(profiles_enum[, classes, drop = FALSE])
  storage.mode(profiles_mat) <- "double"

  # Marginal rows
  M_marg <- t(profiles_mat)                       # n x 2^n
  v_marg <- r_marg[classes]
  marg_names <- paste0("marg_", classes)

  # Pairwise rows
  pairs_mat  <- utils::combn(n, 2L)
  n_pair     <- ncol(pairs_mat)
  d1_idx     <- pairs_mat[1L, ];  d2_idx <- pairs_mat[2L, ]
  c1_names   <- classes[d1_idx];  c2_names <- classes[d2_idx]
  pair_names <- paste0("pair_", c1_names, "_", c2_names)

  M_pair <- t(
    profiles_mat[, d1_idx, drop = FALSE] *
      profiles_mat[, d2_idx, drop = FALSE]
  )                                                # n(n-1)/2 x 2^n

  r1 <- r_marg[c1_names];  r2 <- r_marg[c2_names]

  # Look up observed pairwise values
  co_vals <- rep(NA_real_, n_pair)
  if (!is.null(co_mat) && !is.null(rownames(co_mat))) {
    in_r <- c1_names %in% rownames(co_mat)
    in_c <- c2_names %in% colnames(co_mat)
    ok   <- in_r & in_c
    if (any(ok)) co_vals[ok] <- co_mat[cbind(c1_names[ok], c2_names[ok])]
  }

  # Cap pairwise to min(P(A), P(B))
  cap_vals    <- pmin(r1, r2)
  has_co      <- !is.na(co_vals)
  capped_co   <- pmin(co_vals, cap_vals)
  was_capped  <- has_co & !is.na(capped_co) & (capped_co < co_vals)
  v_pair      <- ifelse(has_co, capped_co, r1 * r2)

  fallback_pairs <- pair_names[!has_co]
  capped_pairs   <- if (any(was_capped)) {
    stats::setNames(
      paste0(round(co_vals[was_capped], 4L), "->", round(capped_co[was_capped], 4L)),
      pair_names[was_capped]
    )
  } else character(0L)

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


# ---------------------------------------------------------------------------
# validate_aggregate_inputs()
# ---------------------------------------------------------------------------

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
    pairwise              = NULL,
    pathogen_col          = "pathogen",
    class_col             = "antibiotic_class",
    rate_col              = "marginal_resistance",
    n_tested_col          = "n_tested",
    n_resistant_col       = "n_resistant",
    class1_col            = "antibiotic_class_1",
    class2_col            = "antibiotic_class_2",
    pairwise_rate_col     = "pairwise_resistance_prevalence",
    min_classes_per_pathogen = 2L
) {
  results <- list()

  # 1. Required columns
  required <- c(pathogen_col, class_col, rate_col)
  missing  <- setdiff(required, names(marginals))
  if (length(missing) > 0L) {
    results <- .add_check(results, "required_cols", "fail",
      sprintf("Required column(s) missing from marginals: %s", paste(missing, collapse = ", ")))
  } else {
    results <- .add_check(results, "required_cols", "pass",
      sprintf("All required columns present (%s).", paste(required, collapse = ", ")))
  }

  if ("fail" %in% sapply(results, `[[`, "status")) {
    stop(dplyr::bind_rows(results)$message[dplyr::bind_rows(results)$status == "fail"][[1]],
         call. = FALSE)
  }

  # 2. Prevalence bounds [0, 1]
  rates    <- marginals[[rate_col]]
  n_oob    <- sum(!is.na(rates) & (rates < 0 | rates > 1))
  if (n_oob > 0L) {
    results <- .add_check(results, "prevalence_bounds", "fail",
      sprintf("%d row(s) have %s outside [0, 1].", n_oob, rate_col), n_oob)
  } else {
    results <- .add_check(results, "prevalence_bounds", "pass",
      sprintf("All %s values in [0, 1].", rate_col))
  }

  # 3. n_tested > 0 (if column present)
  if (!is.null(n_tested_col) && n_tested_col %in% names(marginals)) {
    n_zero <- sum(!is.na(marginals[[n_tested_col]]) & marginals[[n_tested_col]] <= 0L)
    if (n_zero > 0L) {
      results <- .add_check(results, "n_tested_positive", "fail",
        sprintf("%d row(s) have %s <= 0.", n_zero, n_tested_col), n_zero)
    } else {
      results <- .add_check(results, "n_tested_positive", "pass",
        sprintf("All %s > 0.", n_tested_col))
    }

    # n_resistant <= n_tested
    if (!is.null(n_resistant_col) && n_resistant_col %in% names(marginals)) {
      n_bad <- sum(!is.na(marginals[[n_resistant_col]]) &
                   !is.na(marginals[[n_tested_col]]) &
                   marginals[[n_resistant_col]] > marginals[[n_tested_col]])
      if (n_bad > 0L) {
        results <- .add_check(results, "resistant_le_tested", "fail",
          sprintf("%d row(s) have %s > %s.", n_bad, n_resistant_col, n_tested_col), n_bad)
      } else {
        results <- .add_check(results, "resistant_le_tested", "pass",
          sprintf("%s <= %s for all rows.", n_resistant_col, n_tested_col))
      }
    }
  }

  # 4. No duplicate pathogen x class rows
  dup_key <- paste(marginals[[pathogen_col]], marginals[[class_col]], sep = "||")
  n_dup   <- sum(duplicated(dup_key))
  if (n_dup > 0L) {
    results <- .add_check(results, "no_duplicates", "fail",
      sprintf("%d duplicate pathogen x class row(s) found.", n_dup), n_dup)
  } else {
    results <- .add_check(results, "no_duplicates", "pass",
      "No duplicate pathogen x class rows.")
  }

  # 5. Minimum classes per pathogen
  class_counts <- tapply(marginals[[class_col]], marginals[[pathogen_col]], length)
  n_too_few    <- sum(class_counts < min_classes_per_pathogen)
  if (n_too_few > 0L) {
    pathogens_few <- names(class_counts)[class_counts < min_classes_per_pathogen]
    results <- .add_check(results, "min_classes", "warn",
      sprintf("%d pathogen(s) have < %d classes (profiles trivial): %s",
              n_too_few, min_classes_per_pathogen,
              paste(head(pathogens_few, 5L), collapse = ", ")),
      n_too_few)
  } else {
    results <- .add_check(results, "min_classes", "pass",
      sprintf("All pathogens have >= %d classes.", min_classes_per_pathogen))
  }

  # 6. Pairwise checks (if supplied)
  if (!is.null(pairwise)) {
    pw_req <- c(pathogen_col, class1_col, class2_col, pairwise_rate_col)
    pw_miss <- setdiff(pw_req, names(pairwise))
    if (length(pw_miss) > 0L) {
      results <- .add_check(results, "pairwise_cols", "fail",
        sprintf("Pairwise table missing column(s): %s", paste(pw_miss, collapse = ", ")))
    } else {
      # Pairwise rate in [0,1]
      pw_rates <- pairwise[[pairwise_rate_col]]
      n_pw_oob <- sum(!is.na(pw_rates) & (pw_rates < 0 | pw_rates > 1))
      if (n_pw_oob > 0L) {
        results <- .add_check(results, "pairwise_bounds", "fail",
          sprintf("%d pairwise row(s) have %s outside [0, 1].", n_pw_oob, pairwise_rate_col), n_pw_oob)
      }

      # Pairwise <= min(marginal_A, marginal_B)
      pw_joined <- pairwise %>%
        dplyr::left_join(marginals[, c(pathogen_col, class_col, rate_col)],
                         by = stats::setNames(c(pathogen_col, class_col), c(pathogen_col, class1_col))) %>%
        dplyr::rename(rate_A = !!rate_col) %>%
        dplyr::left_join(marginals[, c(pathogen_col, class_col, rate_col)],
                         by = stats::setNames(c(pathogen_col, class_col), c(pathogen_col, class2_col))) %>%
        dplyr::rename(rate_B = !!rate_col)

      n_impossible <- sum(!is.na(pw_joined[[pairwise_rate_col]]) &
                          !is.na(pw_joined$rate_A) & !is.na(pw_joined$rate_B) &
                          pw_joined[[pairwise_rate_col]] > pmin(pw_joined$rate_A, pw_joined$rate_B))
      if (n_impossible > 0L) {
        results <- .add_check(results, "pairwise_le_marginals", "warn",
          sprintf("%d pairwise value(s) exceed min(marginal_A, marginal_B) and will be capped.", n_impossible), n_impossible)
      } else {
        results <- .add_check(results, "pairwise_le_marginals", "pass",
          "All pairwise values <= min(marginal_A, marginal_B).")
      }
    }
  }

  result_tbl <- dplyr::bind_rows(results)
  n_fail <- sum(result_tbl$status == "fail")
  n_warn <- sum(result_tbl$status == "warn")
  n_pass <- sum(result_tbl$status == "pass")

  message(sprintf("[validate_aggregate_inputs] %d passed | %d warnings | %d failed",
                  n_pass, n_warn, n_fail))
  if (n_warn > 0L)
    for (msg in result_tbl$message[result_tbl$status == "warn"]) message("  [!] ", msg)

  if (n_fail > 0L) {
    stop(sprintf("%d check(s) failed:\n%s", n_fail,
                 paste(sprintf("  - %s", result_tbl$message[result_tbl$status == "fail"]),
                       collapse = "\n")), call. = FALSE)
  }

  invisible(result_tbl)
}


# ---------------------------------------------------------------------------
# compute_marginals_from_data()
# ---------------------------------------------------------------------------

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
    stratify_by        = NULL,
    outcome_col        = NULL,
    min_n_tested       = 30L,
    external_marginals = NULL,
    ext_col_map        = list(
      pathogen_col  = "pathogen",
      class_col     = "antibiotic_class",
      geography_col = "geography",
      year_col      = "year",
      rate_col      = "resistance_prevalence"
    )
) {
  pathogen_col  <- .null_default(col_map$pathogen_col,  "pathogen")
  geography_col <- .null_default(col_map$geography_col, "state")
  isolate_col   <- .null_default(col_map$isolate_col,   "isolate_id")

  # Identify antibiotic-class columns: all columns that are not known metadata
  meta_cols <- unique(c(
    isolate_col, pathogen_col,
    col_map$patient_col, col_map$date_col, geography_col,
    col_map$specimen_col, col_map$age_col, col_map$dob_col,
    col_map$location_col, outcome_col, "year"
  ))
  meta_cols  <- meta_cols[!is.null(meta_cols) & meta_cols %in% names(data_wide)]
  class_cols <- setdiff(names(data_wide), meta_cols)
  class_cols <- class_cols[class_cols %in% unlist(panel_map)]   # keep only panel classes

  if (length(class_cols) == 0L)
    stop("No antibiotic-class columns found in data_wide matching panel_map. Check panel_map class names.", call. = FALSE)

  # Build stratification grouping columns
  strat_cols <- character(0)
  if ("geography" %in% stratify_by && geography_col %in% names(data_wide))
    strat_cols <- c(strat_cols, geography_col)
  if ("year" %in% stratify_by && "year" %in% names(data_wide))
    strat_cols <- c(strat_cols, "year")
  if (!is.null(outcome_col) && outcome_col %in% names(data_wide))
    strat_cols <- c(strat_cols, outcome_col)

  group_vars <- c(strat_cols, pathogen_col)

  # Pivot class columns to long, then compute marginals
  data_long <- data_wide %>%
    dplyr::select(dplyr::all_of(c(isolate_col, group_vars, class_cols))) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(class_cols),
      names_to  = "antibiotic_class",
      values_to = "class_result"
    ) %>%
    dplyr::filter(!is.na(class_result))   # only tested isolate-class pairs

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
  n_before  <- nrow(marginals_tbl)
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

    ext_path_col <- .null_default(ext_col_map$pathogen_col,  "pathogen")
    ext_cls_col  <- .null_default(ext_col_map$class_col,     "antibiotic_class")
    ext_rate_col <- .null_default(ext_col_map$rate_col,      "resistance_prevalence")
    ext_geo_col  <- .null_default(ext_col_map$geography_col, "geography")
    ext_yr_col   <- .null_default(ext_col_map$year_col,      "year")

    join_keys <- c(pathogen_col = ext_path_col, antibiotic_class = ext_cls_col)
    if ("geography" %in% stratify_by && ext_geo_col %in% names(external_marginals))
      join_keys <- c(join_keys, stats::setNames(ext_geo_col, geography_col))
    if ("year" %in% stratify_by && ext_yr_col %in% names(external_marginals))
      join_keys <- c(join_keys, stats::setNames(ext_yr_col, "year"))

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


# ---------------------------------------------------------------------------
# compute_pairwise_from_data()
# ---------------------------------------------------------------------------

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
#'   \item{\code{pairwise_prevalence}}{Back-calculated P(A∩B), capped and floored.}
#'   \item{\code{method}}{\code{"pearson_back_calc"} or \code{"independence_fallback"}.}
#' }
#'
#' @export
compute_pairwise_from_data <- function(
    data_wide,
    marginals,
    col_map,
    panel_map,
    stratify_by  = NULL,
    outcome_col  = NULL,
    min_co_tested = 10L
) {
  pathogen_col  <- .null_default(col_map$pathogen_col,  "pathogen")
  geography_col <- .null_default(col_map$geography_col, "state")
  isolate_col   <- .null_default(col_map$isolate_col,   "isolate_id")

  meta_cols <- unique(c(
    isolate_col, pathogen_col,
    col_map$patient_col, col_map$date_col, geography_col,
    col_map$specimen_col, col_map$age_col, col_map$dob_col,
    col_map$location_col, outcome_col, "year"
  ))
  meta_cols  <- meta_cols[!is.null(meta_cols) & meta_cols %in% names(data_wide)]
  class_cols <- setdiff(names(data_wide), meta_cols)
  class_cols <- class_cols[class_cols %in% unlist(panel_map)]

  strat_cols <- character(0)
  if ("geography" %in% stratify_by && geography_col %in% names(data_wide))
    strat_cols <- c(strat_cols, geography_col)
  if ("year" %in% stratify_by && "year" %in% names(data_wide))
    strat_cols <- c(strat_cols, "year")
  if (!is.null(outcome_col) && outcome_col %in% names(data_wide))
    strat_cols <- c(strat_cols, outcome_col)

  group_key_cols <- c(strat_cols, pathogen_col)

  # Get unique strata x pathogen combinations
  strata_df <- data_wide %>%
    dplyr::select(dplyr::all_of(c(group_key_cols))) %>%
    dplyr::distinct()

  all_pairs <- list()

  for (i in seq_len(nrow(strata_df))) {
    row_key  <- strata_df[i, , drop = FALSE]
    path     <- row_key[[pathogen_col]]
    cls_panel <- panel_map[[path]]
    if (is.null(cls_panel)) next
    cls_avail <- intersect(cls_panel, class_cols)
    if (length(cls_avail) < 2L) next

    # Filter data_wide to this stratum
    stratum_filter <- rep(TRUE, nrow(data_wide))
    for (sc in strat_cols)
      stratum_filter <- stratum_filter & (data_wide[[sc]] == row_key[[sc]])
    stratum_filter <- stratum_filter & (data_wide[[pathogen_col]] == path)

    sub_wide <- data_wide[stratum_filter, cls_avail, drop = FALSE]

    # Get marginals for this stratum x pathogen
    marg_filter <- marginals[[pathogen_col]] == path
    for (sc in strat_cols)
      marg_filter <- marg_filter & (marginals[[sc]] == row_key[[sc]])
    marg_k <- marginals[marg_filter, , drop = FALSE]
    p_lookup <- stats::setNames(marg_k$marginal_resistance, marg_k$antibiotic_class)

    # Enumerate class pairs
    pairs_mat <- utils::combn(cls_avail, 2L)
    for (j in seq_len(ncol(pairs_mat))) {
      c1 <- pairs_mat[1L, j];  c2 <- pairs_mat[2L, j]
      if (!c1 %in% names(p_lookup) || !c2 %in% names(p_lookup)) next

      PA <- p_lookup[[c1]];  PB <- p_lookup[[c2]]
      if (is.na(PA) || is.na(PB)) next

      # Co-tested isolates: non-NA in both columns
      v1   <- sub_wide[[c1]];  v2 <- sub_wide[[c2]]
      both <- !is.na(v1) & !is.na(v2)
      n_co <- sum(both)

      if (n_co < min_co_tested) {
        pairwise_val <- PA * PB
        method       <- "independence_fallback"
        rho          <- NA_real_
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
          method       <- "independence_fallback"
        } else {
          # GBD back-calculation formula
          var_A <- PA * (1 - PA)
          var_B <- PB * (1 - PB)
          p_ab_raw <- PA * PB + rho * sqrt(var_A * var_B)
          # Cap to min(PA, PB) and floor to independence product
          p_ab_capped <- pmin(p_ab_raw, min(PA, PB))
          pairwise_val <- max(p_ab_capped, PA * PB * 0)  # floor at 0, not independence
          pairwise_val <- max(pairwise_val, 0)
          method       <- "pearson_back_calc"
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
  n_pearson  <- sum(result$method == "pearson_back_calc")
  n_fallback <- sum(result$method == "independence_fallback")
  message(sprintf(
    "[compute_pairwise_from_data] %d pairs: %d Pearson back-calc | %d independence fallback.",
    nrow(result), n_pearson, n_fallback
  ))
  result
}


# ---------------------------------------------------------------------------
# estimate_profiles_convex()
# ---------------------------------------------------------------------------

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
#' marg   <- compute_marginals_from_data(result$data_wide, result$col_map_resolved, panel_map)
#' pw     <- compute_pairwise_from_data(result$data_wide, marg, result$col_map_resolved, panel_map)
#' out    <- estimate_profiles_convex(marg, pw, panel_map)
#' }
estimate_profiles_convex <- function(
    marginals,
    pairwise          = NULL,
    panel_map         = NULL,
    exclude_near_zero = TRUE,
    zero_threshold    = 0,
    top_n_classes     = NULL,
    lambda            = 1e-8,
    sigma_sq          = 1,
    solver            = c("osqp", "quadprog"),
    n_cores           = 1L,
    pathogen_col      = "pathogen",
    class_col         = "antibiotic_class",
    rate_col          = "marginal_resistance",
    n_tested_col      = "n_tested"
) {
  solver  <- match.arg(solver)
  n_cores <- as.integer(n_cores)

  has_osqp     <- requireNamespace("osqp",     quietly = TRUE) &&
                  requireNamespace("Matrix",   quietly = TRUE)
  has_quadprog <- requireNamespace("quadprog", quietly = TRUE)
  if (!has_osqp && !has_quadprog)
    stop("Install 'osqp' + 'Matrix' (recommended) or 'quadprog'.", call. = FALSE)
  use_osqp <- has_osqp && solver != "quadprog"

  req_cols <- c(pathogen_col, class_col, rate_col)
  miss <- setdiff(req_cols, names(marginals))
  if (length(miss) > 0L)
    stop(sprintf("Column(s) not found in `marginals`: %s", paste(miss, collapse = ", ")), call. = FALSE)

  # Detect stratum columns automatically
  known_non_strat <- c(pathogen_col, class_col, rate_col, n_tested_col,
                       "n_resistant", "marginal_source", "org_group")
  strat_cols <- setdiff(names(marginals), known_non_strat)

  # Build unique (stratum x pathogen) combinations
  key_cols  <- c(strat_cols, pathogen_col)
  strata_df <- marginals %>%
    dplyr::select(dplyr::all_of(key_cols)) %>%
    dplyr::distinct()

  message(sprintf(
    "[estimate_profiles_convex] %d stratum x pathogen combination(s) | solver: %s | n_cores: %d",
    nrow(strata_df), if (use_osqp) "osqp" else "quadprog", n_cores
  ))

  # Build a look-up structure for pairwise prevalence
  .get_co_mat <- function(path, row_key) {
    if (is.null(pairwise) || nrow(pairwise) == 0L) return(NULL)
    pw_f <- pairwise[[pathogen_col]] == path
    for (sc in strat_cols)
      pw_f <- pw_f & (pairwise[[sc]] == row_key[[sc]])
    pw_k <- pairwise[pw_f, , drop = FALSE]
    if (nrow(pw_k) == 0L) return(NULL)

    all_cls <- unique(c(pw_k$antibiotic_class_1, pw_k$antibiotic_class_2))
    mat     <- matrix(NA_real_, length(all_cls), length(all_cls),
                      dimnames = list(all_cls, all_cls))
    for (ri in seq_len(nrow(pw_k))) {
      c1 <- pw_k$antibiotic_class_1[ri];  c2 <- pw_k$antibiotic_class_2[ri]
      pv <- pw_k$pairwise_prevalence[ri]
      mat[c1, c2] <- pv;  mat[c2, c1] <- pv
    }
    mat
  }

  .solve_one <- function(idx) {
    row_key  <- strata_df[idx, , drop = FALSE]
    path     <- row_key[[pathogen_col]]

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
      cls   <- sort(n_ord[[class_col]][seq_len(top_n_classes)])
    }

    if (length(cls) < 2L) {
      message(sprintf("  '%s': fewer than 2 classes after filtering -- skipped.", path))
      return(NULL)
    }

    r_marg   <- stats::setNames(marg_k[[rate_col]][match(cls, marg_k[[class_col]])], cls)
    co_mat_k <- .get_co_mat(path, row_key)

    # Enumerate profiles and build constraint matrix
    profiles_enum <- enumerate_binary_profiles(cls)
    cm            <- build_constraint_matrix(profiles_enum, r_marg, co_mat_k)
    M <- cm$M;  v <- cm$v
    n_profiles <- nrow(profiles_enum)

    # Identifiability: rank(M) vs 2^n
    rank_M            <- qr(M)$rank
    identif_flag      <- rank_M < n_profiles

    # QP
    coef  <- 2.0 / sigma_sq
    H_mat <- coef * crossprod(M)
    diag(H_mat) <- diag(H_mat) + lambda
    d_qp  <- coef * drop(crossprod(M, v))

    converged <- TRUE
    p_hat <- tryCatch({
      if (use_osqp) {
        A_sp <- rbind(
          Matrix::Matrix(1.0, nrow = 1L, ncol = n_profiles, sparse = TRUE),
          Matrix::Diagonal(n_profiles)
        )
        prob <- osqp::osqp(
          P    = Matrix::forceSymmetric(Matrix::Matrix(H_mat)),
          q    = -d_qp,
          A    = A_sp,
          l    = c(1.0, rep(0.0, n_profiles)),
          u    = c(1.0, rep(Inf, n_profiles)),
          pars = osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-8,
                                    eps_rel = 1e-8, max_iter = 10000L, polish = TRUE)
        )
        res <- prob$solve()
        if (!(res$info$status %in% c("solved", "solved_inaccurate"))) stop(res$info$status)
        pmax(res$x, 0.0)
      } else {
        Amat <- cbind(rep(1.0, n_profiles), diag(n_profiles))
        bvec <- c(1.0, rep(0.0, n_profiles))
        pmax(quadprog::solve.QP(H_mat, d_qp, Amat, bvec, meq = 1L)$solution, 0.0)
      }
    }, error = function(e) {
      converged <<- FALSE
      rep(1.0 / n_profiles, n_profiles)
    })
    p_hat <- p_hat / sum(p_hat)

    residuals    <- drop(M %*% p_hat) - v
    max_abs_res  <- max(abs(residuals))
    notes_parts  <- character(0)
    if (length(cm$capped_pairs) > 0L)
      notes_parts <- c(notes_parts, paste0("capped: ", paste(names(cm$capped_pairs), collapse = ", ")))
    if (length(cm$fallback_pairs) > 0L)
      notes_parts <- c(notes_parts, paste0("indep_fallback: ", paste(cm$fallback_pairs, collapse = ", ")))
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
        pathogen         = path,
        profile_set_type = "aggregate_convex",
        profile_class_set = paste(cls, collapse = "|"),
        estimator        = "convex"
      )

    # Attach stratum columns
    for (sc in strat_cols) out_df[[sc]] <- row_key[[sc]]

    # Reorder columns: strat | pathogen | profile meta | probability | flags
    front_cols <- c(strat_cols, "pathogen", "profile_set_type", "profile_class_set",
                    "profile_delta", "profile_probability", "estimator",
                    "convergence_flag", "identifiability_flag", "max_abs_residual", "notes")
    class_indicator_cols <- setdiff(names(out_df), front_cols)
    out_df <- out_df[, c(front_cols, class_indicator_cols), drop = FALSE]
    tibble::as_tibble(out_df)
  }

  if (n_cores > 1L) {
    results_list <- parallel::mclapply(seq_len(nrow(strata_df)), .solve_one,
                                       mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    results_list <- lapply(seq_len(nrow(strata_df)), .solve_one)
  }

  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0L)
    stop("No profiles estimated. Check marginals, panel_map, and class filters.", call. = FALSE)

  final_tbl <- dplyr::bind_rows(results_list)
  message(sprintf(
    "\n[estimate_profiles_convex] Done. %d row(s) | %d pathogen(s) | %d stratum combination(s).",
    nrow(final_tbl),
    dplyr::n_distinct(final_tbl$pathogen),
    if (length(strat_cols) > 0L) dplyr::n_distinct(final_tbl[, strat_cols, drop = FALSE]) else 1L
  ))
  final_tbl
}
