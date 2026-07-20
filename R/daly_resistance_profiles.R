# daly_resistance_profiles.R
#
# Estimates antimicrobial resistance profiles for two complementary data sources.
# All profile probability distributions produced by either pathway feed directly
# into the DALY burden pipeline (YLL and YLD calculations).
#
# Pathway 1 -- Convex optimisation (aggregate surveillance or GBD-style data)
#   Accepts marginal resistance prevalences and optional pairwise co-resistance
#   rates per pathogen, enumerates all 2^n binary resistance profiles, and
#   recovers a valid probability distribution over those profiles by solving a
#   simplex-constrained weighted least-squares quadratic programme.
#   Supports bootstrapped uncertainty intervals, stratification by geography
#   and year, and integration of externally modelled marginals.
#
#   Key functions:
#     validate_profile_inputs()       input checks, format detection
#     preprocess_for_profiles()       full AST preprocessing pipeline
#     compute_marginals_from_data()   marginal resistance from wide isolate data
#     compute_pairwise_from_data()    pairwise co-resistance via Pearson back-calculation
#     validate_aggregate_inputs()     validates pre-computed aggregate inputs
#     enumerate_binary_profiles()     enumerates 2^n binary profiles
#     build_constraint_matrix()       constructs QP constraint matrix and target vector
#     estimate_profiles_convex()      solves the QP and returns profile probabilities
#     bootstrap_profiles_convex()     quantifies uncertainty via binomial resampling
#     check_profile_constraints()     verifies non-negativity, sum-to-one, and residuals
#     compute_marginal_resistance()   class-level marginals with any-R collapse rule
#     compute_pairwise_coresistance() pairwise co-resistance matrices per pathogen
#     compute_resistance_profiles()   QP profile estimation from isolate-level inputs
#     select_resistance_class()       selects one resistance class per event for attribution
#
# Pathway 2 -- Bayesian hierarchical modelling (facility-level AST data)
#   Fits a multivariate probit model to event-level AST records, accounting for
#   partial and selective testing, hospital-level heterogeneity, patient-admission
#   clustering, and correlated resistance outcomes across antibiotic classes.
#   Produces hospital-specific R_ALL and R_NF profile distributions for YLL
#   and YLD calculations respectively.
#
#   Key functions:
#     [Pathway 2 functions follow below]


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

  if (!is.null(facility_col) || !is.null(facility_name)) {
    if (is.null(facility_col) || is.null(facility_name)) {
      stop("Both facility_col and facility_name must be provided together, or both NULL.")
    }
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

  if (!is.null(outcome_col) || !is.null(outcome_value)) {
    if (is.null(outcome_col) || is.null(outcome_value)) {
      stop("Both outcome_col and outcome_value must be provided together, or both NULL.")
    }
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

  if (!is.null(facility_col) || !is.null(facility_name)) {
    if (is.null(facility_col) || is.null(facility_name)) {
      stop("Both facility_col and facility_name must be provided together, or both NULL.")
    }
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

  if (!is.null(outcome_col) || !is.null(outcome_value)) {
    if (is.null(outcome_col) || is.null(outcome_value)) {
      stop("Both outcome_col and outcome_value must be provided together, or both NULL.")
    }
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

  if (!is.null(facility_col) || !is.null(facility_name)) {
    if (is.null(facility_col) || is.null(facility_name)) {
      stop("Both facility_col and facility_name must be provided together, or both NULL.")
    }
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

  if (!is.null(outcome_col) || !is.null(outcome_value)) {
    if (is.null(outcome_col) || is.null(outcome_value)) {
      stop("Both outcome_col and outcome_value must be provided together, or both NULL.")
    }
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
    profiles_mat <- matrix(
      as.integer(
        outer(0L:(n_profiles - 1L), 2L^(0L:(n - 1L)), bitwAnd) > 0L
      ),
      nrow = n_profiles, ncol = n,
      dimnames = list(NULL, classes)
    )

    char_mat <- matrix("S", nrow = n_profiles, ncol = n)
    char_mat[profiles_mat == 1L] <- "R"
    profile_labels <- do.call(paste0, as.data.frame(char_mat))

    # -- Constraint matrix M (m x 2^n) and target vector v -----------------
    M_marg <- t(profiles_mat)
    v_marg <- r_marg

    pairs_mat <- utils::combn(n, 2L)
    n_pair <- ncol(pairs_mat)
    d1_idx <- pairs_mat[1L, ]
    d2_idx <- pairs_mat[2L, ]
    c1_names <- classes[d1_idx]
    c2_names <- classes[d2_idx]
    pair_names <- paste0(c1_names, "_", c2_names)

    M_pair <- t(
      profiles_mat[, d1_idx, drop = FALSE] *
        profiles_mat[, d2_idx, drop = FALSE]
    )

    r1 <- r_marg[c1_names]
    r2 <- r_marg[c2_names]

    if (!is.null(co_mat) &&
      !is.null(rownames(co_mat)) && !is.null(colnames(co_mat))) {
      in_rows <- c1_names %in% rownames(co_mat)
      in_cols <- c2_names %in% colnames(co_mat)
      can_look <- in_rows & in_cols
      co_vals <- rep(NA_real_, n_pair)
      if (any(can_look)) {
        co_vals[can_look] <- co_mat[cbind(
          c1_names[can_look],
          c2_names[can_look]
        )]
      }
    } else {
      co_vals <- rep(NA_real_, n_pair)
    }

    cap_vals <- pmin(r1, r2)
    has_co <- !is.na(co_vals)
    capped_co <- pmin(co_vals, cap_vals)
    was_capped <- has_co & !is.na(capped_co) & (capped_co < co_vals)
    v_pair <- ifelse(has_co, capped_co, r1 * r2)

    if (any(was_capped)) {
      message(sprintf(
        "'%s': %d pairwise value(s) capped to min(marginal): %s",
        path, sum(was_capped),
        paste(
          sprintf(
            "(%s,%s) %.4f->%.4f",
            c1_names[was_capped], c2_names[was_capped],
            co_vals[was_capped], capped_co[was_capped]
          ),
          collapse = "; "
        )
      ))
    }
    if (any(!has_co)) {
      message(sprintf(
        "'%s': %d pair(s) used independence fallback: %s",
        path, sum(!has_co),
        paste(sprintf("(%s,%s)", c1_names[!has_co], c2_names[!has_co]),
          collapse = "; "
        )
      ))
    }

    M <- rbind(M_marg, M_pair)
    v <- c(v_marg, v_pair)
    storage.mode(M) <- "double"

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
    names(residuals) <- c(
      paste0("marg_", classes),
      paste0("pair_", pair_names)
    )

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
        constraint_targets    = setNames(v, names(residuals)),
        constraint_names      = names(residuals)
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


# ---------------------------------------------------------------------------
# Resistance class selection (moved from prep_ast_and_syndrome.R)
# These are analysis functions that inform DALY burden attribution - they
# do not belong in the preprocessing layer.
# ---------------------------------------------------------------------------

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
                                    event_col          = "event_id",
                                    class_col          = "antibiotic_class",
                                    susceptibility_col = "antibiotic_value",
                                    rr_col             = "rr_value",
                                    hierarchy          = NULL,
                                    filter_resistant   = TRUE) {
  required_cols <- c(event_col, class_col, susceptibility_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

  if (is.null(hierarchy)) hierarchy <- get_beta_lactam_hierarchy()

  use_rr <- rr_col %in% names(data)
  if (!use_rr)
    message(sprintf("RR column '%s' not found. Using hierarchy only for selection.", rr_col))

  n_before        <- nrow(data)
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  message(sprintf("Selecting resistance classes using hierarchy%s...",
                  ifelse(use_rr, " + RR", "")))

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
                                  rr_col    = NULL,
                                  hierarchy) {
  # Accept either a named vector (name = class, value = rank) or an unnamed
  # character vector (position = rank) -- get_beta_lactam_hierarchy() returns
  # the latter.
  hier_classes <- if (!is.null(names(hierarchy))) names(hierarchy) else as.character(hierarchy)
  hierarchy_df <- data.frame(
    class          = hier_classes,
    hierarchy_rank = seq_along(hierarchy),
    stringsAsFactors = FALSE
  )
  names(hierarchy_df)[1] <- class_col

  data     <- data %>% dplyr::left_join(hierarchy_df, by = class_col)
  max_rank <- max(hierarchy_df$hierarchy_rank, na.rm = TRUE)
  data     <- data %>%
    dplyr::mutate(hierarchy_rank = dplyr::coalesce(hierarchy_rank, max_rank + 1L))

  if (!is.null(rr_col) && rr_col %in% names(data)) {
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(hierarchy_rank, dplyr::desc(!!rlang::sym(rr_col)),
                     !!rlang::sym(class_col)) %>%
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

  if (nrow(multi_class_events) > 0)
    message(sprintf("Applied selection to %d events with multiple resistant classes.",
                    nrow(multi_class_events)))

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




# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Check pairwise co-testing overlap per hospital x class-pair
#' @keywords internal
.check_pairwise_cotesting <- function(event_class_data, class_cols, upper_re_col,
                                       min_pairwise_cotested = 30L) {
  rows <- list()
  hospitals <- sort(unique(event_class_data[[upper_re_col]]))
  for (h in hospitals) {
    sub <- event_class_data[event_class_data[[upper_re_col]] == h, class_cols,
                             drop = FALSE]
    pairs <- utils::combn(class_cols, 2L, simplify = FALSE)
    for (pr in pairs) {
      c1 <- pr[1L]; c2 <- pr[2L]
      n_cotested <- sum(!is.na(sub[[c1]]) & !is.na(sub[[c2]]))
      rows[[length(rows) + 1L]] <- tibble::tibble(
        !!upper_re_col   := h,
        class_1           = c1,
        class_2           = c2,
        n_cotested        = n_cotested,
        sufficient        = n_cotested >= min_pairwise_cotested
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
    min_tested      = 30L,
    min_resistant   = 5L,
    min_susceptible = 5L
) {
  rows <- list()
  pathogens <- sort(unique(event_class_data[[pathogen_col]]))
  hospitals <- sort(unique(event_class_data[[upper_re_col]]))

  for (h in hospitals) {
    for (g in pathogens) {
      sub <- event_class_data[
        event_class_data[[upper_re_col]] == h &
        event_class_data[[pathogen_col]]  == g, , drop = FALSE
      ]
      if (nrow(sub) == 0L) next
      for (cc in class_cols) {
        vals <- sub[[cc]]
        n_tested <- sum(!is.na(vals))
        n_res    <- sum(!is.na(vals) & vals == 1L)
        n_sus    <- sum(!is.na(vals) & vals == 0L)
        eligible <- (n_tested >= min_tested) &
                    (n_res    >= min_resistant) &
                    (n_sus    >= min_susceptible)
        rows[[length(rows) + 1L]] <- tibble::tibble(
          !!upper_re_col := h,
          !!pathogen_col  := g,
          antibiotic_class = cc,
          n_tested         = n_tested,
          n_resistant      = n_res,
          n_susceptible    = n_sus,
          eligible         = eligible
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}


# ---------------------------------------------------------------------------
# Internal: Stan model code -- three variants by number of RE levels
# ---------------------------------------------------------------------------
# Priors are passed as data so the model compiles once and prior values
# can be changed freely between runs without recompilation.
#
# KEY CHANGE vs prior skeleton:
#   - eta[N_events, D] in parameters: per-event latent noise whose prior is
#     N(0,I). Because z_{ed} = mu_{ed} + hospital_effect[d,h] + (L_Omega * eta_e)[d],
#     L_Omega participates in the likelihood. Previously L_Omega only had a
#     prior and was never updated by data.
#   - Hospital (and lower-level) effects use diag_pre_multiply(tau, L_corr) * z_raw
#     so they are correlated across antibiotic classes.

.amr_probit_stan_1re <- function() {
  r"(
// Multivariate probit with data augmentation -- one grouping level (hospital)
//
// Correct implementation: z_aug[e] ~ MVN(mu_e, Omega).
// Sign constraints on observed outcomes are imposed via an exp-transformation
// of the unconstrained z_free parameters. The Jacobian correction accounts for
// the change-of-variables from z_free to z_aug. There is NO additional
// Bernoulli/probit residual -- the MVN supplies the entire residual covariance.
data {
  int<lower=1> N;                               // observed (event x class) pairs
  int<lower=1> N_events;                        // unique events
  int<lower=2> D;                               // antibiotic classes
  int<lower=1> H;                               // hospitals
  int<lower=1> K;                               // FE columns (including intercept)

  matrix[N_events, K] X_event;                  // event-level FE design matrix
  array[N_events] int<lower=1,upper=H> h_ev;    // hospital index per event

  // Wide AST arrays for data augmentation
  array[N_events, D] int<lower=0,upper=1> obs_mask; // 1 = tested, 0 = not tested
  array[N_events, D] int<lower=0,upper=1> y_mat;    // outcome if tested, else 0

  // Long-format indices (for Jacobian only)
  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;

  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[D, H] z_hospital;
  vector<lower=0>[D] tau_hospital;
  cholesky_factor_corr[D] L_corr_hospital;
  cholesky_factor_corr[D] L_Omega;
  matrix[N_events, D] z_free; // unconstrained; sign-constrained in TP block
}
transformed parameters {
  matrix[D, H] hospital_effect =
    diag_pre_multiply(tau_hospital, L_corr_hospital) * z_hospital;

  // Sign-constrained latent utilities
  matrix[N_events, D] z_aug;
  for (e in 1:N_events) {
    for (d in 1:D) {
      if (obs_mask[e, d] == 1) {
        z_aug[e, d] = (y_mat[e, d] == 1) ? exp(z_free[e, d]) : -exp(z_free[e, d]);
      } else {
        z_aug[e, d] = z_free[e, d];
      }
    }
  }

  // Event-level linear predictor (N_events x D)
  matrix[N_events, D] mu_mat = X_event * beta;
  for (e in 1:N_events) {
    for (d in 1:D) {
      mu_mat[e, d] += hospital_effect[d, h_ev[e]];
    }
  }
}
model {
  to_vector(beta)         ~ normal(0, prior_beta_sd);
  tau_hospital            ~ normal(0, prior_tau_sd);
  to_vector(z_hospital)   ~ std_normal();
  L_corr_hospital         ~ lkj_corr_cholesky(lkj_eta);
  L_Omega                 ~ lkj_corr_cholesky(lkj_eta);

  // Multivariate-probit data augmentation: z_aug[e] ~ MVN(mu_mat[e], Omega)
  for (e in 1:N_events) {
    target += multi_normal_cholesky_lpdf(z_aug[e]' | mu_mat[e]', L_Omega);
  }

  // Jacobian for sign-constrained observations: log|dz_aug/dz_free| = z_free
  for (n in 1:N) {
    target += z_free[ev_idx[n], d_idx[n]];
  }
}
generated quantities {
  corr_matrix[D] Omega      = multiply_lower_tri_self_transpose(L_Omega);
  corr_matrix[D] R_hospital = multiply_lower_tri_self_transpose(L_corr_hospital);
}
)"
}

.amr_probit_stan_2re <- function() {
  r"(
// Multivariate probit with data augmentation -- two grouping levels
// (hospital + patient or admission)
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> Pt;
  int<lower=1> K;

  matrix[N_events, K] X_event;
  array[N_events] int<lower=1,upper=H>  h_ev;
  array[N_events] int<lower=1,upper=Pt> p_ev;

  array[N_events, D] int<lower=0,upper=1> obs_mask;
  array[N_events, D] int<lower=0,upper=1> y_mat;

  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;

  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[D, H] z_hospital;
  vector<lower=0>[D] tau_hospital;
  cholesky_factor_corr[D] L_corr_hospital;
  matrix[D, Pt] z_patient;
  vector<lower=0>[D] tau_patient;
  cholesky_factor_corr[D] L_corr_patient;
  cholesky_factor_corr[D] L_Omega;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, H]  hospital_effect =
    diag_pre_multiply(tau_hospital, L_corr_hospital) * z_hospital;
  matrix[D, Pt] patient_effect  =
    diag_pre_multiply(tau_patient,  L_corr_patient)  * z_patient;

  matrix[N_events, D] z_aug;
  for (e in 1:N_events) {
    for (d in 1:D) {
      if (obs_mask[e, d] == 1) {
        z_aug[e, d] = (y_mat[e, d] == 1) ? exp(z_free[e, d]) : -exp(z_free[e, d]);
      } else {
        z_aug[e, d] = z_free[e, d];
      }
    }
  }

  matrix[N_events, D] mu_mat = X_event * beta;
  for (e in 1:N_events) {
    for (d in 1:D) {
      mu_mat[e, d] += hospital_effect[d, h_ev[e]] + patient_effect[d, p_ev[e]];
    }
  }
}
model {
  to_vector(beta)        ~ normal(0, prior_beta_sd);
  tau_hospital           ~ normal(0, prior_tau_sd);
  tau_patient            ~ normal(0, prior_tau_sd);
  to_vector(z_hospital)  ~ std_normal();
  to_vector(z_patient)   ~ std_normal();
  L_corr_hospital        ~ lkj_corr_cholesky(lkj_eta);
  L_corr_patient         ~ lkj_corr_cholesky(lkj_eta);
  L_Omega                ~ lkj_corr_cholesky(lkj_eta);

  for (e in 1:N_events) {
    target += multi_normal_cholesky_lpdf(z_aug[e]' | mu_mat[e]', L_Omega);
  }
  for (n in 1:N) {
    target += z_free[ev_idx[n], d_idx[n]];
  }
}
generated quantities {
  corr_matrix[D] Omega      = multiply_lower_tri_self_transpose(L_Omega);
  corr_matrix[D] R_hospital = multiply_lower_tri_self_transpose(L_corr_hospital);
  corr_matrix[D] R_patient  = multiply_lower_tri_self_transpose(L_corr_patient);
}
)"
}

.amr_probit_stan_3re <- function() {
  r"(
// Multivariate probit with data augmentation -- three grouping levels
// (hospital + patient + admission)
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> Pt;
  int<lower=1> Adm;
  int<lower=1> K;

  matrix[N_events, K] X_event;
  array[N_events] int<lower=1,upper=H>   h_ev;
  array[N_events] int<lower=1,upper=Pt>  p_ev;
  array[N_events] int<lower=1,upper=Adm> a_ev;

  array[N_events, D] int<lower=0,upper=1> obs_mask;
  array[N_events, D] int<lower=0,upper=1> y_mat;

  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;

  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[D, H]   z_hospital;
  vector<lower=0>[D] tau_hospital;
  cholesky_factor_corr[D] L_corr_hospital;
  matrix[D, Pt]  z_patient;
  vector<lower=0>[D] tau_patient;
  cholesky_factor_corr[D] L_corr_patient;
  matrix[D, Adm] z_admission;
  vector<lower=0>[D] tau_admission;
  cholesky_factor_corr[D] L_corr_admission;
  cholesky_factor_corr[D] L_Omega;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, H]   hospital_effect  =
    diag_pre_multiply(tau_hospital,  L_corr_hospital)  * z_hospital;
  matrix[D, Pt]  patient_effect   =
    diag_pre_multiply(tau_patient,   L_corr_patient)   * z_patient;
  matrix[D, Adm] admission_effect =
    diag_pre_multiply(tau_admission, L_corr_admission) * z_admission;

  matrix[N_events, D] z_aug;
  for (e in 1:N_events) {
    for (d in 1:D) {
      if (obs_mask[e, d] == 1) {
        z_aug[e, d] = (y_mat[e, d] == 1) ? exp(z_free[e, d]) : -exp(z_free[e, d]);
      } else {
        z_aug[e, d] = z_free[e, d];
      }
    }
  }

  matrix[N_events, D] mu_mat = X_event * beta;
  for (e in 1:N_events) {
    for (d in 1:D) {
      mu_mat[e, d] += hospital_effect[d, h_ev[e]]
                    + patient_effect[d, p_ev[e]]
                    + admission_effect[d, a_ev[e]];
    }
  }
}
model {
  to_vector(beta)          ~ normal(0, prior_beta_sd);
  tau_hospital             ~ normal(0, prior_tau_sd);
  tau_patient              ~ normal(0, prior_tau_sd);
  tau_admission            ~ normal(0, prior_tau_sd);
  to_vector(z_hospital)    ~ std_normal();
  to_vector(z_patient)     ~ std_normal();
  to_vector(z_admission)   ~ std_normal();
  L_corr_hospital          ~ lkj_corr_cholesky(lkj_eta);
  L_corr_patient           ~ lkj_corr_cholesky(lkj_eta);
  L_corr_admission         ~ lkj_corr_cholesky(lkj_eta);
  L_Omega                  ~ lkj_corr_cholesky(lkj_eta);

  for (e in 1:N_events) {
    target += multi_normal_cholesky_lpdf(z_aug[e]' | mu_mat[e]', L_Omega);
  }
  for (n in 1:N) {
    target += z_free[ev_idx[n], d_idx[n]];
  }
}
generated quantities {
  corr_matrix[D] Omega        = multiply_lower_tri_self_transpose(L_Omega);
  corr_matrix[D] R_hospital   = multiply_lower_tri_self_transpose(L_corr_hospital);
  corr_matrix[D] R_patient    = multiply_lower_tri_self_transpose(L_corr_patient);
  corr_matrix[D] R_admission  = multiply_lower_tri_self_transpose(L_corr_admission);
}
)"
}


# ---------------------------------------------------------------------------
# Stan model variants -- identity residual structure
# Classes are conditionally independent given fixed and random effects.
# L_Omega is NOT estimated; residual covariance = I_D.
# ---------------------------------------------------------------------------

.amr_probit_stan_1re_identity <- function() {
  r"(
// Multivariate probit -- one grouping level -- identity residual structure
// Classes are conditionally independent; L_Omega is not estimated.
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> K;
  matrix[N_events, K] X_event;
  array[N_events] int<lower=1,upper=H> h_ev;
  array[N_events, D] int<lower=0,upper=1> obs_mask;
  array[N_events, D] int<lower=0,upper=1> y_mat;
  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[D, H] z_hospital;
  vector<lower=0>[D] tau_hospital;
  cholesky_factor_corr[D] L_corr_hospital;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, H] hospital_effect =
    diag_pre_multiply(tau_hospital, L_corr_hospital) * z_hospital;

  matrix[N_events, D] z_aug;
  for (e in 1:N_events) {
    for (d in 1:D) {
      if (obs_mask[e, d] == 1) {
        z_aug[e, d] = (y_mat[e, d] == 1) ? exp(z_free[e, d]) : -exp(z_free[e, d]);
      } else {
        z_aug[e, d] = z_free[e, d];
      }
    }
  }

  matrix[N_events, D] mu_mat = X_event * beta;
  for (e in 1:N_events)
    for (d in 1:D)
      mu_mat[e, d] += hospital_effect[d, h_ev[e]];
}
model {
  to_vector(beta)         ~ normal(0, prior_beta_sd);
  tau_hospital            ~ normal(0, prior_tau_sd);
  to_vector(z_hospital)   ~ std_normal();
  L_corr_hospital         ~ lkj_corr_cholesky(lkj_eta);

  for (e in 1:N_events)
    for (d in 1:D)
      target += normal_lpdf(z_aug[e, d] | mu_mat[e, d], 1.0);

  for (n in 1:N)
    target += z_free[ev_idx[n], d_idx[n]];
}
generated quantities {
  corr_matrix[D] R_hospital = multiply_lower_tri_self_transpose(L_corr_hospital);
}
)"
}

.amr_probit_stan_2re_identity <- function() {
  r"(
// Multivariate probit -- two grouping levels -- identity residual structure
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> Pt;
  int<lower=1> K;
  matrix[N_events, K] X_event;
  array[N_events] int<lower=1,upper=H>  h_ev;
  array[N_events] int<lower=1,upper=Pt> p_ev;
  array[N_events, D] int<lower=0,upper=1> obs_mask;
  array[N_events, D] int<lower=0,upper=1> y_mat;
  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[D, H] z_hospital;
  vector<lower=0>[D] tau_hospital;
  cholesky_factor_corr[D] L_corr_hospital;
  matrix[D, Pt] z_patient;
  vector<lower=0>[D] tau_patient;
  cholesky_factor_corr[D] L_corr_patient;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, H]  hospital_effect =
    diag_pre_multiply(tau_hospital, L_corr_hospital) * z_hospital;
  matrix[D, Pt] patient_effect  =
    diag_pre_multiply(tau_patient,  L_corr_patient)  * z_patient;

  matrix[N_events, D] z_aug;
  for (e in 1:N_events) {
    for (d in 1:D) {
      if (obs_mask[e, d] == 1) {
        z_aug[e, d] = (y_mat[e, d] == 1) ? exp(z_free[e, d]) : -exp(z_free[e, d]);
      } else {
        z_aug[e, d] = z_free[e, d];
      }
    }
  }

  matrix[N_events, D] mu_mat = X_event * beta;
  for (e in 1:N_events)
    for (d in 1:D)
      mu_mat[e, d] += hospital_effect[d, h_ev[e]] + patient_effect[d, p_ev[e]];
}
model {
  to_vector(beta)        ~ normal(0, prior_beta_sd);
  tau_hospital           ~ normal(0, prior_tau_sd);
  tau_patient            ~ normal(0, prior_tau_sd);
  to_vector(z_hospital)  ~ std_normal();
  to_vector(z_patient)   ~ std_normal();
  L_corr_hospital        ~ lkj_corr_cholesky(lkj_eta);
  L_corr_patient         ~ lkj_corr_cholesky(lkj_eta);

  for (e in 1:N_events)
    for (d in 1:D)
      target += normal_lpdf(z_aug[e, d] | mu_mat[e, d], 1.0);

  for (n in 1:N)
    target += z_free[ev_idx[n], d_idx[n]];
}
generated quantities {
  corr_matrix[D] R_hospital = multiply_lower_tri_self_transpose(L_corr_hospital);
  corr_matrix[D] R_patient  = multiply_lower_tri_self_transpose(L_corr_patient);
}
)"
}

.amr_probit_stan_3re_identity <- function() {
  r"(
// Multivariate probit -- three grouping levels -- identity residual structure
data {
  int<lower=1> N;
  int<lower=1> N_events;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> Pt;
  int<lower=1> Adm;
  int<lower=1> K;
  matrix[N_events, K] X_event;
  array[N_events] int<lower=1,upper=H>   h_ev;
  array[N_events] int<lower=1,upper=Pt>  p_ev;
  array[N_events] int<lower=1,upper=Adm> a_ev;
  array[N_events, D] int<lower=0,upper=1> obs_mask;
  array[N_events, D] int<lower=0,upper=1> y_mat;
  array[N] int<lower=1,upper=N_events> ev_idx;
  array[N] int<lower=1,upper=D> d_idx;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[D, H]   z_hospital;
  vector<lower=0>[D] tau_hospital;
  cholesky_factor_corr[D] L_corr_hospital;
  matrix[D, Pt]  z_patient;
  vector<lower=0>[D] tau_patient;
  cholesky_factor_corr[D] L_corr_patient;
  matrix[D, Adm] z_admission;
  vector<lower=0>[D] tau_admission;
  cholesky_factor_corr[D] L_corr_admission;
  matrix[N_events, D] z_free;
}
transformed parameters {
  matrix[D, H]   hospital_effect  =
    diag_pre_multiply(tau_hospital,  L_corr_hospital)  * z_hospital;
  matrix[D, Pt]  patient_effect   =
    diag_pre_multiply(tau_patient,   L_corr_patient)   * z_patient;
  matrix[D, Adm] admission_effect =
    diag_pre_multiply(tau_admission, L_corr_admission) * z_admission;

  matrix[N_events, D] z_aug;
  for (e in 1:N_events) {
    for (d in 1:D) {
      if (obs_mask[e, d] == 1) {
        z_aug[e, d] = (y_mat[e, d] == 1) ? exp(z_free[e, d]) : -exp(z_free[e, d]);
      } else {
        z_aug[e, d] = z_free[e, d];
      }
    }
  }

  matrix[N_events, D] mu_mat = X_event * beta;
  for (e in 1:N_events)
    for (d in 1:D)
      mu_mat[e, d] += hospital_effect[d, h_ev[e]]
                    + patient_effect[d, p_ev[e]]
                    + admission_effect[d, a_ev[e]];
}
model {
  to_vector(beta)          ~ normal(0, prior_beta_sd);
  tau_hospital             ~ normal(0, prior_tau_sd);
  tau_patient              ~ normal(0, prior_tau_sd);
  tau_admission            ~ normal(0, prior_tau_sd);
  to_vector(z_hospital)    ~ std_normal();
  to_vector(z_patient)     ~ std_normal();
  to_vector(z_admission)   ~ std_normal();
  L_corr_hospital          ~ lkj_corr_cholesky(lkj_eta);
  L_corr_patient           ~ lkj_corr_cholesky(lkj_eta);
  L_corr_admission         ~ lkj_corr_cholesky(lkj_eta);

  for (e in 1:N_events)
    for (d in 1:D)
      target += normal_lpdf(z_aug[e, d] | mu_mat[e, d], 1.0);

  for (n in 1:N)
    target += z_free[ev_idx[n], d_idx[n]];
}
generated quantities {
  corr_matrix[D] R_hospital   = multiply_lower_tri_self_transpose(L_corr_hospital);
  corr_matrix[D] R_patient    = multiply_lower_tri_self_transpose(L_corr_patient);
  corr_matrix[D] R_admission  = multiply_lower_tri_self_transpose(L_corr_admission);
}
)"
}


# ---------------------------------------------------------------------------
# fit_bayesian_multivariate_probit()
# ---------------------------------------------------------------------------

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
#' \strong{Random effects} (\code{random_effects}): 1, 2, or 3 grouping column
#' names defining the clustering hierarchy (e.g. hospital; hospital + patient;
#' hospital + patient + admission). The first element is the upper-most level;
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
#' @param random_effects Character vector of length 1, 2, or 3. Grouping
#'   column names. First element: hospital (upper); second (optional): patient;
#'   third (optional): admission. Required.
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
#'   \code{tau_sd} (default 1.0), \code{lkj_eta} (default 2.0).
#' @param sampler_config Named list. Sampler settings forwarded to
#'   \code{cmdstanr::sample()}: \code{chains} (4), \code{iter_warmup} (1000),
#'   \code{iter_sampling} (1000), \code{adapt_delta} (NULL, uses Stan default),
#'   \code{max_treedepth} (NULL), \code{seed} (123), \code{parallel_chains}
#'   (NULL), \code{max_param_count} (NULL -- set to a positive integer to stop
#'   if approximate parameter count exceeds the threshold). Any additional entries
#'   are forwarded via \code{...}.
#' @param show_messages Logical. Print sampling progress. Default \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{cmdstanr::sample()}.
#'
#' @return Named list with elements: \code{draws}, \code{diagnostics},
#'   \code{diagnostics_detail}, \code{fit}, \code{data_long}, \code{index_maps},
#'   \code{X_design}, \code{X_event}, \code{event_re_idx}, \code{class_cols},
#'   \code{event_metadata}, \code{n_re_levels}, \code{upper_re_col},
#'   \code{middle_re_col}, \code{lower_re_col}, \code{patient_key_col},
#'   \code{admission_key_col}, \code{pathogen_col}, \code{pathogen_fitted},
#'   \code{residual_structure}, \code{estimand}, \code{prior_config_used},
#'   \code{sampler_config_used}, and \code{eligibility_report}.
#'
#'   \code{diagnostics} is a one-row monitored summary. The main diagnostic
#'   fields \code{max_rhat}, \code{min_ess_bulk}, and \code{min_ess_tail}
#'   are computed over the monitored parameters retained in \code{draws}.
#'   These exclude the latent \code{z_free} nuisance parameters used for
#'   probit data augmentation.
#'
#'   Full-fit diagnostic fields, including \code{max_rhat_all},
#'   \code{min_ess_bulk_all}, and \code{min_ess_tail_all}, are reported
#'   separately for visibility and include all Stan parameters, including
#'   \code{z_free}.
#'
#'   \code{diagnostics_detail} contains monitored-parameter, all-parameter,
#'   and chain-level diagnostic tables for detailed inspection.
#' @export
fit_bayesian_multivariate_probit <- function(
    event_class_data,
    class_cols,
    fixed_effects,
    random_effects,
    pathogen            = NULL,
    pathogen_col        = "pathogen",
    event_id_col        = "event_id",
    eligible_pairs      = NULL,
    outcome_col         = NULL,
    reserve_drug_cols   = NULL,
    panel_eligibility   = list(),
    residual_structure  = c("identity", "correlated"),
    estimand            = "observed_stewardship_event_mix",
    prior_config        = list(),
    sampler_config      = list(),
    show_messages       = TRUE,
    ...
) {
  residual_structure <- match.arg(residual_structure)

  # -- Package checks ---------------------------------------------------------
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop(paste0("Package 'cmdstanr' is required for Pathway 2.\n",
                "  install.packages('cmdstanr', ",
                "repos = c('https://mc-stan.org/r-packages/', getOption('repos')))"),
         call. = FALSE)
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed automatically with cmdstanr).",
         call. = FALSE)

  # -- Data checks ------------------------------------------------------------
  if (!is.data.frame(event_class_data) || nrow(event_class_data) == 0L)
    stop("`event_class_data` must be a non-empty data frame.", call. = FALSE)

  # -- Validate class_cols ----------------------------------------------------
  if (missing(class_cols) || length(class_cols) == 0L)
    stop("`class_cols` is required.", call. = FALSE)
  miss_cls <- setdiff(class_cols, names(event_class_data))
  if (length(miss_cls) > 0L)
    stop(sprintf("class_cols not found in event_class_data: %s",
                 paste(miss_cls, collapse = ", ")), call. = FALSE)

  # -- Validate fixed_effects -------------------------------------------------
  if (missing(fixed_effects) || is.null(fixed_effects) || length(fixed_effects) == 0L)
    stop("`fixed_effects` is required (no default).", call. = FALSE)
  miss_fe <- setdiff(fixed_effects, names(event_class_data))
  if (length(miss_fe) > 0L)
    stop(sprintf("fixed_effects column(s) not found: %s",
                 paste(miss_fe, collapse = ", ")), call. = FALSE)

  # -- Validate random_effects ------------------------------------------------
  if (missing(random_effects) || is.null(random_effects) || length(random_effects) == 0L)
    stop("`random_effects` is required (no default). Supply 1-3 grouping column names.",
         call. = FALSE)
  if (length(random_effects) > 3L)
    stop("`random_effects` accepts 1, 2, or 3 grouping columns (hospital; +patient; +admission).",
         call. = FALSE)
  miss_re <- setdiff(random_effects, names(event_class_data))
  if (length(miss_re) > 0L)
    stop(sprintf("random_effects column(s) not found: %s",
                 paste(miss_re, collapse = ", ")), call. = FALSE)

  n_re_levels  <- length(random_effects)
  upper_re_col <- random_effects[1L]
  middle_re_col <- if (n_re_levels >= 2L) random_effects[2L] else NULL
  lower_re_col  <- if (n_re_levels == 3L) random_effects[3L] else NULL

  # -- Validate pathogen_col and optionally filter to one pathogen -----------
  if (!pathogen_col %in% names(event_class_data))
    stop(sprintf("pathogen_col '%s' not found in event_class_data.", pathogen_col),
         call. = FALSE)

  event_data <- event_class_data
  pathogen_fitted <- NULL
  if (!is.null(pathogen)) {
    if (!pathogen %in% event_data[[pathogen_col]])
      stop(sprintf("pathogen '%s' not found in column '%s'.", pathogen, pathogen_col),
           call. = FALSE)
    event_data <- event_data[event_data[[pathogen_col]] == pathogen, , drop = FALSE]
    pathogen_fitted <- pathogen
    message(sprintf("[fit_bayesian_multivariate_probit] Filtered to pathogen '%s': %d events.",
                    pathogen, nrow(event_data)))
  }

  # -- Enforce single-pathogen fit (mandatory) --------------------------------
  n_pathogens <- dplyr::n_distinct(event_data[[pathogen_col]])
  if (n_pathogens != 1L)
    stop(sprintf(
      paste0("Bayesian multivariate-probit fitting requires exactly one pathogen ",
             "(found %d). Supply `pathogen = '<name>'` or pre-filter event_class_data."),
      n_pathogens), call. = FALSE)

  # -- Resolve prior_config ---------------------------------------------------
  pc <- list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0)
  for (nm in names(prior_config)) pc[[nm]] <- prior_config[[nm]]
  if (!is.numeric(pc$beta_sd) || pc$beta_sd <= 0)
    stop("`prior_config$beta_sd` must be a positive number.", call. = FALSE)
  if (!is.numeric(pc$tau_sd) || pc$tau_sd <= 0)
    stop("`prior_config$tau_sd` must be a positive number.", call. = FALSE)
  if (!is.numeric(pc$lkj_eta) || pc$lkj_eta < 1)
    stop("`prior_config$lkj_eta` must be >= 1.", call. = FALSE)

  # -- Resolve sampler_config -------------------------------------------------
  sc_defaults <- list(chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
                      seed = 123L, adapt_delta = NULL, max_treedepth = NULL,
                      parallel_chains = NULL, max_param_count = NULL)
  for (nm in names(sampler_config)) sc_defaults[[nm]] <- sampler_config[[nm]]
  sc <- sc_defaults

  message(sprintf(
    "[fit_bayesian_multivariate_probit] Priors: beta~N(0,%.2g) | tau~HN(0,%.2g) | LKJ(%.2g)",
    pc$beta_sd, pc$tau_sd, pc$lkj_eta))

  # -- Remove reserve drug columns --------------------------------------------
  if (!is.null(reserve_drug_cols)) {
    class_cols <- setdiff(class_cols, reserve_drug_cols)
    message(sprintf("[fit_bayesian_multivariate_probit] %d reserve drug column(s) excluded.",
                    length(reserve_drug_cols)))
  }
  if (length(class_cols) < 2L)
    stop("At least 2 class columns are required for the multivariate probit.", call. = FALSE)

  # -- Filter to eligible hospital-pathogen pairs -----------------------------
  if (!is.null(eligible_pairs)) {
    ep_req <- c(upper_re_col, pathogen_col)
    miss_ep <- setdiff(ep_req, names(eligible_pairs))
    if (length(miss_ep) > 0L)
      stop(sprintf("`eligible_pairs` must have columns: %s", paste(ep_req, collapse = ", ")),
           call. = FALSE)
    event_data <- dplyr::semi_join(event_data, eligible_pairs, by = ep_req)
    if (nrow(event_data) == 0L)
      stop("No events remain after filtering to eligible_pairs.", call. = FALSE)
    message(sprintf("[fit_bayesian_multivariate_probit] %d events after eligible_pairs filter.",
                    nrow(event_data)))
  }

  # -- Drop all-NA event rows (events with no observed AST in any class_col) --
  has_any_obs <- rowSums(!is.na(event_data[, class_cols, drop = FALSE])) > 0L
  n_all_na <- sum(!has_any_obs)
  if (n_all_na > 0L) {
    warning(sprintf(
      paste0("%d event(s) have no observed AST values across all selected class ",
             "columns and will be removed before model construction. Ensure this ",
             "is expected (e.g. events with only reserve-drug results)."),
      n_all_na), call. = FALSE)
    event_data <- event_data[has_any_obs, , drop = FALSE]
    if (nrow(event_data) == 0L)
      stop("No events remain after removing all-NA event rows.", call. = FALSE)
    message(sprintf("[fit_bayesian_multivariate_probit] %d events remain after all-NA filter.",
                    nrow(event_data)))
  }

  # -- Panel eligibility report (warn, do not filter) -------------------------
  pe <- list(min_tested = 30L, min_resistant = 5L, min_susceptible = 5L,
             min_pairwise_cotested = 30L)
  for (nm in names(panel_eligibility)) pe[[nm]] <- panel_eligibility[[nm]]

  eligibility_report <- .validate_panel_eligibility(
    event_data, class_cols, upper_re_col, pathogen_col,
    min_tested = pe$min_tested, min_resistant = pe$min_resistant,
    min_susceptible = pe$min_susceptible
  )

  # Pairwise co-testing overlap check
  elig_action <- .null_default(panel_eligibility$eligibility_action, "warn")
  co_test_report <- .check_pairwise_cotesting(
    event_data, class_cols, upper_re_col,
    min_pairwise_cotested = pe$min_pairwise_cotested
  )
  n_low_cotested <- sum(!co_test_report$sufficient, na.rm = TRUE)
  if (n_low_cotested > 0L) {
    if (identical(residual_structure, "correlated")) {
      msg <- sprintf(
        paste0("%d hospital x class-pair combination(s) have fewer than %d co-tested events. ",
               "Omega_{jk} for these pairs may be informed mainly by the LKJ prior. ",
               "See $eligibility_report$pairwise. Consider identity residuals or dropping classes with low overlap."),
        n_low_cotested, pe$min_pairwise_cotested)
    } else {
      msg <- sprintf(
        paste0("%d hospital x class-pair combination(s) have fewer than %d co-tested events. ",
               "Residual Omega is not estimated because residual_structure='identity'. ",
               "Joint profile estimates for weakly co-tested class pairs may be less data-supported. ",
               "See $eligibility_report$pairwise."),
        n_low_cotested, pe$min_pairwise_cotested)
    }
    if (elig_action == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  n_ineligible <- sum(!eligibility_report$eligible)
  if (n_ineligible > 0L)
    warning(sprintf(
      paste0("%d hospital x class cell(s) below marginal eligibility thresholds ",
             "(min_tested=%d, min_resistant=%d, min_susceptible=%d). ",
             "Results for these cells may be prior-dominated. See $eligibility_report."),
      n_ineligible, pe$min_tested, pe$min_resistant, pe$min_susceptible),
      call. = FALSE)

  # -- Build globally unique composite keys for nested RE --------------------
  # Prevents patient ID collisions across hospitals.
  if (n_re_levels >= 2L) {
    event_data$.patient_key <- paste(
      event_data[[upper_re_col]], event_data[[middle_re_col]], sep = "::"
    )
  }
  if (n_re_levels == 3L) {
    event_data$.admission_key <- paste(
      event_data[[upper_re_col]], event_data[[middle_re_col]],
      event_data[[lower_re_col]], sep = "::"
    )
  }

  # -- Validate event identifier (real uniqueness check, not seq_len trick) ---
  if (!event_id_col %in% names(event_data)) {
    # No external event ID: assign a warning-level fallback.
    warning(sprintf(
      paste0("event_id_col '%s' not found. Using row position as event key. ",
             "Ensure event_class_data has one row per organism-event -- no duplicates."),
      event_id_col), call. = FALSE)
    event_data[[event_id_col]] <- seq_len(nrow(event_data))
  } else {
    if (anyDuplicated(event_data[[event_id_col]]))
      stop(sprintf(
        paste0("event_id_col '%s' contains duplicate values. ",
               "Each row must correspond to exactly one organism-event."),
        event_id_col), call. = FALSE)
  }

  # -- Canonical event index (assigned BEFORE pivoting to long) ---------------
  event_data$.event_idx <- seq_len(nrow(event_data))
  N_events <- nrow(event_data)

  if (N_events > 1000L)
    warning(sprintf(
      paste0("N_events = %d with D = %d classes adds %d latent parameters (z_free). ",
             "Sampling may be slow. Use parallel_chains in sampler_config."),
      N_events, length(class_cols), N_events * length(class_cols)), call. = FALSE)

  # -- Columns to carry into long format ---------------------------------------
  patient_col_use   <- if (n_re_levels >= 2L) ".patient_key"   else NULL
  admission_col_use <- if (n_re_levels == 3L) ".admission_key" else NULL

  meta_carry <- unique(c(
    ".event_idx", upper_re_col, patient_col_use, admission_col_use, pathogen_col,
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

  if (nrow(data_long) == 0L)
    stop("No observed (event x class) pairs after removing NAs.", call. = FALSE)

  # -- Fail on fixed-effect missingness (do NOT impute; stop with clear msg) --
  fe_df_check <- data_long[, fixed_effects, drop = FALSE]
  na_cols <- vapply(fixed_effects, function(cc) sum(is.na(fe_df_check[[cc]])), integer(1L))
  if (any(na_cols > 0L)) {
    bad <- paste(sprintf("'%s' (%d NA)", names(na_cols)[na_cols > 0L], na_cols[na_cols > 0L]),
                 collapse = ", ")
    stop(sprintf(
      paste0("Fixed-effect column(s) contain NA values: %s. ",
             "Resolve missingness before calling fit_bayesian_multivariate_probit(). ",
             "Options: complete-case filter, median + indicator, explicit 'Unknown' level, ",
             "or multiple imputation -- but the decision belongs in the analysis config."),
      bad), call. = FALSE)
  }
  for (cc in fixed_effects)
    if (is.character(fe_df_check[[cc]])) fe_df_check[[cc]] <- factor(fe_df_check[[cc]])

  # -- Build integer index maps -----------------------------------------------
  class_levels <- class_cols
  upper_levels <- sort(unique(data_long[[upper_re_col]]))

  data_long <- data_long %>%
    dplyr::mutate(
      d_idx  = match(.data$antibiotic_class, class_levels),
      h_idx  = match(.data[[upper_re_col]],  upper_levels),
      ev_idx = as.integer(.data$.event_idx)
    )

  patient_levels   <- NULL
  admission_levels <- NULL
  Pt  <- NULL
  Adm <- NULL

  if (n_re_levels >= 2L) {
    patient_levels <- sort(unique(data_long[[patient_col_use]]))
    data_long <- data_long %>%
      dplyr::mutate(p_idx = match(.data[[patient_col_use]], patient_levels))
    Pt <- length(patient_levels)
    if (Pt * length(class_levels) > 10000L)
      warning(sprintf(
        "Large patient-RE block: %d patients x %d classes = %d parameters.",
        Pt, length(class_levels), Pt * length(class_levels)), call. = FALSE)
  }

  if (n_re_levels == 3L) {
    admission_levels <- sort(unique(data_long[[admission_col_use]]))
    data_long <- data_long %>%
      dplyr::mutate(a_idx = match(.data[[admission_col_use]], admission_levels))
    Adm <- length(admission_levels)
  }

  # -- Validate ev_idx (must span 1..N_events without gaps) -------------------
  stopifnot(all(data_long$ev_idx >= 1L & data_long$ev_idx <= N_events))
  stopifnot(!anyDuplicated(event_data$.event_idx))

  # -- Fixed-effects design matrix (event-level, not observation-level) -------
  X_long <- stats::model.matrix(~ ., data = fe_df_check)

  D <- length(class_levels)
  H <- length(upper_levels)
  K <- ncol(X_long)

  # Derive event-level arrays from the observation-level long table.
  # ev_idx in data_long is 1-based; one unique row per event.
  first_obs_per_event <- match(seq_len(N_events), data_long$ev_idx)
  stopifnot(!anyNA(first_obs_per_event))
  X_event_mat <- X_long[first_obs_per_event, , drop = FALSE]   # N_events x K
  h_ev_arr    <- as.integer(data_long$h_idx[first_obs_per_event])

  # Build wide AST arrays for multivariate-probit data augmentation
  obs_mask_mat <- matrix(0L, nrow = N_events, ncol = D)
  y_mat_mat    <- matrix(0L, nrow = N_events, ncol = D)
  for (i in seq_len(nrow(data_long))) {
    e <- data_long$ev_idx[i]
    d <- data_long$d_idx[i]
    obs_mask_mat[e, d] <- 1L
    y_mat_mat[e, d]    <- as.integer(data_long$resistance_binary[i])
  }

  message(sprintf(
    "[fit_bayesian_multivariate_probit] %d obs | %d events | D=%d | H=%d | RE levels=%d",
    nrow(data_long), N_events, D, H, n_re_levels))

  # -- Build Stan data list ---------------------------------------------------
  stan_data <- list(
    N             = nrow(data_long),         # observation pairs (for Jacobian)
    N_events      = N_events,
    D             = D,
    H             = H,
    K             = K,
    X_event       = unname(X_event_mat),     # N_events x K
    h_ev          = h_ev_arr,                # N_events
    obs_mask      = obs_mask_mat,            # N_events x D
    y_mat         = y_mat_mat,               # N_events x D
    ev_idx        = as.integer(data_long$ev_idx),
    d_idx         = as.integer(data_long$d_idx),
    prior_beta_sd = as.numeric(pc$beta_sd),
    prior_tau_sd  = as.numeric(pc$tau_sd),
    lkj_eta       = as.numeric(pc$lkj_eta)
  )
  if (n_re_levels >= 2L) {
    p_ev_arr     <- as.integer(data_long$p_idx[first_obs_per_event])
    stan_data$Pt  <- as.integer(Pt)
    stan_data$p_ev <- p_ev_arr
  }
  if (n_re_levels == 3L) {
    a_ev_arr     <- as.integer(data_long$a_idx[first_obs_per_event])
    stan_data$Adm <- as.integer(Adm)
    stan_data$a_ev <- a_ev_arr
  }

  # -- Parameter-count preflight -----------------------------------------------
  n_z_free    <- N_events * D
  n_hosp_raw  <- H * D
  n_pt_raw    <- if (n_re_levels >= 2L) Pt * D else 0L
  n_adm_raw   <- if (n_re_levels == 3L) Adm * D else 0L
  n_beta      <- K * D
  n_tau       <- D * n_re_levels
  n_L_corr    <- (D * (D - 1L) %/% 2L) * n_re_levels
  n_L_omega   <- if (residual_structure == "correlated") D * (D - 1L) %/% 2L else 0L
  n_param_approx <- n_z_free + n_hosp_raw + n_pt_raw + n_adm_raw +
                    n_beta + n_tau + n_L_corr + n_L_omega
  message(sprintf(paste0(
    "[fit_bayesian_multivariate_probit] Approximate parameter count:\n",
    "  N_events=%d x D=%d z_free              = %d\n",
    "  H=%d x D=%d hospital RE               = %d\n",
    "  %s\n",
    "  %s\n",
    "  beta(%dx%d) + tau(%d) + L_corr(%d) + L_Omega(%d)\n",
    "  Total approx: %d"),
    N_events, D, n_z_free,
    H, D, n_hosp_raw,
    if (n_re_levels >= 2L)
      sprintf("Pt=%d x D=%d level-2 RE              = %d", Pt, D, n_pt_raw)
    else "level-2 RE                            = 0 (not used)",
    if (n_re_levels == 3L)
      sprintf("Adm=%d x D=%d level-3 RE             = %d", Adm, D, n_adm_raw)
    else "level-3 RE                            = 0 (not used)",
    K, D, n_tau, n_L_corr, n_L_omega,
    n_param_approx
  ))
  if (!is.null(sc$max_param_count) && n_param_approx > as.integer(sc$max_param_count))
    stop(sprintf(
      paste0("Approximate parameter count (%d) exceeds sampler_config$max_param_count (%d). ",
             "Reduce N_events, D, or number of RE groups, or increase the threshold."),
      n_param_approx, as.integer(sc$max_param_count)), call. = FALSE)

  # -- Select and compile Stan model ------------------------------------------
  stan_code <- switch(paste0(n_re_levels, "_", residual_structure),
    "1_correlated" = .amr_probit_stan_1re(),
    "2_correlated" = .amr_probit_stan_2re(),
    "3_correlated" = .amr_probit_stan_3re(),
    "1_identity"   = .amr_probit_stan_1re_identity(),
    "2_identity"   = .amr_probit_stan_2re_identity(),
    "3_identity"   = .amr_probit_stan_3re_identity()
  )
  if (is.null(stan_code))
    stop(sprintf("Unsupported combination: n_re_levels=%d, residual_structure='%s'.",
                 n_re_levels, residual_structure), call. = FALSE)
  stan_file <- cmdstanr::write_stan_file(stan_code)
  message(sprintf("[fit_bayesian_multivariate_probit] Compiling %d-RE Stan model...",
                  n_re_levels))
  mod <- cmdstanr::cmdstan_model(stan_file, compile = TRUE)

  # -- Build sample() call args -----------------------------------------------
  sample_args <- list(
    data          = stan_data,
    seed          = as.integer(sc$seed),
    chains        = as.integer(sc$chains),
    iter_warmup   = as.integer(sc$iter_warmup),
    iter_sampling = as.integer(sc$iter_sampling),
    refresh       = if (show_messages) 200L else 0L
  )
  if (!is.null(sc$adapt_delta))   sample_args$adapt_delta   <- sc$adapt_delta
  if (!is.null(sc$max_treedepth)) sample_args$max_treedepth <- sc$max_treedepth
  if (!is.null(sc$parallel_chains)) sample_args$parallel_chains <- as.integer(sc$parallel_chains)
  extra_args <- list(...)
  sample_args <- c(sample_args, extra_args)

  # -- Sample -----------------------------------------------------------------
  message(sprintf(
    "[fit_bayesian_multivariate_probit] Sampling: %d chains x (%d warmup + %d sampling)...",
    sc$chains, sc$iter_warmup, sc$iter_sampling))
  fit <- do.call(mod$sample, sample_args)

  # -- Extract draws (exclude z_free; it is large and not needed downstream) --
  keep_vars <- c("beta", "hospital_effect",
                 if (n_re_levels >= 2L) "patient_effect",
                 if (n_re_levels == 3L) "admission_effect",
                 if (residual_structure == "correlated") c("L_Omega", "Omega"),
                 if (n_re_levels >= 1L) "R_hospital",
                 if (n_re_levels >= 2L) "R_patient",
                 if (n_re_levels == 3L) "R_admission",
                 "tau_hospital",
                 if (n_re_levels >= 2L) "tau_patient",
                 if (n_re_levels == 3L) "tau_admission",
                 "lp__")
  keep_vars <- unique(keep_vars[!is.null(keep_vars)])
  draws_arr <- tryCatch(
    fit$draws(variables = keep_vars, format = "draws_array"),
    error = function(e) fit$draws(format = "draws_array")
  )

  .probit_parameter_group <- function(parameter) {
    out <- rep("other", length(parameter))
    out[grepl("^beta(\\[|$)", parameter)] <- "beta"
    out[grepl("^hospital_effect(\\[|$)", parameter)] <- "hospital_effect"
    out[grepl("^patient_effect(\\[|$)", parameter)] <- "patient_effect"
    out[grepl("^admission_effect(\\[|$)", parameter)] <- "admission_effect"
    out[grepl("^R_hospital(\\[|$)", parameter)] <- "R_hospital"
    out[grepl("^R_patient(\\[|$)", parameter)] <- "R_patient"
    out[grepl("^R_admission(\\[|$)", parameter)] <- "R_admission"
    out[grepl("^(Omega|L_Omega)(\\[|$)", parameter)] <- "Omega"
    out[grepl("^tau_hospital(\\[|$)", parameter)] <- "tau_hospital"
    out[grepl("^tau_patient(\\[|$)", parameter)] <- "tau_patient"
    out[grepl("^tau_admission(\\[|$)", parameter)] <- "tau_admission"
    out[grepl("^z_free(\\[|$)", parameter)] <- "z_free"
    out[parameter == "lp__"] <- "lp__"
    out
  }
  .combine_diag_tables <- function(rhat_tbl, ess_tbl) {
    out <- merge(rhat_tbl, ess_tbl, by = "variable", all = TRUE)
    out$parameter_group <- .probit_parameter_group(out$variable)
    out[, c("variable", "parameter_group",
            setdiff(names(out), c("variable", "parameter_group"))),
        drop = FALSE]
  }
  .diag_extreme <- function(tbl, value_col, direction = c("max", "min")) {
    direction <- match.arg(direction)
    if (is.null(tbl) || !value_col %in% names(tbl))
      return(list(value = NA_real_, parameter = NA_character_, group = NA_character_))
    vals <- suppressWarnings(as.numeric(tbl[[value_col]]))
    ok <- !is.na(vals)
    if (!any(ok))
      return(list(value = NA_real_, parameter = NA_character_, group = NA_character_))
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
  ess_tbl_monitored  <- fit$summary(variables = keep_vars, "ess_bulk", "ess_tail")
  rhat_tbl_all       <- fit$summary(variables = NULL, "rhat")
  ess_tbl_all        <- fit$summary(variables = NULL, "ess_bulk", "ess_tail")
  monitored_diag_tbl <- .combine_diag_tables(rhat_tbl_monitored, ess_tbl_monitored)
  all_diag_tbl       <- .combine_diag_tables(rhat_tbl_all, ess_tbl_all)
  samp_diag <- fit$sampler_diagnostics(format = "matrix")

  n_divergent         <- sum(samp_diag[, "divergent__"], na.rm = TRUE)
  effective_max_td    <- as.integer(.null_default(sc$max_treedepth, 10L))
  n_treedepth         <- if ("treedepth__" %in% colnames(samp_diag))
                           sum(samp_diag[, "treedepth__"] >= effective_max_td,
                               na.rm = TRUE)
                         else NA_integer_
  max_rhat_monitored <- .diag_extreme(monitored_diag_tbl, "rhat", "max")
  min_bulk_monitored <- .diag_extreme(monitored_diag_tbl, "ess_bulk", "min")
  min_tail_monitored <- .diag_extreme(monitored_diag_tbl, "ess_tail", "min")
  max_rhat_all <- .diag_extreme(all_diag_tbl, "rhat", "max")
  min_bulk_all <- .diag_extreme(all_diag_tbl, "ess_bulk", "min")
  min_tail_all <- .diag_extreme(all_diag_tbl, "ess_tail", "min")

  max_rhat     <- max_rhat_monitored$value
  min_ess_bulk <- min_bulk_monitored$value
  min_ess_tail <- min_tail_monitored$value

  # E-BFMI per chain (robust to missing chain__ column)
  ebfmi_vals <- tryCatch({
    if (!"energy__" %in% colnames(samp_diag)) {
      NA_real_
    } else {
      e_vals   <- samp_diag[, "energy__"]
      n_chains <- as.integer(sc$chains)
      chain_id <- if ("chain__" %in% colnames(samp_diag)) {
        samp_diag[, "chain__"]
      } else {
        n_per_chain <- length(e_vals) %/% n_chains
        rep(seq_len(n_chains), each = n_per_chain)[seq_along(e_vals)]
      }
      vapply(unique(chain_id), function(ch) {
        e_ch <- e_vals[chain_id == ch]
        if (length(e_ch) < 2L) return(NA_real_)
        var_e <- stats::var(e_ch)
        if (is.na(var_e) || var_e < .Machine$double.eps) return(NA_real_)
        mean(diff(e_ch)^2) / var_e
      }, numeric(1L))
    }
  }, error = function(e) NA_real_)
  ebfmi_min <- min(ebfmi_vals, na.rm = TRUE)
  ebfmi_mean <- mean(ebfmi_vals, na.rm = TRUE)
  if (is.infinite(ebfmi_min)) ebfmi_min <- NA_real_
  if (is.nan(ebfmi_mean)) ebfmi_mean <- NA_real_
  ebfmi_by_chain <- if (all(is.na(ebfmi_vals))) {
    NA_character_
  } else {
    paste(sprintf("chain%s=%.3f", seq_along(ebfmi_vals), ebfmi_vals), collapse = "; ")
  }
  chain_diag_tbl <- tibble::tibble(
    chain = seq_along(ebfmi_vals),
    ebfmi = as.numeric(ebfmi_vals)
  )

  converged_structural <-
    isTRUE(max_rhat_structural < 1.01) &&
    isTRUE(min_ess_bulk_structural >= 100) &&
    isTRUE(is.na(min_ess_tail_structural) || min_ess_tail_structural >= 100) &&
    isTRUE(n_divergent == 0L) &&
    isTRUE(is.na(ebfmi) || ebfmi >= 0.3)
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
    isTRUE(is.na(ebfmi) || ebfmi >= 0.3)
  # Narrower, more diagnostic signal than "converged_full is FALSE": TRUE only
  # when the structural parameters have converged cleanly AND the full scope
  # (i.e. the z_free latent block specifically) has not -- isolating "the
  # augmented-data nuisance variables mix less well" from "the model itself
  # has not converged" (the latter already shows up as converged_structural
  # = FALSE and needs no separate flag).
  latent_diagnostic_warning <-
    isTRUE(converged_structural) &&
    (isTRUE(max_rhat_full > 1.01) || isTRUE(min_ess_bulk_full < 100))

  # -- Warnings are based on the STRUCTURAL scope; z_free stragglers are ------
  # -- reported separately below so they don't masquerade as model failure. --
  if (max_rhat_structural > 1.01)
    warning(sprintf("Convergence concern: structural max Rhat = %.3f (> 1.01). %s",
                    max_rhat_structural,
                    "Increase iter_warmup or adapt_delta."), call. = FALSE)
  if (min_ess_bulk_structural < 100L)
    warning(sprintf("Low structural bulk ESS: %.0f (< 100). Profile probabilities may be unreliable.",
                    min_ess_bulk_structural), call. = FALSE)
  if (!is.na(min_ess_tail_structural) && min_ess_tail_structural < 100L)
    warning(sprintf("Low structural tail ESS: %.0f (< 100).", min_ess_tail_structural), call. = FALSE)
  if (!is.na(n_treedepth) && n_treedepth > 0L)
    warning(sprintf("%d iteration(s) saturated max treedepth. Increase max_treedepth.",
                    n_treedepth), call. = FALSE)
  if (!is.na(ebfmi_min) && ebfmi_min < 0.3)
    warning(sprintf("Low minimum chain E-BFMI = %.3f (< 0.3). Possible poor geometry.",
                    ebfmi_min), call. = FALSE)
  if (n_divergent > 0L)
    warning(sprintf("%d divergent transition(s). Increase adapt_delta.", n_divergent),
            call. = FALSE)
  if (latent_diagnostic_warning)
    message(sprintf(
      paste0("[fit_bayesian_multivariate_probit] NOTE: latent z_free block (%d parameters) has ",
             "diagnostic stragglers (full-scope max Rhat = %.3f, min ESS bulk = %.0f) even though ",
             "structural parameters converged cleanly. This is informational only and does NOT ",
             "affect converged_structural -- see $diagnostics$latent_diagnostic_warning."),
      n_z_free, max_rhat_full, min_ess_bulk_full))

  diag_tbl <- tibble::tibble(
    n_chains          = as.integer(sc$chains),
    iter_warmup       = as.integer(sc$iter_warmup),
    iter_sampling     = as.integer(sc$iter_sampling),
    n_re_levels       = as.integer(n_re_levels),
    n_observed_pairs  = nrow(data_long),
    n_events          = N_events,
    n_classes         = D,
    n_upper_groups    = H,
    n_lower_groups    = if (n_re_levels >= 2L) as.integer(Pt) else NA_integer_,
    n_admission_groups = if (n_re_levels == 3L) as.integer(Adm) else NA_integer_,
    n_divergent       = as.integer(n_divergent),
    n_treedepth_sat   = as.integer(n_treedepth),
    ebfmi             = round(ebfmi_min, 4L),
    ebfmi_min         = round(ebfmi_min, 4L),
    ebfmi_mean        = round(ebfmi_mean, 4L),
    ebfmi_by_chain    = ebfmi_by_chain,
    max_rhat          = round(max_rhat, 4L),
    min_ess_bulk      = round(min_ess_bulk, 1L),
    min_ess_tail      = round(min_ess_tail, 1L),
    parameter_with_max_rhat = max_rhat_monitored$parameter,
    parameter_group_with_max_rhat = max_rhat_monitored$group,
    parameter_with_min_ess_bulk = min_bulk_monitored$parameter,
    parameter_group_with_min_ess_bulk = min_bulk_monitored$group,
    parameter_with_min_ess_tail = min_tail_monitored$parameter,
    parameter_group_with_min_ess_tail = min_tail_monitored$group,
    max_rhat_monitored     = round(max_rhat_monitored$value, 4L),
    min_ess_bulk_monitored = round(min_bulk_monitored$value, 1L),
    min_ess_tail_monitored = round(min_tail_monitored$value, 1L),
    max_rhat_all           = round(max_rhat_all$value, 4L),
    min_ess_bulk_all       = round(min_bulk_all$value, 1L),
    min_ess_tail_all       = round(min_tail_all$value, 1L)
  )

  message(sprintf(
    paste0("[fit_bayesian_multivariate_probit] Done. structural max_Rhat=%.3f (full incl. z_free=%.3f) | ",
           "structural min_ESS_bulk=%.0f (full=%.0f) | divergent=%d | converged_structural=%s"),
    max_rhat_structural, max_rhat_full,
    min_ess_bulk_structural, min_ess_bulk_full,
    n_divergent, converged_structural))

  # Store event-level arrays needed by compute_event_profile_probabilities()
  event_re_idx <- list(h_ev = h_ev_arr)
  if (n_re_levels >= 2L) event_re_idx$p_ev <- p_ev_arr
  if (n_re_levels == 3L) event_re_idx$a_ev <- a_ev_arr

  list(
    draws              = draws_arr,
    diagnostics        = diag_tbl,
    diagnostics_detail = list(
      monitored_parameters = tibble::as_tibble(monitored_diag_tbl),
      all_parameters       = tibble::as_tibble(all_diag_tbl),
      chains               = chain_diag_tbl
    ),
    fit                = fit,
    data_long          = tibble::as_tibble(data_long),
    index_maps         = list(
      class_levels      = class_levels,
      upper_levels      = upper_levels,
      patient_levels    = patient_levels,
      admission_levels  = admission_levels
    ),
    X_design           = X_long,   # kept for back-compat; same as X_event aligned to obs
    X_event            = X_event_mat,   # N_events x K (use this in simulation)
    event_re_idx       = event_re_idx,  # h_ev, p_ev, a_ev per event
    class_cols         = class_cols,
    event_metadata     = tibble::as_tibble(event_data),
    n_re_levels        = n_re_levels,
    upper_re_col       = upper_re_col,
    middle_re_col      = middle_re_col,
    lower_re_col       = lower_re_col,
    patient_key_col    = patient_col_use,
    admission_key_col  = admission_col_use,
    pathogen_col        = pathogen_col,
    pathogen_fitted     = pathogen_fitted,
    residual_structure  = residual_structure,
    estimand            = estimand,
    prior_config_used   = pc,
    sampler_config_used = sc,
    eligibility_report  = list(
      marginal  = eligibility_report,
      pairwise  = co_test_report
    )
  )
}


# ---------------------------------------------------------------------------
# compute_event_profile_probabilities()
# ---------------------------------------------------------------------------

#' Compute Posterior Resistance Profile Probabilities via MVN Simulation
#'
#' For each posterior draw, constructs the event-level linear predictor
#' \eqn{\mu_e} using fixed effects and correlated random effects from the
#' fitted model, simulates latent
#' \eqn{Z_e = \mu_e + L_\Omega\,\varepsilon_e}, \eqn{\varepsilon_e \sim N(0,I_D)},
#' converts to binary resistance profiles, and accumulates profile
#' probabilities. All \eqn{2^D} profiles appear in the output; profiles not
#' observed in simulation receive probability 0.
#'
#' \strong{Estimand:} The posterior predictive distribution over the observed
#' event case-mix -- the profile distribution you would see if you drew a new
#' event uniformly from the set of observed events (same covariate distribution
#' as the data). This is labelled \code{"observed_stewardship_event_mix"}.
#'
#' @param fitted_model List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param n_posterior_draws_for_profiles Integer. Number of posterior draws to
#'   use for profile simulation. Subsampled without replacement when total draws
#'   exceed this value. Each draw generates one simulated latent profile per
#'   event; rare profiles may be underestimated when this is small. For finer
#'   Monte Carlo resolution of rare profiles, consider increasing this value or
#'   adding \code{n_predictive_replicates_per_draw} in a future extension.
#'   Default \code{2000L}.
#' @param outcome_col Character or \code{NULL}. Patient outcome column in
#'   \code{fitted_model$event_metadata}. When \code{NULL}, all events are
#'   treated as having a known outcome and R_NF is not computed separately.
#' @param nonfatal_values Character vector. Outcome values for the non-fatal
#'   cohort (R_NF). Default covers common discharge/survival labels.
#' @param seed Integer. Random seed for MVN simulation. Default \code{123L}.
#'
#' @return Named list: \code{event_profiles} (event-level posterior means) and
#'   \code{aggregate_draws} (per-draw R_ALL and R_NF per hospital x pathogen x
#'   profile, used by \code{aggregate_profiles_for_daly()} for credible
#'   intervals). Both tibbles contain all \eqn{2^D} profiles per
#'   hospital-pathogen pair.
#' @export
compute_event_profile_probabilities <- function(
    fitted_model,
    n_posterior_draws_for_profiles = 2000L,
    outcome_col                    = NULL,
    nonfatal_values                = c("Discharged", "Survived", "Alive",
                                       "discharged", "survived", "alive"),
    seed                           = 123L
) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)

  set.seed(as.integer(seed))

  draws               <- fitted_model$draws
  class_cols          <- fitted_model$class_cols
  idx_maps            <- fitted_model$index_maps
  X_long              <- fitted_model$X_design
  data_long           <- fitted_model$data_long
  event_meta          <- fitted_model$event_metadata
  n_re_levels         <- fitted_model$n_re_levels
  upper_re_col        <- fitted_model$upper_re_col
  pathogen_col        <- fitted_model$pathogen_col
  residual_structure  <- .null_default(fitted_model$residual_structure, "identity")

  D <- length(idx_maps$class_levels)
  H <- length(idx_maps$upper_levels)

  # Prefer event-level arrays stored by fit function; fall back to deriving them
  X_event_sim  <- if (!is.null(fitted_model$X_event))
                    fitted_model$X_event
                  else fitted_model$X_design  # legacy fallback
  event_re_idx <- fitted_model$event_re_idx   # list: h_ev, [p_ev], [a_ev]

  K <- ncol(X_event_sim)

  # -- Thin draws to n_posterior_draws_for_profiles ---------------------------
  draws_mat <- posterior::as_draws_matrix(draws)
  n_total   <- nrow(draws_mat)
  S         <- min(as.integer(n_posterior_draws_for_profiles), n_total)
  draw_idx  <- if (S < n_total) sort(sample.int(n_total, S)) else seq_len(n_total)
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  # -- Helper: extract [S, d1, d2] array from draws ---------------------------
  .arr <- function(prefix, d1, d2) {
    cols <- as.vector(outer(seq_len(d1), seq_len(d2),
                            function(a, b) sprintf("%s[%d,%d]", prefix, a, b)))
    array(draws_mat[, cols, drop = FALSE], dim = c(S, d1, d2))
  }

  beta_arr     <- .arr("beta", K, D)                 # S x K x D
  hosp_eff_arr <- .arr("hospital_effect", D, H)      # S x D x H
  L_omega_arr  <- if (residual_structure == "correlated")
                    .arr("L_Omega", D, D)             # S x D x D (lower triangular)
                  else NULL

  has_patient   <- n_re_levels >= 2L
  has_admission <- n_re_levels == 3L
  Pt  <- if (has_patient)   length(idx_maps$patient_levels)   else 0L
  Adm <- if (has_admission) length(idx_maps$admission_levels) else 0L

  if (has_patient)   pt_eff_arr  <- .arr("patient_effect",   D, Pt)
  if (has_admission) adm_eff_arr <- .arr("admission_effect", D, Adm)

  # -- Canonical event table --------------------------------------------------
  stopifnot(".event_idx" %in% names(event_meta))
  stopifnot(!anyDuplicated(event_meta$.event_idx))

  # Use event_meta to identify which events appear in data_long (have observations)
  has_obs       <- event_meta$.event_idx %in% data_long$ev_idx
  event_meta_obs <- event_meta[has_obs, , drop = FALSE]
  N_ev           <- nrow(event_meta_obs)

  # Event-level design matrix and RE indices (from stored event_re_idx)
  ev_row_idx <- event_meta_obs$.event_idx           # 1-based canonical event indices
  X_event    <- X_event_sim[ev_row_idx, , drop = FALSE]    # N_ev x K
  h_events   <- as.integer(event_re_idx$h_ev[ev_row_idx])
  p_events   <- if (has_patient)   as.integer(event_re_idx$p_ev[ev_row_idx])   else NULL
  a_events   <- if (has_admission) as.integer(event_re_idx$a_ev[ev_row_idx])   else NULL

  # -- Outcome flags per event ------------------------------------------------
  oc_col_use <- if (!is.null(outcome_col) && outcome_col %in% names(event_meta_obs))
                  outcome_col
                else NULL

  if (!is.null(oc_col_use)) {
    ov               <- event_meta_obs[[oc_col_use]]
    is_known_outcome <- !is.na(ov)
    is_nonfatal      <- !is.na(ov) & (ov %in% nonfatal_values)
  } else {
    is_known_outcome <- rep(TRUE, N_ev)
    is_nonfatal      <- rep(TRUE, N_ev)
  }

  # -- Hospital-pathogen eligible class sets ----------------------------------
  hp_pairs <- unique(event_meta_obs[, c(upper_re_col, pathogen_col), drop = FALSE])
  hp_keys  <- paste(hp_pairs[[upper_re_col]], hp_pairs[[pathogen_col]], sep = "||")

  hp_class_d <- stats::setNames(lapply(seq_len(nrow(hp_pairs)), function(r) {
    h_nm <- hp_pairs[[upper_re_col]][r]
    k_nm <- hp_pairs[[pathogen_col]][r]
    sub  <- event_meta_obs[
      event_meta_obs[[upper_re_col]] == h_nm & event_meta_obs[[pathogen_col]] == k_nm,
      class_cols, drop = FALSE
    ]
    which(vapply(class_cols, function(cc) {
      v <- sub[[cc]]
      any(!is.na(v))
    }, logical(1L)))
  }), hp_keys)

  # Map hp_keys -> row indices in event_meta_obs (canonical, not event_meta!)
  hp_ev_idx <- stats::setNames(lapply(hp_keys, function(key) {
    parts <- strsplit(key, "||", fixed = TRUE)[[1L]]
    which(event_meta_obs[[upper_re_col]] == parts[1L] &
          event_meta_obs[[pathogen_col]]  == parts[2L])
  }), hp_keys)

  # Storage: per-hp, per-draw profile labels matrix [N_hp x S]
  hp_store <- stats::setNames(lapply(hp_keys, function(key) {
    n_hp <- length(hp_ev_idx[[key]])
    list(labels    = matrix(NA_character_, n_hp, S),
         known_out = is_known_outcome[hp_ev_idx[[key]]],
         nonfatal  = is_nonfatal[hp_ev_idx[[key]]])
  }), hp_keys)

  message(sprintf(
    "[compute_event_profile_probabilities] %d draws | %d events | %d hp-pairs | %d RE level(s)",
    S, N_ev, length(hp_keys), n_re_levels))

  # -- Main simulation loop ---------------------------------------------------
  for (s in seq_len(S)) {
    beta_s     <- matrix(beta_arr[s, , ],     nrow = K, ncol = D)
    hosp_eff_s <- matrix(hosp_eff_arr[s, , ], nrow = D, ncol = H)

    # mu_all: N_ev x D
    mu_all <- X_event %*% beta_s          # N_ev x D (fixed effects)
    for (ev_i in seq_len(N_ev)) {
      mu_all[ev_i, ] <- mu_all[ev_i, ] + hosp_eff_s[, h_events[ev_i]]
    }
    if (has_patient) {
      pt_eff_s <- matrix(pt_eff_arr[s, , ], nrow = D, ncol = Pt)
      for (ev_i in seq_len(N_ev))
        mu_all[ev_i, ] <- mu_all[ev_i, ] + pt_eff_s[, p_events[ev_i]]
    }
    if (has_admission) {
      adm_eff_s <- matrix(adm_eff_arr[s, , ], nrow = D, ncol = Adm)
      for (ev_i in seq_len(N_ev))
        mu_all[ev_i, ] <- mu_all[ev_i, ] + adm_eff_s[, a_events[ev_i]]
    }

    # Reconstruct full Omega for this draw (only needed for correlated structure)
    Omega_s_full <- if (residual_structure == "correlated") {
      L_omega_s <- matrix(L_omega_arr[s, , ], nrow = D, ncol = D)
      tcrossprod(L_omega_s)
    } else NULL

    for (key in hp_keys) {
      d_hp   <- hp_class_d[[key]]
      ev_idx <- hp_ev_idx[[key]]
      if (length(ev_idx) == 0L || length(d_hp) < 1L) next

      mu_hp <- mu_all[ev_idx, d_hp, drop = FALSE]   # N_hp x |d_hp|

      if (residual_structure == "correlated") {
        # Extract the correct correlation sub-matrix and its Cholesky factor.
        # A submatrix of L_Omega is NOT the Cholesky factor of Omega[d_hp, d_hp].
        Omega_hp <- Omega_s_full[d_hp, d_hp, drop = FALSE]
        Omega_hp <- (Omega_hp + t(Omega_hp)) / 2 + diag(1e-9, length(d_hp))
        L_hp <- tryCatch(t(chol(Omega_hp)), error = function(e) NULL)
        if (is.null(L_hp)) next
        Z_std <- matrix(stats::rnorm(nrow(mu_hp) * length(d_hp)), nrow = nrow(mu_hp))
        Z_hp  <- mu_hp + Z_std %*% t(L_hp)
      } else {
        # Identity residual structure: classes are conditionally independent
        Z_hp <- mu_hp + matrix(stats::rnorm(nrow(mu_hp) * length(d_hp)), nrow = nrow(mu_hp))
      }

      hp_store[[key]]$labels[, s] <- apply(
        (Z_hp > 0) * 1L, 1L,
        function(row) paste(ifelse(row == 1L, "R", "S"), collapse = "")
      )
    }
  }

  # -- Enumerate all 2^|d_hp| profiles per hp pair ----------------------------
  # All profiles appear in the output. Profiles not sampled in finite Monte Carlo
  # receive an estimated probability of 0, which may reflect simulation error
  # for rare profiles rather than true structural zero probability.
  event_profile_rows  <- list()
  aggregate_draw_rows <- list()

  for (key in hp_keys) {
    parts  <- strsplit(key, "||", fixed = TRUE)[[1L]]
    h_nm   <- parts[1L]; k_nm <- parts[2L]
    d_hp   <- hp_class_d[[key]]
    ev_idx <- hp_ev_idx[[key]]
    store  <- hp_store[[key]]
    n_hp   <- nrow(store$labels)
    if (n_hp == 0L) next

    valid_s <- which(colSums(!is.na(store$labels)) == n_hp)
    if (length(valid_s) == 0L) next
    lbl_v <- store$labels[, valid_s, drop = FALSE]

    class_set <- paste(class_cols[d_hp], collapse = "|")

    # Complete set of 2^K profiles for this class panel
    all_lbls <- if (length(d_hp) >= 1L) {
      enum_df <- enumerate_binary_profiles(class_cols[d_hp])
      enum_df$profile_delta
    } else {
      sort(unique(as.vector(lbl_v)))
    }

    # Event-level posterior mean profile probabilities
    for (ev_i in seq_len(n_hp)) {
      ev_lbl <- lbl_v[ev_i, ]
      probs  <- vapply(all_lbls, function(lbl) mean(ev_lbl == lbl, na.rm = TRUE),
                       numeric(1L))
      event_id_val <- event_meta_obs$.event_idx[ev_idx[ev_i]]
      for (li in seq_along(all_lbls)) {
        event_profile_rows[[length(event_profile_rows) + 1L]] <- tibble::tibble(
          !!upper_re_col    := h_nm,
          !!pathogen_col    := k_nm,
          event_idx          = event_id_val,
          profile_class_set  = class_set,
          profile_delta      = all_lbls[li],
          profile_probability = probs[li]
        )
      }
    }

    # Draw-level aggregate R_ALL and R_NF (complete profile enumeration)
    for (si in seq_along(valid_s)) {
      s_lbl <- lbl_v[, si]
      for (lbl in all_lbls) {
        aggregate_draw_rows[[length(aggregate_draw_rows) + 1L]] <- tibble::tibble(
          !!upper_re_col   := h_nm,
          !!pathogen_col   := k_nm,
          profile_class_set = class_set,
          profile_delta     = lbl,
          draw_s            = valid_s[si],
          R_ALL_s = if (any(store$known_out))
                      mean(s_lbl[store$known_out] == lbl, na.rm = TRUE)
                    else NA_real_,
          R_NF_s  = if (any(store$nonfatal))
                      mean(s_lbl[store$nonfatal]  == lbl, na.rm = TRUE)
                    else NA_real_
        )
      }
    }
  }

  message(sprintf(
    "[compute_event_profile_probabilities] Done. %d event-profile rows | %d draw-aggregate rows.",
    length(event_profile_rows), length(aggregate_draw_rows)))

  list(
    event_profiles  = dplyr::bind_rows(event_profile_rows),
    aggregate_draws = dplyr::bind_rows(aggregate_draw_rows)
  )
}

#' Aggregate Posterior Profile Draws into R_ALL and R_NF Summaries
#'
#' Summarises draw-level R_ALL and R_NF values from
#' \code{compute_event_profile_probabilities()} into posterior mean and
#' credible interval per hospital x pathogen x profile combination.
#' All \eqn{2^D} profiles are present in the output (inherited from the
#' complete enumeration in the simulation step).
#'
#' @param profile_output List returned by \code{compute_event_profile_probabilities()}.
#' @param hospital_col Character. Hospital column in \code{aggregate_draws}.
#'   Must match the upper RE column used during fitting. Default \code{"hospital"}.
#' @param pathogen_col Character. Pathogen column. Default \code{"pathogen"}.
#' @param estimand Character. Estimand label to attach to output.
#'   Default \code{"observed_stewardship_event_mix"}.
#' @param ci_level Numeric. Credible interval coverage. Default \code{0.95}.
#'
#' @return Tibble with one row per hospital x pathogen x profile:
#'   R_ALL and R_NF posterior mean and credible interval.
#' @export
aggregate_profiles_for_daly <- function(
    profile_output,
    hospital_col = "hospital",
    pathogen_col = "pathogen",
    estimand     = "observed_stewardship_event_mix",
    ci_level     = 0.95
) {
  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q

  if (!is.list(profile_output) ||
      !all(c("aggregate_draws", "event_profiles") %in% names(profile_output)))
    stop("`profile_output` must be the list returned by compute_event_profile_probabilities().",
         call. = FALSE)

  agg_draws <- profile_output$aggregate_draws
  if (nrow(agg_draws) == 0L) {
    warning("aggregate_draws is empty. Returning empty tibble.", call. = FALSE)
    return(tibble::tibble())
  }

  # Detect the actual hospital column name from the draws tibble
  if (!hospital_col %in% names(agg_draws)) {
    candidate <- setdiff(names(agg_draws),
                         c(pathogen_col, "profile_class_set", "profile_delta",
                           "draw_s", "R_ALL_s", "R_NF_s"))
    if (length(candidate) == 1L) {
      hospital_col <- candidate[1L]
      message(sprintf("[aggregate_profiles_for_daly] Using '%s' as hospital column.",
                      hospital_col))
    } else {
      stop(sprintf("hospital_col '%s' not found in aggregate_draws.", hospital_col),
           call. = FALSE)
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
      R_ALL_mean   = mean(.data$R_ALL_s,  na.rm = TRUE),
      R_ALL_lower  = stats::quantile(.data$R_ALL_s, lo_q, na.rm = TRUE),
      R_ALL_upper  = stats::quantile(.data$R_ALL_s, hi_q, na.rm = TRUE),
      R_NF_mean    = mean(.data$R_NF_s,   na.rm = TRUE),
      R_NF_lower   = stats::quantile(.data$R_NF_s,  lo_q, na.rm = TRUE),
      R_NF_upper   = stats::quantile(.data$R_NF_s,  hi_q, na.rm = TRUE),
      n_draws_used = dplyr::n(),
      .groups      = "drop"
    ) %>%
    dplyr::mutate(
      profile_label        = .data$profile_delta,
      profile_set_type     = "facility_bayesian_probit",
      estimand             = estimand,
      used_for_YLL         = TRUE,
      used_for_YLD         = TRUE,
      estimator            = "bayesian_multivariate_probit",
      low_draws_flag       = .data$n_draws_used < 100L,
      dplyr::across(dplyr::starts_with("R_"), ~ round(.x, 6L))
    ) %>%
    dplyr::select(
      dplyr::all_of(c(hospital_col, pathogen_col)),
      profile_set_type, estimand, profile_class_set, profile_delta, profile_label,
      R_ALL_mean, R_ALL_lower, R_ALL_upper,
      R_NF_mean,  R_NF_lower,  R_NF_upper,
      n_draws_used, used_for_YLL, used_for_YLD, estimator, low_draws_flag
    )

  message(sprintf(
    "[aggregate_profiles_for_daly] %d hospital-pathogen-profile rows.",
    nrow(result)))

  result
}


# ---------------------------------------------------------------------------
# estimate_resistance_profiles()  -- top-level dispatcher for both pathways
# ---------------------------------------------------------------------------

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
#' @param random_effects Character vector (1-3 elements). Pathway 2 only.
#'   Required. Elements: hospital; [+patient]; [+admission].
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
    method                         = c("convex", "bayesian"),
    # Pathway 1
    pairwise                       = NULL,
    panel_map                      = NULL,
    # Pathway 2 -- required with no defaults
    class_cols                     = NULL,
    fixed_effects                  = NULL,
    random_effects                 = NULL,
    pathogen                       = NULL,
    pathogen_col                   = "pathogen",
    eligible_pairs                 = NULL,
    outcome_col                    = NULL,
    nonfatal_values                = c("Discharged", "Survived", "Alive",
                                       "discharged", "survived", "alive"),
    panel_eligibility              = list(),
    residual_structure             = c("identity", "correlated"),
    estimand                       = "observed_stewardship_event_mix",
    prior_config                   = list(),
    sampler_config                 = list(),
    n_posterior_draws_for_profiles = 2000L
) {
  method             <- match.arg(method)
  residual_structure <- match.arg(residual_structure)

  config_used <- list(
    method                         = method,
    pathogen_col                   = pathogen_col,
    pathogen                       = pathogen,
    residual_structure             = residual_structure,
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
        pct_converged       = mean(.data$convergence_flag,    na.rm = TRUE),
        max_abs_residual    = max(.data$max_abs_residual,     na.rm = TRUE),
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
  if (!is.null(pathogen))
    message(sprintf("  Pathogen: %s", pathogen))

  if (is.null(class_cols))
    stop("`class_cols` is required for method = 'bayesian'.", call. = FALSE)
  if (is.null(fixed_effects))
    stop("`fixed_effects` is required for method = 'bayesian'.", call. = FALSE)
  if (is.null(random_effects))
    stop("`random_effects` is required for method = 'bayesian'.", call. = FALSE)

  fitted_mod <- fit_bayesian_multivariate_probit(
    event_class_data   = data,
    class_cols         = class_cols,
    fixed_effects      = fixed_effects,
    random_effects     = random_effects,
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
    hospital_col   = fitted_mod$upper_re_col,
    pathogen_col   = pathogen_col,
    estimand       = estimand
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
  } else NULL

  message("[estimate_resistance_profiles] Pathway 2 complete.")

  list(
    profiles      = profiles_tbl,
    eligibility   = eligibility_tbl,
    diagnostics   = fitted_mod$diagnostics,
    fitted_models = list(bayesian_probit = fitted_mod),
    config_used   = config_used
  )
}
