# daly_resistance_profiles.R
#
# Estimates antimicrobial resistance profiles for two complementary data sources.
# All profile probability distributions produced by either pathway feed directly
# into the DALY burden pipeline (YLL and YLD calculations).
#
# Pathway 1 — Convex optimisation (aggregate surveillance or GBD-style data)
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
# Pathway 2 — Bayesian hierarchical modelling (facility-level AST data)
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


# ===========================================================================
# Pathway 1 — Isolate-level engine
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
# Pathway 2 — Bayesian hierarchical multivariate probit
#
# Accepts wide-format event-level AST data (one row per organism-event, one
# column per antibiotic class) and estimates hospital-specific resistance
# profile probability distributions using a hierarchical Bayesian multivariate
# probit model. The model is fitted once across all hospitals simultaneously,
# enabling partial pooling of information through shared correlation structure
# and hyperpriors while still producing hospital-specific profile estimates.
#
# The model captures:
#   - Class-specific baseline resistance (fixed effects per class)
#   - Patient covariates: ICU/ward, HAI/CAI, specimen type, year (shared across classes)
#   - Hospital heterogeneity: hospital x class random effects (partial pooling)
#   - Patient-admission clustering: patient-admission random effects
#   - Correlated resistance across classes: residual LKJ correlation matrix Omega
#
# Hospitals with fewer eligible classes contribute narrower profiles; the
# function uses each hospital's observed class set (classes with at least one
# non-NA value in the data) to define hospital-specific profile dimensions.
# Profile labels are always anchored to the hospital-specific ordered class set.
#
# Outputs feed directly into the DALY burden pipeline:
#   R_ALL  all eligible infected events with known outcome  -> YLL
#   R_NF   non-fatal / discharged events only              -> YLD
# ===========================================================================




# ---------------------------------------------------------------------------
# Internal: Stan model code — two variants based on number of RE levels
# ---------------------------------------------------------------------------
# Priors are passed as data so the model compiles once and prior values
# can be changed freely between runs without recompilation.

.amr_probit_stan_1re <- function() {
  r"(
// Multivariate probit — one grouping level (upper clustering only)
data {
  int<lower=1> N;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> K;
  array[N] int<lower=0,upper=1> y;
  matrix[N, K] X;
  array[N] int<lower=1,upper=D> d_idx;
  array[N] int<lower=1,upper=H> h_idx;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[H, D] u_raw;
  vector<lower=0>[D] tau_h;
  cholesky_factor_corr[D] L_Omega;
}
transformed parameters {
  matrix[H, D] u;
  for (d in 1:D)
    u[, d] = tau_h[d] * u_raw[, d];
}
model {
  to_vector(beta)  ~ normal(0, prior_beta_sd);
  tau_h            ~ normal(0, prior_tau_sd);
  to_vector(u_raw) ~ std_normal();
  L_Omega          ~ lkj_corr_cholesky(lkj_eta);
  {
    vector[N] mu;
    for (n in 1:N)
      mu[n] = dot_product(X[n], beta[, d_idx[n]])
              + u[h_idx[n], d_idx[n]];
    y ~ bernoulli(Phi(mu));
  }
}
generated quantities {
  corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
)"
}

.amr_probit_stan_2re <- function() {
  r"(
// Multivariate probit — two grouping levels (upper + lower clustering)
data {
  int<lower=1> N;
  int<lower=2> D;
  int<lower=1> H;
  int<lower=1> Pt;
  int<lower=1> K;
  array[N] int<lower=0,upper=1> y;
  matrix[N, K] X;
  array[N] int<lower=1,upper=D> d_idx;
  array[N] int<lower=1,upper=H> h_idx;
  array[N] int<lower=1,upper=Pt> p_idx;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
  real<lower=1> lkj_eta;
}
parameters {
  matrix[K, D] beta;
  matrix[H, D] u_raw;
  matrix[Pt, D] v_raw;
  vector<lower=0>[D] tau_h;
  vector<lower=0>[D] tau_pt;
  cholesky_factor_corr[D] L_Omega;
}
transformed parameters {
  matrix[H, D] u;
  matrix[Pt, D] v;
  for (d in 1:D) {
    u[, d] = tau_h[d]  * u_raw[, d];
    v[, d] = tau_pt[d] * v_raw[, d];
  }
}
model {
  to_vector(beta)  ~ normal(0, prior_beta_sd);
  tau_h            ~ normal(0, prior_tau_sd);
  tau_pt           ~ normal(0, prior_tau_sd);
  to_vector(u_raw) ~ std_normal();
  to_vector(v_raw) ~ std_normal();
  L_Omega          ~ lkj_corr_cholesky(lkj_eta);
  {
    vector[N] mu;
    for (n in 1:N)
      mu[n] = dot_product(X[n], beta[, d_idx[n]])
              + u[h_idx[n], d_idx[n]]
              + v[p_idx[n], d_idx[n]];
    y ~ bernoulli(Phi(mu));
  }
}
generated quantities {
  corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L_Omega);
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
#' \strong{Model:}
#' \deqn{Y_{id} = \mathbf{1}(Z_{id} > 0), \quad
#'   Z_i \sim \text{MVN}(\mu_i,\,\Omega)}
#' \deqn{\mu_{id} = \mathbf{x}_i^\top\beta_d + u_{\text{upper}(i),\,d}
#'   \;[+\; v_{\text{lower}(i),\,d}]}
#' where \eqn{\Omega} is the residual correlation matrix across classes
#' estimated via an LKJ prior.
#'
#' \strong{Fixed effects} (\code{fixed_effects}): event-level covariates
#' supplied by the user — e.g. age, ICU/ward, HAI/CAI, specimen type, year.
#' All named columns must be present in \code{event_class_data}; the function
#' stops with an informative error if any are missing.
#'
#' \strong{Random effects} (\code{random_effects}): 1 or 2 grouping column
#' names. The first defines the upper clustering level (typically hospital or
#' site); the second, if supplied, defines the lower level (typically patient
#' or admission). Isolate-level residual correlation is captured by \eqn{\Omega}
#' and does not need a third random-effect level. Supplying more than 2 columns
#' raises an error.
#'
#' \strong{Priors} (\code{prior_config}): passed as Stan data so the model
#' compiles once and prior values can be changed without recompilation. Defaults
#' are weakly informative on the probit scale:
#' \itemize{
#'   \item \code{beta_sd = 1.5}  Normal(0, 1.5) for fixed effects —
#'     \eqn{\Phi(\pm1.5) \approx 7\%\text{–}93\%} resistance.
#'   \item \code{tau_sd = 1.0}   Half-Normal(0, 1) for all random-effect SDs —
#'     allows ~30–35 pp between-group variation.
#'   \item \code{lkj_eta = 2.0}  LKJ(2) mildly favours near-zero correlations.
#' }
#'
#' @param event_class_data Data frame. One row per organism-event. Antibiotic
#'   class columns hold 0 (susceptible), 1 (resistant), or \code{NA} (not
#'   tested). Metadata columns hold covariates and grouping variables.
#' @param class_cols Character vector. Names of the antibiotic class columns.
#'   Required — no default.
#' @param fixed_effects Character vector. Event-level covariate column names
#'   to include as fixed effects. Required — no default. All columns must be
#'   present in \code{event_class_data}.
#' @param random_effects Character vector of length 1 or 2. Grouping column
#'   names. Required — no default. First element is the upper clustering level;
#'   second (optional) is the lower clustering level.
#' @param pathogen_col Character. Column identifying the pathogen. Used for
#'   eligibility filtering and output labelling. Default \code{"pathogen"}.
#' @param eligible_pairs Tibble or \code{NULL}. Hospital x pathogen pairs to
#'   include. \code{NULL} uses all pairs present in the data.
#' @param outcome_col Character or \code{NULL}. Patient outcome column. Only
#'   used downstream in \code{compute_event_profile_probabilities()} to split
#'   R_ALL and R_NF cohorts — does not enter the probit likelihood. Default
#'   \code{NULL}.
#' @param reserve_drug_cols Character vector or \code{NULL}. Class columns to
#'   exclude from the main model; handled by \code{summarize_reserve_drugs()}.
#' @param prior_config Named list. Prior hyperparameters. Any subset of
#'   \code{beta_sd}, \code{tau_sd}, \code{lkj_eta} can be supplied; missing
#'   entries fall back to the defaults listed above.
#' @param chains Integer. MCMC chains. Default \code{4L}.
#' @param iter_warmup Integer. Warmup iterations per chain. Default \code{1000L}.
#' @param iter_sampling Integer. Post-warmup iterations. Default \code{1000L}.
#' @param seed Integer. Random seed. Default \code{123L}.
#' @param show_messages Logical. Print sampling progress. Default \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{cmdstanr::sample()}.
#'
#' @return Named list with elements: \code{draws}, \code{diagnostics},
#'   \code{fit}, \code{data_long}, \code{index_maps}, \code{X_design},
#'   \code{class_cols}, \code{event_metadata}, \code{n_re_levels},
#'   \code{upper_re_col}, \code{lower_re_col}, \code{pathogen_col},
#'   \code{prior_config_used}.
#' @export
fit_bayesian_multivariate_probit <- function(
    event_class_data,
    class_cols,
    fixed_effects,
    random_effects,
    pathogen_col      = "pathogen",
    eligible_pairs    = NULL,
    outcome_col       = NULL,
    reserve_drug_cols = NULL,
    prior_config      = list(),
    chains            = 4L,
    iter_warmup       = 1000L,
    iter_sampling     = 1000L,
    seed              = 123L,
    show_messages     = TRUE,
    ...
) {
  # -- Package checks ---------------------------------------------------------
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop(paste0("Package 'cmdstanr' is required for Pathway 2.\n",
                "  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))"),
         call. = FALSE)
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed automatically with cmdstanr).", call. = FALSE)

  # -- Data checks ------------------------------------------------------------
  if (!is.data.frame(event_class_data) || nrow(event_class_data) == 0L)
    stop("`event_class_data` must be a non-empty data frame.", call. = FALSE)

  # -- Validate class_cols ----------------------------------------------------
  if (missing(class_cols) || length(class_cols) == 0L)
    stop("`class_cols` is required. Supply a character vector of antibiotic class column names.",
         call. = FALSE)
  miss_cls <- setdiff(class_cols, names(event_class_data))
  if (length(miss_cls) > 0L)
    stop(sprintf("class_cols not found in event_class_data: %s",
                 paste(miss_cls, collapse = ", ")), call. = FALSE)

  # -- Validate fixed_effects -------------------------------------------------
  if (missing(fixed_effects) || is.null(fixed_effects) || length(fixed_effects) == 0L)
    stop(paste0("`fixed_effects` is required with no default.\n",
                "  Supply a character vector of covariate column names, e.g.:\n",
                "  fixed_effects = c('Age', 'ICU_or_ward', 'HAI_or_CAI', 'specimen_type', 'year')"),
         call. = FALSE)
  miss_fe <- setdiff(fixed_effects, names(event_class_data))
  if (length(miss_fe) > 0L)
    stop(sprintf("fixed_effects column(s) not found in event_class_data: %s",
                 paste(miss_fe, collapse = ", ")), call. = FALSE)

  # -- Validate random_effects ------------------------------------------------
  if (missing(random_effects) || is.null(random_effects) || length(random_effects) == 0L)
    stop(paste0("`random_effects` is required with no default.\n",
                "  Supply 1 or 2 grouping column names, e.g.:\n",
                "  random_effects = c('hospital', 'patient_id')"),
         call. = FALSE)
  if (length(random_effects) > 2L)
    stop(paste0("`random_effects` accepts 1 or 2 grouping columns. ",
                "Isolate-level residual correlation is captured by Omega; ",
                "a third RE level is not identifiable in this model."),
         call. = FALSE)
  miss_re <- setdiff(random_effects, names(event_class_data))
  if (length(miss_re) > 0L)
    stop(sprintf("random_effects column(s) not found in event_class_data: %s",
                 paste(miss_re, collapse = ", ")), call. = FALSE)

  upper_re_col <- random_effects[1L]
  lower_re_col <- if (length(random_effects) == 2L) random_effects[2L] else NULL
  n_re_levels  <- length(random_effects)

  # -- Validate pathogen_col --------------------------------------------------
  if (!pathogen_col %in% names(event_class_data))
    stop(sprintf("pathogen_col '%s' not found in event_class_data.", pathogen_col), call. = FALSE)

  # -- Resolve prior_config ---------------------------------------------------
  pc <- list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0)
  for (nm in names(prior_config)) pc[[nm]] <- prior_config[[nm]]

  if (!is.numeric(pc$beta_sd) || pc$beta_sd <= 0)
    stop("`prior_config$beta_sd` must be a positive number.", call. = FALSE)
  if (!is.numeric(pc$tau_sd) || pc$tau_sd <= 0)
    stop("`prior_config$tau_sd` must be a positive number.", call. = FALSE)
  if (!is.numeric(pc$lkj_eta) || pc$lkj_eta < 1)
    stop("`prior_config$lkj_eta` must be >= 1.", call. = FALSE)

  message(sprintf(
    "[fit_bayesian_multivariate_probit] Priors: beta ~ N(0, %.2g) | tau ~ HN(0, %.2g) | LKJ(%.2g)",
    pc$beta_sd, pc$tau_sd, pc$lkj_eta))

  # -- Remove reserve drug columns from main model ----------------------------
  if (!is.null(reserve_drug_cols)) {
    class_cols <- setdiff(class_cols, reserve_drug_cols)
    if (length(class_cols) < 2L)
      stop("Fewer than 2 class columns remain after removing reserve_drug_cols.", call. = FALSE)
    message(sprintf("[fit_bayesian_multivariate_probit] %d reserve drug column(s) excluded.",
                    length(reserve_drug_cols)))
  }
  if (length(class_cols) < 2L)
    stop("At least 2 class columns are required for the multivariate probit.", call. = FALSE)

  # -- Filter to eligible hospital-pathogen pairs -----------------------------
  event_data <- event_class_data
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

  # -- Internal event identifier ----------------------------------------------
  if ("event_id" %in% names(event_data)) {
    event_data$.eid <- event_data[["event_id"]]
  } else {
    event_data$.eid <- seq_len(nrow(event_data))
  }
  event_id_col <- ".eid"

  # -- Columns to carry into long format --------------------------------------
  meta_carry <- unique(c(
    event_id_col, upper_re_col, lower_re_col, pathogen_col,
    if (!is.null(outcome_col) && outcome_col %in% names(event_data)) outcome_col,
    fixed_effects
  ))
  meta_carry <- meta_carry[!is.null(meta_carry) & meta_carry %in% names(event_data)]

  # -- Pivot to long format (only observed outcomes) --------------------------
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

  message(sprintf(
    "[fit_bayesian_multivariate_probit] %d observed pairs | %d events | %d classes | %d %s(s)",
    nrow(data_long),
    dplyr::n_distinct(data_long[[event_id_col]]),
    length(class_cols),
    dplyr::n_distinct(data_long[[upper_re_col]]),
    upper_re_col))

  # -- Build integer index maps (1-based for Stan) ----------------------------
  class_levels  <- class_cols
  upper_levels  <- sort(unique(data_long[[upper_re_col]]))

  data_long <- data_long %>%
    dplyr::mutate(
      d_idx = match(.data$antibiotic_class, class_levels),
      h_idx = match(.data[[upper_re_col]],  upper_levels)
    )

  lower_levels <- NULL
  if (n_re_levels == 2L) {
    lower_levels <- sort(unique(data_long[[lower_re_col]]))
    data_long <- data_long %>%
      dplyr::mutate(p_idx = match(.data[[lower_re_col]], lower_levels))
    Pt <- length(lower_levels)

    if (Pt * length(class_levels) > 5000L)
      warning(sprintf(
        "Large lower-RE block: %d %s(s) x %d classes = %d parameters. Sampling may be slow.",
        Pt, lower_re_col, length(class_levels), Pt * length(class_levels)), call. = FALSE)
  }

  # -- Fixed-effects design matrix (intercept always included) ----------------
  fe_df <- data_long[, fixed_effects, drop = FALSE]
  for (cc in fixed_effects) {
    if (is.character(fe_df[[cc]])) fe_df[[cc]] <- factor(fe_df[[cc]])
    if (anyNA(fe_df[[cc]])) {
      if (is.factor(fe_df[[cc]])) {
        mode_v <- names(sort(table(fe_df[[cc]]), decreasing = TRUE))[1L]
        fe_df[[cc]][is.na(fe_df[[cc]])] <- mode_v
        warning(sprintf("NAs in fixed effect '%s' imputed with mode ('%s').", cc, mode_v),
                call. = FALSE)
      } else {
        med_v <- stats::median(fe_df[[cc]], na.rm = TRUE)
        fe_df[[cc]][is.na(fe_df[[cc]])] <- med_v
        warning(sprintf("NAs in fixed effect '%s' imputed with median (%g).", cc, med_v),
                call. = FALSE)
      }
    }
  }
  X <- stats::model.matrix(~ ., data = fe_df)

  D <- length(class_levels)
  H <- length(upper_levels)
  K <- ncol(X)

  # -- Build Stan data list ---------------------------------------------------
  stan_data <- list(
    N             = nrow(data_long),
    D             = D,
    H             = H,
    K             = K,
    y             = as.integer(data_long$resistance_binary),
    X             = unname(X),
    d_idx         = as.integer(data_long$d_idx),
    h_idx         = as.integer(data_long$h_idx),
    prior_beta_sd = as.numeric(pc$beta_sd),
    prior_tau_sd  = as.numeric(pc$tau_sd),
    lkj_eta       = as.numeric(pc$lkj_eta)
  )
  if (n_re_levels == 2L) {
    stan_data$Pt    <- as.integer(Pt)
    stan_data$p_idx <- as.integer(data_long$p_idx)
  }

  # -- Compile Stan model (cached by cmdstanr via content hash) ---------------
  stan_code <- if (n_re_levels == 1L) .amr_probit_stan_1re() else .amr_probit_stan_2re()
  stan_file  <- cmdstanr::write_stan_file(stan_code)
  message(sprintf("[fit_bayesian_multivariate_probit] Compiling %d-RE Stan model...",
                  n_re_levels))
  mod <- cmdstanr::cmdstan_model(stan_file, compile = TRUE)

  # -- Sample -----------------------------------------------------------------
  message(sprintf(
    "[fit_bayesian_multivariate_probit] Sampling: %d chains x (%d warmup + %d sampling)...",
    chains, iter_warmup, iter_sampling))
  fit <- mod$sample(
    data          = stan_data,
    seed          = as.integer(seed),
    chains        = as.integer(chains),
    iter_warmup   = as.integer(iter_warmup),
    iter_sampling = as.integer(iter_sampling),
    refresh       = if (show_messages) 200L else 0L,
    ...
  )

  # -- Diagnostics ------------------------------------------------------------
  draws_arr    <- fit$draws(format = "draws_array")
  rhat_tbl     <- fit$summary(variables = NULL, "rhat")
  ess_tbl      <- fit$summary(variables = NULL, "ess_bulk", "ess_tail")
  n_divergent  <- sum(fit$sampler_diagnostics(format = "matrix")[, "divergent__"])
  max_rhat     <- max(rhat_tbl$rhat,    na.rm = TRUE)
  min_ess_bulk <- min(ess_tbl$ess_bulk, na.rm = TRUE)

  if (max_rhat > 1.05)
    warning(sprintf("Convergence concern: max Rhat = %.3f (> 1.05).", max_rhat), call. = FALSE)
  if (n_divergent > 0L)
    warning(sprintf("%d divergent transition(s). Consider higher adapt_delta.", n_divergent),
            call. = FALSE)

  diag_tbl <- tibble::tibble(
    n_chains          = as.integer(chains),
    iter_warmup       = as.integer(iter_warmup),
    iter_sampling     = as.integer(iter_sampling),
    n_re_levels       = as.integer(n_re_levels),
    n_observed_pairs  = nrow(data_long),
    n_classes         = D,
    n_upper_groups    = H,
    n_lower_groups    = if (n_re_levels == 2L) as.integer(Pt) else NA_integer_,
    n_divergent       = as.integer(n_divergent),
    max_rhat          = round(max_rhat, 4L),
    min_ess_bulk      = round(min_ess_bulk, 1L)
  )

  message(sprintf(
    "[fit_bayesian_multivariate_probit] Done. max_Rhat=%.3f | min_ESS=%.0f | divergent=%d",
    max_rhat, min_ess_bulk, n_divergent))

  list(
    draws            = draws_arr,
    diagnostics      = diag_tbl,
    fit              = fit,
    data_long        = tibble::as_tibble(data_long),
    index_maps       = list(
      class_levels  = class_levels,
      upper_levels  = upper_levels,
      lower_levels  = lower_levels
    ),
    X_design         = X,
    class_cols       = class_cols,
    event_metadata   = tibble::as_tibble(event_data),
    n_re_levels      = n_re_levels,
    upper_re_col     = upper_re_col,
    lower_re_col     = lower_re_col,
    pathogen_col     = pathogen_col,
    prior_config_used = pc
  )
}


# ---------------------------------------------------------------------------
# compute_event_profile_probabilities()
# ---------------------------------------------------------------------------

#' Compute Posterior Resistance Profile Probabilities via MVN Simulation
#'
#' For each posterior draw, constructs the linear predictor \eqn{\mu_i} for
#' every eligible event using the fixed effects and random effects from the
#' fitted model, simulates latent \eqn{Z_i \sim \text{MVN}(\mu_i, \Omega)}
#' over the hospital-specific eligible class set, converts to binary outcomes,
#' and accumulates profile probabilities. Returns both event-level posterior
#' means and draw-level aggregate R_ALL / R_NF values per hospital x pathogen
#' pair for computing proper Bayesian credible intervals downstream.
#'
#' @param fitted_model List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param n_profile_simulations Integer. Number of posterior draws to use.
#'   Subsampled without replacement when total draws exceed this value.
#'   Default \code{2000L}.
#' @param outcome_col Character or \code{NULL}. Patient outcome column in
#'   \code{fitted_model$event_metadata}. When \code{NULL}, all events are
#'   treated as having a known outcome and R_NF is not computed separately.
#' @param nonfatal_values Character vector. Outcome values that define the
#'   non-fatal cohort for R_NF. Default covers common discharge/survival labels.
#' @param seed Integer. Random seed for MVN simulation. Default \code{123L}.
#'
#' @return Named list: \code{event_profiles} (event-level posterior means) and
#'   \code{aggregate_draws} (per-draw R_ALL and R_NF per hospital x pathogen x
#'   profile, used by \code{aggregate_profiles_for_daly()} for credible
#'   intervals).
#' @export
compute_event_profile_probabilities <- function(
    fitted_model,
    n_profile_simulations = 2000L,
    outcome_col           = NULL,
    nonfatal_values       = c("Discharged", "Survived", "Alive",
                              "discharged", "survived", "alive"),
    seed                  = 123L
) {
  if (!requireNamespace("mvtnorm",   quietly = TRUE))
    stop("Package 'mvtnorm' is required. Install with: install.packages('mvtnorm')", call. = FALSE)
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)

  set.seed(as.integer(seed))

  draws        <- fitted_model$draws
  class_cols   <- fitted_model$class_cols
  idx_maps     <- fitted_model$index_maps
  X_long       <- fitted_model$X_design
  data_long    <- fitted_model$data_long
  event_meta   <- fitted_model$event_metadata
  n_re_levels  <- fitted_model$n_re_levels
  upper_re_col <- fitted_model$upper_re_col
  lower_re_col <- fitted_model$lower_re_col
  pathogen_col <- fitted_model$pathogen_col

  D <- length(idx_maps$class_levels)
  H <- length(idx_maps$upper_levels)
  K <- ncol(X_long)

  # -- Thin draws to n_profile_simulations ------------------------------------
  draws_mat <- posterior::as_draws_matrix(draws)
  n_total   <- nrow(draws_mat)
  S         <- min(as.integer(n_profile_simulations), n_total)
  draw_idx  <- if (S < n_total) sort(sample.int(n_total, S)) else seq_len(n_total)
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  # -- Extract parameter arrays [S, dim1, dim2] -------------------------------
  .arr <- function(prefix, d1, d2) {
    cols <- as.vector(outer(seq_len(d1), seq_len(d2),
                            function(a, b) sprintf("%s[%d,%d]", prefix, a, b)))
    array(draws_mat[, cols, drop = FALSE], dim = c(S, d1, d2))
  }
  beta_arr  <- .arr("beta",  K, D)
  u_arr     <- .arr("u",     H, D)
  omega_arr <- .arr("Omega", D, D)

  has_lower_re <- n_re_levels == 2L
  if (has_lower_re) {
    Pt    <- length(idx_maps$lower_levels)
    v_arr <- .arr("v", Pt, D)
  }

  # -- Event-level design and index vectors -----------------------------------
  event_id_col <- if ("event_id" %in% names(data_long)) "event_id" else ".eid"
  event_rows   <- data_long %>%
    dplyr::group_by(.data[[event_id_col]]) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()

  row_match  <- match(event_rows[[event_id_col]], data_long[[event_id_col]])
  X_event    <- X_long[row_match, , drop = FALSE]
  h_events   <- as.integer(event_rows$h_idx)
  p_events   <- if (has_lower_re) as.integer(event_rows$p_idx) else NULL
  N_ev       <- nrow(event_rows)

  # -- Outcome flags per event ------------------------------------------------
  if (!is.null(outcome_col) && outcome_col %in% names(event_rows)) {
    ov               <- event_rows[[outcome_col]]
    is_known_outcome <- !is.na(ov)
    is_nonfatal      <- !is.na(ov) & (ov %in% nonfatal_values)
  } else {
    is_known_outcome <- rep(TRUE, N_ev)
    is_nonfatal      <- rep(TRUE, N_ev)
  }

  # -- Hospital-pathogen eligible class sets ----------------------------------
  hp_pairs <- unique(event_meta[, c(upper_re_col, pathogen_col), drop = FALSE])
  hp_keys  <- paste(hp_pairs[[upper_re_col]], hp_pairs[[pathogen_col]], sep = "‖")

  hp_class_d <- stats::setNames(lapply(seq_len(nrow(hp_pairs)), function(r) {
    h_nm <- hp_pairs[[upper_re_col]][r]
    k_nm <- hp_pairs[[pathogen_col]][r]
    sub  <- event_meta[event_meta[[upper_re_col]] == h_nm &
                         event_meta[[pathogen_col]]  == k_nm, class_cols, drop = FALSE]
    which(vapply(class_cols, function(cc) any(!is.na(sub[[cc]])), logical(1L)))
  }), hp_keys)

  hp_ev_idx <- stats::setNames(lapply(hp_keys, function(key) {
    parts <- strsplit(key, "‖", fixed = TRUE)[[1L]]
    which(event_meta[[upper_re_col]] == parts[1L] &
            event_meta[[pathogen_col]]  == parts[2L])
  }), hp_keys)

  # Storage: per-hp, per-draw profile labels matrix [N_hp × S]
  hp_store <- stats::setNames(lapply(hp_keys, function(key) {
    n_hp <- length(hp_ev_idx[[key]])
    list(labels    = matrix(NA_character_, n_hp, S),
         known_out = is_known_outcome[hp_ev_idx[[key]]],
         nonfatal  = is_nonfatal[hp_ev_idx[[key]]])
  }), hp_keys)

  message(sprintf(
    "[compute_event_profile_probabilities] %d draws | %d events | %d hp-pairs | %d RE level(s)",
    S, N_ev, length(hp_keys), n_re_levels))

  # -- Main simulation loop (vectorised over events within each hp pair) ------
  for (s in seq_len(S)) {
    beta_s  <- matrix(beta_arr[s, , ], nrow = K, ncol = D)
    u_s     <- matrix(u_arr[s, , ],    nrow = H, ncol = D)
    Omega_s <- matrix(omega_arr[s, , ], nrow = D, ncol = D)

    mu_all <- X_event %*% beta_s + u_s[h_events, , drop = FALSE]
    if (has_lower_re) {
      v_s     <- matrix(v_arr[s, , ], nrow = Pt, ncol = D)
      mu_all  <- mu_all + v_s[p_events, , drop = FALSE]
    }

    for (key in hp_keys) {
      d_hp   <- hp_class_d[[key]]
      ev_idx <- hp_ev_idx[[key]]
      if (length(ev_idx) == 0L || length(d_hp) < 1L) next

      mu_hp    <- mu_all[ev_idx, d_hp, drop = FALSE]
      Omega_hp <- Omega_s[d_hp, d_hp, drop = FALSE]
      Omega_hp <- (Omega_hp + t(Omega_hp)) / 2 + diag(1e-9, length(d_hp))

      L_hp <- tryCatch(t(chol(Omega_hp)), error = function(e) NULL)
      if (is.null(L_hp)) next

      Z_std <- matrix(stats::rnorm(nrow(mu_hp) * length(d_hp)), nrow = nrow(mu_hp))
      Z_hp  <- mu_hp + Z_std %*% t(L_hp)

      hp_store[[key]]$labels[, s] <- apply(
        (Z_hp > 0) * 1L, 1L,
        function(row) paste(ifelse(row == 1L, "R", "S"), collapse = "")
      )
    }
  }

  # -- Compute outputs --------------------------------------------------------
  event_profile_rows  <- list()
  aggregate_draw_rows <- list()

  for (key in hp_keys) {
    parts   <- strsplit(key, "‖", fixed = TRUE)[[1L]]
    h_nm    <- parts[1L]; k_nm <- parts[2L]
    d_hp    <- hp_class_d[[key]]
    ev_idx  <- hp_ev_idx[[key]]
    store   <- hp_store[[key]]
    n_hp    <- nrow(store$labels)
    if (n_hp == 0L) next

    valid_s <- which(colSums(!is.na(store$labels)) == n_hp)
    if (length(valid_s) == 0L) next
    lbl_v   <- store$labels[, valid_s, drop = FALSE]

    class_set <- paste(class_cols[d_hp], collapse = "|")
    all_lbls  <- sort(unique(as.vector(lbl_v)))

    # Event-level posterior means
    for (ev_i in seq_len(n_hp)) {
      ev_lbl <- lbl_v[ev_i, ]
      for (lbl in all_lbls) {
        event_profile_rows[[length(event_profile_rows) + 1L]] <- tibble::tibble(
          !!upper_re_col   := h_nm,
          !!pathogen_col   := k_nm,
          event_id          = event_rows[[event_id_col]][ev_idx[ev_i]],
          profile_class_set = class_set,
          profile_delta     = lbl,
          profile_probability = mean(ev_lbl == lbl, na.rm = TRUE)
        )
      }
    }

    # Draw-level aggregates
    for (s in seq_along(valid_s)) {
      s_lbl <- lbl_v[, s]
      for (lbl in all_lbls) {
        aggregate_draw_rows[[length(aggregate_draw_rows) + 1L]] <- tibble::tibble(
          !!upper_re_col   := h_nm,
          !!pathogen_col   := k_nm,
          profile_class_set = class_set,
          profile_delta     = lbl,
          draw_s            = valid_s[s],
          R_ALL_s = if (any(store$known_out)) mean(s_lbl[store$known_out] == lbl, na.rm = TRUE) else NA_real_,
          R_NF_s  = if (any(store$nonfatal))  mean(s_lbl[store$nonfatal]  == lbl, na.rm = TRUE) else NA_real_
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

aggregate_profiles_for_daly <- function(
    profile_output,
    event_class_data  = NULL,
    hospital_col      = "hospital",
    patient_col       = "patient_id",
    admission_col     = "admission_id",
    pathogen_col      = "pathogen",
    weighting         = c("equal_event", "patient_admission"),
    ci_level          = 0.95
) {
  weighting <- match.arg(weighting)
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

  # Compute per-event admission weights if requested
  if (weighting == "patient_admission") {
    if (is.null(event_class_data))
      stop("`event_class_data` must be supplied when weighting = 'patient_admission'.",
           call. = FALSE)
    adm_present <- !is.null(admission_col) && admission_col %in% names(event_class_data)
    adm_key_col <- if (adm_present) admission_col else patient_col
    # For each admission, weight = 1 / n_events_in_admission — applied in aggregate_draws
    # (draw-level R values already aggregate over events; re-weighting would require
    #  event-level draw storage; here we apply equal_event and note the limitation)
    warning(paste0("patient_admission weighting at the draw level requires event-level draw storage. ",
                   "Falling back to equal_event weighting. ",
                   "To use patient_admission weighting, re-aggregate from event_profiles directly."),
            call. = FALSE)
    weighting <- "equal_event"
  }

  # Summarise draw-level R_ALL and R_NF to posterior mean + CI
  result <- agg_draws %>%
    dplyr::group_by(
      .data[[hospital_col]],
      .data[[pathogen_col]],
      .data$profile_class_set,
      .data$profile_delta
    ) %>%
    dplyr::summarise(
      R_ALL_mean  = mean(.data$R_ALL_s,  na.rm = TRUE),
      R_ALL_lower = stats::quantile(.data$R_ALL_s, lo_q, na.rm = TRUE),
      R_ALL_upper = stats::quantile(.data$R_ALL_s, hi_q, na.rm = TRUE),
      R_NF_mean   = mean(.data$R_NF_s,   na.rm = TRUE),
      R_NF_lower  = stats::quantile(.data$R_NF_s,  lo_q, na.rm = TRUE),
      R_NF_upper  = stats::quantile(.data$R_NF_s,  hi_q, na.rm = TRUE),
      n_draws_used = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      profile_label        = .data$profile_delta,
      profile_set_type     = "facility_bayesian_probit",
      used_for_YLL         = TRUE,
      used_for_YLD         = TRUE,
      estimator            = "bayesian_multivariate_probit",
      identifiability_flag = n_draws_used < 100L,
      dplyr::across(dplyr::starts_with("R_"), ~ round(.x, 6L))
    ) %>%
    dplyr::select(
      dplyr::all_of(c(hospital_col, pathogen_col)),
      profile_set_type, profile_class_set, profile_delta, profile_label,
      R_ALL_mean, R_ALL_lower, R_ALL_upper,
      R_NF_mean,  R_NF_lower,  R_NF_upper,
      used_for_YLL, used_for_YLD, estimator, identifiability_flag
    )

  message(sprintf(
    "[aggregate_profiles_for_daly] %d hospital-pathogen-profile rows | weighting: %s",
    nrow(result), weighting))

  result
}


# ---------------------------------------------------------------------------
# estimate_resistance_profiles()  — top-level dispatcher for both pathways
# ---------------------------------------------------------------------------

#' Estimate Resistance Profiles: Pathway 1 (Convex) or Pathway 2 (Bayesian)
#'
#' Top-level dispatcher. Runs either the convex optimisation pathway
#' (Pathway 1, aggregate surveillance data) or the Bayesian hierarchical
#' multivariate probit pathway (Pathway 2, facility-level AST data) and
#' returns a standardised named list compatible with the DALY burden pipeline.
#'
#' The caller is responsible for all upstream decisions: which antibiotic
#' classes to include for each hospital-pathogen pair, which events are
#' eligible, and which events to exclude. These are dataset-specific analysis
#' choices that belong in the calling script, not in this function.
#' This function receives a pre-prepared dataset and performs estimation only.
#'
#' @param data For Pathway 1: marginals tibble from
#'   \code{compute_marginals_from_data()} or a validated aggregate table.
#'   For Pathway 2: wide-format event-level tibble. One row per organism-event;
#'   antibiotic class columns hold 0 (susceptible), 1 (resistant), or
#'   \code{NA} (not tested). Only eligible classes and events should be
#'   present — subsetting is the caller's responsibility.
#' @param method Character. \code{"convex"} or \code{"bayesian"}.
#' @param pairwise Tibble or \code{NULL}. Pathway 1 only. Pairwise
#'   co-resistance rates from \code{compute_pairwise_from_data()}.
#' @param panel_map Named list or \code{NULL}. Pathway 1 only.
#' @param class_cols Character vector. Pathway 2 only. Names of the antibiotic
#'   class columns in \code{data}. Required — no default.
#' @param fixed_effects Character vector. Pathway 2 only. Event-level covariate
#'   column names (e.g. age, ICU/ward, HAI/CAI, specimen type, year).
#'   Required — no default.
#' @param random_effects Character vector (1 or 2 elements). Pathway 2 only.
#'   Grouping column names (e.g. \code{c("hospital", "patient_id")}).
#'   Required — no default.
#' @param pathogen_col Character. Column identifying the pathogen.
#'   Default \code{"pathogen"}.
#' @param eligible_pairs Tibble or \code{NULL}. Pathway 2 only. When supplied,
#'   restricts fitting to these hospital x pathogen combinations. Must contain
#'   the upper random-effect column and \code{pathogen_col}.
#' @param outcome_col Character or \code{NULL}. Pathway 2 only. Patient outcome
#'   column used to split R_ALL and R_NF cohorts. Does not enter the
#'   likelihood. Default \code{NULL}.
#' @param nonfatal_values Character vector. Outcome values defining the
#'   non-fatal cohort (R_NF). Ignored when \code{outcome_col = NULL}.
#' @param prior_config Named list. Pathway 2 only. Any subset of
#'   \code{beta_sd} (default 1.5), \code{tau_sd} (default 1.0),
#'   \code{lkj_eta} (default 2.0). Missing entries use defaults.
#' @param weighting Character. Pathway 2 aggregation weighting.
#'   \code{"equal_event"} (default) or \code{"patient_admission"}.
#' @param n_profile_simulations Integer. Pathway 2 posterior draws used for
#'   MVN simulation. Default \code{2000L}.
#' @param chains Integer. MCMC chains. Default \code{4L}.
#' @param iter_warmup Integer. Warmup iterations per chain. Default \code{1000L}.
#' @param iter_sampling Integer. Post-warmup iterations. Default \code{1000L}.
#' @param seed Integer. Random seed. Default \code{123L}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{profiles}}{Profile probability tibble (Pathway 1) or
#'     R_ALL / R_NF tibble (Pathway 2).}
#'   \item{\code{eligibility}}{Hospital-level summary when
#'     \code{eligible_pairs} is supplied; \code{NULL} otherwise.}
#'   \item{\code{diagnostics}}{Constraint residual summary (Pathway 1) or
#'     MCMC convergence diagnostics (Pathway 2).}
#'   \item{\code{fitted_models}}{The raw fitted model object(s).}
#'   \item{\code{config_used}}{Named list of all resolved configuration values.}
#' }
#' @export
estimate_resistance_profiles <- function(
    data,
    method                = c("convex", "bayesian"),
    # Pathway 1
    pairwise              = NULL,
    panel_map             = NULL,
    # Pathway 2 — required with no defaults
    class_cols            = NULL,
    fixed_effects         = NULL,
    random_effects        = NULL,
    pathogen_col          = "pathogen",
    eligible_pairs        = NULL,
    outcome_col           = NULL,
    nonfatal_values       = c("Discharged", "Survived", "Alive",
                              "discharged", "survived", "alive"),
    prior_config          = list(),
    weighting             = c("equal_event", "patient_admission"),
    n_profile_simulations = 2000L,
    chains                = 4L,
    iter_warmup           = 1000L,
    iter_sampling         = 1000L,
    seed                  = 123L
) {
  method    <- match.arg(method)
  weighting <- match.arg(weighting)

  config_used <- list(
    method                = method,
    pathogen_col          = pathogen_col,
    seed                  = seed,
    n_profile_simulations = n_profile_simulations,
    weighting             = weighting,
    chains                = chains,
    iter_warmup           = iter_warmup,
    iter_sampling         = iter_sampling,
    prior_config          = prior_config
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

  if (is.null(class_cols))
    stop("`class_cols` is required for method = 'bayesian'.", call. = FALSE)
  if (is.null(fixed_effects))
    stop(paste0("`fixed_effects` is required for method = 'bayesian'.\n",
                "  Supply a character vector of covariate column names."), call. = FALSE)
  if (is.null(random_effects))
    stop(paste0("`random_effects` is required for method = 'bayesian'.\n",
                "  Supply 1 or 2 grouping column names."), call. = FALSE)

  fitted_mod <- fit_bayesian_multivariate_probit(
    event_class_data = data,
    class_cols       = class_cols,
    fixed_effects    = fixed_effects,
    random_effects   = random_effects,
    pathogen_col     = pathogen_col,
    eligible_pairs   = eligible_pairs,
    outcome_col      = outcome_col,
    prior_config     = prior_config,
    chains           = chains,
    iter_warmup      = iter_warmup,
    iter_sampling    = iter_sampling,
    seed             = seed
  )

  profile_probs <- compute_event_profile_probabilities(
    fitted_model          = fitted_mod,
    n_profile_simulations = as.integer(n_profile_simulations),
    outcome_col           = outcome_col,
    nonfatal_values       = nonfatal_values,
    seed                  = seed
  )

  profiles_tbl <- aggregate_profiles_for_daly(
    profile_output   = profile_probs,
    event_class_data = data,
    hospital_col     = fitted_mod$upper_re_col,
    patient_col      = if (!is.null(fitted_mod$lower_re_col))
                         fitted_mod$lower_re_col
                       else
                         fitted_mod$upper_re_col,
    pathogen_col     = pathogen_col,
    weighting        = weighting
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
