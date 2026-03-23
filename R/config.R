# config.R
# Configuration system for AMR preprocessing pipeline

#' Create AMR Preprocessing Configuration
#'
#' Creates a configuration object that controls the behavior of the
#' preprocessing pipeline. Users can override default settings for their
#' specific use case.
#'
#' @param column_mappings Named list of column name mappings. Default uses
#'   \code{default_column_mappings}.
#' @param fuzzy_match Logical. Enable fuzzy matching for column names.
#'   Default TRUE.
#' @param strict_validation Logical. If TRUE, stops execution if required
#'   fields are missing. If FALSE, issues warnings. Default FALSE.
#' @param date_columns Character vector of date column names to parse.
#'   Default includes standard date fields.
#' @param hai_cutoff Numeric. Number of days after admission to classify
#'   as Hospital-Acquired Infection (HAI). Default 2.
#' @param infer_department Logical. Attempt to infer hospital department
#'   from other variables if missing. Default TRUE.
#' @param event_gap_days Numeric. Minimum days between events to create
#'   new event_id. Default 14.
#' @param mortality_window Numeric. Days after culture to classify death
#'   as infection-related. Default 14.
#' @param age_bins Character. Age binning strategy: "GBD_standard",
#'   "pediatric", or "geriatric". Can also be a custom vector. Default
#'   "GBD_standard".
#' @param contaminant_method Character. Method for contaminant classification:
#'   "auto" (cascade through available methods), "device_based",
#'   "heuristic", "provided". Default "auto".
#' @param mdr_definition Character or numeric. MDR definition: "CDC", "WHO",
#'   or numeric threshold for number of resistant classes. Default "CDC".
#' @param xdr_definition Character or numeric. XDR definition. Default "CDC".
#' @param map_icd10 Logical. Attempt ICD-10 code mapping. Default TRUE.
#' @param rr_table Character. Name of built-in RR table to use, or path
#'   to custom RR table. Default "GBD_2021".
#' @param organism_map Character or named vector. "default" for built-in
#'   mapping, or custom named vector. Default "default".
#' @param antibiotic_map Character. Built-in antibiotic classification:
#'   "WHO_2023" or path to custom. Default "WHO_2023".
#' @param intermediate_as_resistant Logical. Treat Intermediate (I) as
#'   Resistant (R) except for special cases like Colistin. Default TRUE.
#' @param verbose Logical. Print progress messages during processing.
#'   Default TRUE.
#'
#' @return An \code{amr_config} object (list with class)
#'
#' @export
#' @examples
#' \dontrun{
#' # Use defaults
#' config <- amr_config()
#'
#' # Customize for specific hospital
#' config <- amr_config(
#'   hai_cutoff = 3,
#'   mdr_definition = 5,
#'   strict_validation = FALSE
#' )
#' }
amr_config <- function(column_mappings = NULL,
                       fuzzy_match = TRUE,
                       strict_validation = FALSE,
                       date_columns = c(
                         "date_of_admission", "date_of_culture",
                         "date_of_final_outcome", "DOB"
                       ),
                       hai_cutoff = 2,
                       infer_department = TRUE,
                       event_gap_days = 14,
                       mortality_window = 14,
                       age_bins = "GBD_standard",
                       contaminant_method = "auto",
                       mdr_definition = "CDC",
                       xdr_definition = "CDC",
                       map_icd10 = TRUE,
                       rr_table = "GBD_2021",
                       organism_map = "default",
                       antibiotic_map = "WHO_2023",
                       intermediate_as_resistant = TRUE,
                       verbose = TRUE) {
  # Use default column mappings if none provided
  if (is.null(column_mappings)) {
    column_mappings <- default_column_mappings
  }

  # Convert age_bins if character
  if (is.character(age_bins) && length(age_bins) == 1) {
    age_bins <- get_age_bins(age_bins)
  }

  # Create config list
  config <- list(
    # Phase 1: Standardization
    column_mappings = column_mappings,
    fuzzy_match = fuzzy_match,
    strict_validation = strict_validation,
    date_columns = date_columns,

    # Phase 2: Enrichment
    hai_cutoff = hai_cutoff,
    infer_department = infer_department,

    # Phase 3: Derivation
    event_gap_days = event_gap_days,
    mortality_window = mortality_window,
    age_bins = age_bins,
    contaminant_method = contaminant_method,
    mdr_definition = mdr_definition,
    xdr_definition = xdr_definition,
    map_icd10 = map_icd10,

    # Reference data
    rr_table = rr_table,
    organism_map = organism_map,
    antibiotic_map = antibiotic_map,

    # Processing rules
    intermediate_as_resistant = intermediate_as_resistant,
    verbose = verbose
  )

  # Set class
  class(config) <- c("amr_config", "list")

  return(config)
}


#' Print AMR Configuration
#'
#' @param x An amr_config object
#' @param ... Additional arguments (not used)
#'
#' @export
print.amr_config <- function(x, ...) {
  cat("AMR Preprocessing Configuration\n")
  cat("================================\n\n")

  cat("Phase 1 - Standardization:\n")
  cat(sprintf("  Fuzzy column matching: %s\n", x$fuzzy_match))
  cat(sprintf("  Strict validation: %s\n", x$strict_validation))
  cat(sprintf("  Date columns: %s\n", paste(x$date_columns, collapse = ", ")))

  cat("\nPhase 2 - Enrichment:\n")
  cat(sprintf("  HAI cutoff: %d days\n", x$hai_cutoff))
  cat(sprintf("  Infer department: %s\n", x$infer_department))

  cat("\nPhase 3 - Derivation:\n")
  cat(sprintf("  Event gap: %d days\n", x$event_gap_days))
  cat(sprintf("  Mortality window: %d days\n", x$mortality_window))
  cat(sprintf(
    "  Age bins: %s\n",
    if (is.character(x$age_bins)) "custom" else paste(length(x$age_bins), "bins")
  ))
  cat(sprintf("  Contaminant method: %s\n", x$contaminant_method))
  cat(sprintf("  MDR definition: %s\n", x$mdr_definition))
  cat(sprintf("  XDR definition: %s\n", x$xdr_definition))

  cat("\nReference Data:\n")
  cat(sprintf("  RR table: %s\n", x$rr_table))
  cat(sprintf("  Organism map: %s\n", x$organism_map))
  cat(sprintf("  Antibiotic map: %s\n", x$antibiotic_map))

  cat("\nProcessing Rules:\n")
  cat(sprintf("  Intermediate as Resistant: %s\n", x$intermediate_as_resistant))
  cat(sprintf("  Verbose: %s\n", x$verbose))

  invisible(x)
}


#' Validate AMR Configuration
#'
#' Checks that configuration parameters are valid.
#'
#' @param config An amr_config object
#'
#' @return Logical. TRUE if valid, stops with error if invalid.
#' @export
validate_config <- function(config) {
  if (!inherits(config, "amr_config")) {
    stop("config must be an amr_config object created with amr_config()")
  }

  # Check numeric parameters
  if (config$hai_cutoff < 0) {
    stop("hai_cutoff must be >= 0")
  }

  if (config$event_gap_days < 1) {
    stop("event_gap_days must be >= 1")
  }

  if (config$mortality_window < 1) {
    stop("mortality_window must be >= 1")
  }

  # Check contaminant method
  valid_methods <- c("auto", "device_based", "heuristic", "provided")
  if (!config$contaminant_method %in% valid_methods) {
    stop(sprintf(
      "contaminant_method must be one of: %s",
      paste(valid_methods, collapse = ", ")
    ))
  }

  # Check MDR definition
  if (is.character(config$mdr_definition)) {
    if (!config$mdr_definition %in% c("CDC", "WHO")) {
      stop("mdr_definition must be 'CDC', 'WHO', or a numeric threshold")
    }
  } else if (is.numeric(config$mdr_definition)) {
    if (config$mdr_definition < 1) {
      stop("mdr_definition (numeric) must be >= 1")
    }
  } else {
    stop("mdr_definition must be character or numeric")
  }

  # Check age bins
  if (!is.character(config$age_bins) || length(config$age_bins) < 2) {
    stop("age_bins must be a character vector with at least 2 bins")
  }

  TRUE
}
