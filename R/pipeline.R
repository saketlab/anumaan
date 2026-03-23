# pipeline.R
# Master preprocessing pipeline

#' AMR Data Preprocessing Pipeline
#'
#' Master function that orchestrates the complete AMR data preprocessing
#' pipeline. Applies standardization, enrichment, and derivation in sequence.
#'
#' Pipeline phases:
#' 1. **Standardization**: Column mapping, value normalization, date parsing
#' 2. **Enrichment**: Derive missing optional variables (Age, LOS, infection type)
#' 3. **Derivation**: Create analytical variables (event IDs, MDR/XDR, weights)
#'
#' @param data Data frame. Raw AMR dataset.
#' @param config Configuration object from amr_config(). If NULL, uses defaults.
#' @param phases Character vector. Which phases to run: "standardize", "enrich",
#'   "derive", or "all" (default). Allows partial pipeline execution.
#' @param verbose Logical. Print detailed progress messages. Default TRUE.
#' @param validate Logical. Run validation checks before and after. Default TRUE.
#' @param generate_report Logical. Generate preprocessing report. Default TRUE.
#'
#' @return List with class "amr_result":
#'   - data: Preprocessed data frame
#'   - config: Configuration used
#'   - log: List of processing logs and summaries
#'   - report: Preprocessing report (if generate_report = TRUE)
#'   - metadata: Pipeline execution metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Full pipeline with defaults
#' result <- amr_preprocess(raw_data)
#' clean_data <- result$data
#'
#' # Custom configuration
#' config <- amr_config(
#'   hai_cutoff = 3,
#'   mdr_definition = "CDC",
#'   fuzzy_match = TRUE
#' )
#' result <- amr_preprocess(raw_data, config = config)
#'
#' # Run only standardization phase
#' result <- amr_preprocess(raw_data, phases = "standardize")
#'
#' # Summary
#' summary(result)
#' }
amr_preprocess <- function(data,
                           config = NULL,
                           phases = "all",
                           verbose = TRUE,
                           validate = TRUE,
                           generate_report = TRUE) {
  # Start timer
  start_time <- Sys.time()

  if (verbose) {
    message("═══════════════════════════════════════════════════════")
    message("  anumaan Preprocessing Pipeline")
    message("═══════════════════════════════════════════════════════")
    message(sprintf("Started: %s", format(start_time, "%Y-%m-%d %H:%M:%S")))
    message(sprintf("Input: %d rows × %d columns", nrow(data), ncol(data)))
  }

  # Initialize config
  if (is.null(config)) {
    config <- amr_config()
    if (verbose) message("\nUsing default configuration")
  }

  # Expand "all" phases
  all_phases <- c("standardize", "enrich", "derive")
  if ("all" %in% phases) {
    phases <- all_phases
  }

  # Initialize log
  log <- list(
    phases_run = phases,
    column_mapping = NULL,
    transformations = NULL,
    data_quality = NULL,
    validation = NULL
  )

  # Store original data summary
  original_summary <- summarize_raw_data(data)

  # ========================================
  # PRE-VALIDATION
  # ========================================
  if (validate && verbose) {
    message("\n─── Pre-Processing Validation ───")
  }

  if (validate) {
    # Basic quality check
    quality <- validate_data_quality(
      data,
      min_rows = 1,
      max_missing_pct = 90,
      required_cols = character(0), # Don't require any columns yet
      stop_on_failure = FALSE
    )

    log$validation$pre_quality <- quality

    # Remove exact duplicates
    n_before_dedup <- nrow(data)
    data <- remove_duplicate_rows(data, report = verbose)
    log$validation$duplicates_removed <- n_before_dedup - nrow(data)
  }

  # ========================================
  # PHASE 1: STANDARDIZATION
  # ========================================
  if ("standardize" %in% phases) {
    if (verbose) {
      message("\n═══════════════════════════════════════════════════════")
      message("  PHASE 1: STANDARDIZATION")
      message("═══════════════════════════════════════════════════════")
    }

    # Step 1.1: Column name standardization
    if (verbose) message("\n[1.1] Standardizing column names...")
    mapping_result <- standardize_column_names(
      data,
      fuzzy_match = config$fuzzy_match,
      mapping = config$custom_column_mappings
    )
    data <- mapping_result$data
    log$column_mapping <- mapping_result$mapping_log

    # Step 1.2: Parse dates
    if (verbose) message("\n[1.2] Parsing dates...")
    date_cols <- c("date_of_admission", "date_of_culture", "date_of_final_outcome", "DOB")
    existing_date_cols <- intersect(date_cols, names(data))
    for (col in existing_date_cols) {
      data <- parse_dates(data, date_col = col)
    }

    # Step 1.3: Standardize values
    if (verbose) message("\n[1.3] Standardizing values...")
    if ("gender" %in% names(data)) {
      data <- standardize_gender(data)
    }
    if ("final_outcome" %in% names(data)) {
      data <- standardize_outcome(data)
    }
    if ("antibiotic_value" %in% names(data)) {
      data <- standardize_susceptibility(data)
    }

    # Step 1.4: Normalize organisms and antibiotics
    if (verbose) message("\n[1.4] Normalizing organisms and antibiotics...")
    if ("organism_name" %in% names(data)) {
      data <- normalize_organism(data, organism_col = "organism_name")
    }
    if ("antibiotic_name" %in% names(data)) {
      data <- normalize_antibiotic(data, antibiotic_col = "antibiotic_name")
    }
  }

  # ========================================
  # PHASE 2: ENRICHMENT
  # ========================================
  if ("enrich" %in% phases) {
    if (verbose) {
      message("\n═══════════════════════════════════════════════════════")
      message("  PHASE 2: ENRICHMENT")
      message("═══════════════════════════════════════════════════════")
    }

    # Step 2.1: Enrich Age
    if (verbose) message("\n[2.1] Enriching Age...")
    data <- enrich_age(data)

    # Step 2.2: Enrich Length of Stay
    if (verbose) message("\n[2.2] Enriching Length of Stay...")
    data <- enrich_los(data)

    # Step 2.3: Enrich infection type (CAI/HAI)
    if (verbose) message("\n[2.3] Enriching infection type...")
    data <- enrich_infection_type(data, hai_cutoff = config$hai_cutoff)

    # Step 2.4: Enrich hospital department (if enabled)
    if (config$infer_department) {
      if (verbose) message("\n[2.4] Enriching hospital department...")
      data <- enrich_hospital_department(data)
    }

    # Step 2.5: Groom optional columns
    if (verbose) message("\n[2.5] Grooming optional columns...")
    data <- groom_optional_columns(data)
  }

  # ========================================
  # PHASE 3: DERIVATION
  # ========================================
  if ("derive" %in% phases) {
    if (verbose) {
      message("\n═══════════════════════════════════════════════════════")
      message("  PHASE 3: DERIVATION")
      message("═══════════════════════════════════════════════════════")
    }

    # Step 3.1: Assign age bins
    if ("Age" %in% names(data)) {
      if (verbose) message("\n[3.1] Assigning age bins...")
      data <- assign_age_bins(data, bins = config$age_bins)
    }

    # Step 3.2: Calculate Length of Stay (if not already done)
    if (all(c("date_of_admission", "date_of_final_outcome") %in% names(data))) {
      if (verbose) message("\n[3.2] Calculating Length of Stay...")
      data <- calculate_los(data)
    }

    # Step 3.3: Extract organism taxonomy
    if ("organism_normalized" %in% names(data)) {
      if (verbose) message("\n[3.3] Extracting organism taxonomy...")
      data <- extract_genus(data)
      data <- extract_species(data)
      data <- classify_org_group(data)
    }

    # Step 3.4: Classify antibiotic classes
    if ("antibiotic_normalized" %in% names(data)) {
      if (verbose) message("\n[3.4] Classifying antibiotic classes...")
      data <- classify_antibiotic_class(data)
      data <- classify_aware(data)
    }

    # Step 3.5: Create event IDs (patient → isolate level)
    if (all(c("patient_id", "date_of_culture", "organism_normalized") %in% names(data))) {
      if (verbose) message("\n[3.5] Creating event IDs...")
      data <- create_event_ids(data, gap_days = config$event_gap_days)
    }

    # Step 3.6: Flag contaminants
    if (verbose) message("\n[3.6] Flagging contaminants...")
    data <- flag_contaminants(
      data,
      method = config$contaminant_method
    )

    # Step 3.7: Classify mortality
    if ("final_outcome" %in% names(data)) {
      if (verbose) message("\n[3.7] Classifying mortality...")
      data <- classify_mortality(
        data,
        window = config$mortality_window
      )
    }

    # Step 3.8: Polymicrobial identification and weighting
    poly_required_cols <- c("patient_id", "specimen_type", "date_of_culture", "organism_normalized")
    if (all(poly_required_cols %in% names(data))) {
      if (verbose) message("\n[3.8] Identifying polymicrobial infections...")
      data <- flag_polymicrobial(data)

      if (verbose) message("\n[3.9] Computing polymicrobial weights...")
      data <- compute_polymicrobial_weight(data, method = "monomicrobial_proportion")
    }

    # Step 3.9: Apply death weights (placeholder — not yet implemented)
    if ("mortality_infection" %in% names(data)) {
      if (verbose) message("\n[3.10] Death weights: skipped (not yet implemented)")
    }

    # Step 3.10: MDR/XDR classification
    if (all(c("organism_normalized", "antibiotic_class", "antibiotic_value") %in% names(data))) {
      if (verbose) message("\n[3.11] Classifying MDR/XDR...")
      data <- classify_mdr(data, definition = config$mdr_definition)
      data <- classify_xdr(data, definition = config$xdr_definition)
    }

    # Step 3.11: Map to RR pathogens and classes (if enabled)
    if ("organism_normalized" %in% names(data)) {
      if (verbose) message("\n[3.12] Mapping to RR categories...")
      data <- map_to_rr_pathogen(data)
    }
    if ("antibiotic_class" %in% names(data)) {
      data <- map_class_to_rr(data)
    }
  }

  # ========================================
  # POST-VALIDATION
  # ========================================
  if (validate && verbose) {
    message("\n─── Post-Processing Validation ───")
  }

  if (validate) {
    # Check logical consistency
    consistency <- check_logical_consistency(data, checks = "all", stop_on_failure = FALSE)
    log$validation$post_consistency <- consistency

    # Final quality check
    final_quality <- validate_data_quality(
      data,
      min_rows = 1,
      max_missing_pct = 95,
      required_cols = character(0),
      stop_on_failure = FALSE
    )
    log$validation$post_quality <- final_quality
  }

  # ========================================
  # GENERATE REPORT
  # ========================================
  report <- NULL
  if (generate_report) {
    if (verbose) message("\n─── Generating Preprocessing Report ───")
    report <- generate_preprocessing_report(
      raw_data = original_summary,
      processed_data = data,
      preprocessing_log = log,
      config = config
    )
  }

  # ========================================
  # FINALIZE
  # ========================================
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  metadata <- list(
    pipeline_version = "0.1.0",
    start_time = start_time,
    end_time = end_time,
    elapsed_seconds = elapsed,
    phases_run = phases,
    n_rows_input = nrow(original_summary$data),
    n_rows_output = nrow(data),
    n_cols_input = ncol(original_summary$data),
    n_cols_output = ncol(data)
  )

  # Create result object
  result <- list(
    data = data,
    config = config,
    log = log,
    report = report,
    metadata = metadata
  )
  class(result) <- c("amr_result", "list")

  if (verbose) {
    message("\n═══════════════════════════════════════════════════════")
    message("  PIPELINE COMPLETE")
    message("═══════════════════════════════════════════════════════")
    message(sprintf("Completed: %s", format(end_time, "%Y-%m-%d %H:%M:%S")))
    message(sprintf("Elapsed time: %.2f seconds", elapsed))
    message(sprintf("Output: %d rows × %d columns", nrow(data), ncol(data)))
    message("═══════════════════════════════════════════════════════")
  }

  return(result)
}


#' Summary Method for AMR Preprocessing Results
#'
#' Prints a formatted summary of preprocessing pipeline results.
#'
#' @param object amr_result object from amr_preprocess()
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the object
#' @export
#'
#' @examples
#' \dontrun{
#' result <- amr_preprocess(data)
#' summary(result)
#' }
summary.amr_result <- function(object, ...) {
  cat("═══════════════════════════════════════════════════════\n")
  cat("  anumaan Preprocessing Summary\n")
  cat("═══════════════════════════════════════════════════════\n\n")

  # Metadata
  cat("Pipeline Information:\n")
  cat(sprintf("  Version: %s\n", object$metadata$pipeline_version))
  cat(sprintf("  Phases run: %s\n", paste(object$metadata$phases_run, collapse = ", ")))
  cat(sprintf("  Execution time: %.2f seconds\n", object$metadata$elapsed_seconds))
  cat("\n")

  # Data dimensions
  cat("Data Transformation:\n")
  cat(sprintf(
    "  Input:  %d rows × %d columns\n",
    object$metadata$n_rows_input,
    object$metadata$n_cols_input
  ))
  cat(sprintf(
    "  Output: %d rows × %d columns\n",
    object$metadata$n_rows_output,
    object$metadata$n_cols_output
  ))
  cat(sprintf(
    "  Rows retained: %.1f%%\n",
    100 * object$metadata$n_rows_output / object$metadata$n_rows_input
  ))
  cat("\n")

  # Column mapping
  if (!is.null(object$log$column_mapping)) {
    n_mapped <- sum(!is.na(object$log$column_mapping$standard_name))
    cat(sprintf("Column Mapping: %d columns mapped\n", n_mapped))
    cat("\n")
  }

  # Validation
  if (!is.null(object$log$validation)) {
    cat("Validation:\n")

    if (!is.null(object$log$validation$duplicates_removed)) {
      cat(sprintf("  Duplicates removed: %d\n", object$log$validation$duplicates_removed))
    }

    if (!is.null(object$log$validation$post_consistency)) {
      if (object$log$validation$post_consistency$consistent) {
        cat("  Logical consistency: PASSED ✓\n")
      } else {
        cat(sprintf(
          "  Logical consistency: FAILED (%d issues)\n",
          object$log$validation$post_consistency$n_issues
        ))
      }
    }

    if (!is.null(object$log$validation$post_quality)) {
      cat(sprintf(
        "  Data completeness: %.1f%%\n",
        object$log$validation$post_quality$overall_completeness * 100
      ))
    }
    cat("\n")
  }

  # Key derived variables
  data <- object$data
  cat("Derived Variables:\n")

  if ("event_id" %in% names(data)) {
    cat(sprintf("  Events created: %d\n", dplyr::n_distinct(data$event_id)))
  }
  if ("is_polymicrobial" %in% names(data)) {
    # Count polymicrobial episodes (unique patient-specimen-date combinations)
    poly_data <- data[data$is_polymicrobial == 1, ]
    if (nrow(poly_data) > 0 && all(c("patient_id", "specimen_type", "date_of_culture") %in% names(poly_data))) {
      n_poly <- nrow(unique(poly_data[, c("patient_id", "specimen_type", "date_of_culture")]))
    } else {
      n_poly <- sum(data$is_polymicrobial == 1 & !duplicated(data$event_id), na.rm = TRUE)
    }
    cat(sprintf("  Polymicrobial episodes: %d\n", n_poly))
  }
  if ("mdr" %in% names(data)) {
    n_mdr <- sum(data$mdr == "MDR" & !duplicated(data$event_id), na.rm = TRUE)
    cat(sprintf("  MDR infections: %d\n", n_mdr))
  }
  if ("mortality_infection" %in% names(data)) {
    n_deaths <- sum(data$mortality_infection == "Yes" & !duplicated(data$event_id), na.rm = TRUE)
    cat(sprintf("  Infection-related deaths: %d\n", n_deaths))
  }

  cat("\n═══════════════════════════════════════════════════════\n")
  cat("Access components:\n")
  cat("  result$data    - Preprocessed data frame\n")
  cat("  result$config  - Configuration used\n")
  cat("  result$log     - Processing logs\n")
  cat("  result$report  - Preprocessing report\n")
  cat("═══════════════════════════════════════════════════════\n")

  invisible(object)
}


#' Print Method for AMR Preprocessing Results
#'
#' @param x amr_result object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the object
#' @export
print.amr_result <- function(x, ...) {
  cat("anumaan preprocessing result\n")
  cat(sprintf("  %d rows × %d columns\n", nrow(x$data), ncol(x$data)))
  cat(sprintf("  Phases: %s\n", paste(x$metadata$phases_run, collapse = ", ")))
  cat("\nUse summary() for detailed information\n")
  invisible(x)
}
