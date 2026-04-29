# prep_schema.R
# Layer 2: Column mapping + schema alignment
# Purpose: map raw dataset columns into the minimum standard AMR schema.
#
# Canonical staged API:
#   prep_build_column_map()  -> build/validate the rename map
#   prep_apply_column_map()  -> apply the rename
#   prep_assert_standard_names() -> confirm required columns exist
#
# Convenience wrapper (calls staged API internally):
#   prep_standardize_column_names()  -> one-step rename with optional fuzzy match
#
# Other functions:
#   default_column_mappings   - generic alias list
#   detect_preprocessing_capabilities
#   prep_report_capabilities
#   prep_clean_optional_columns  (moved from prep_clean_and_standardize.R)
#
# Dataset-specific column maps belong in the analysis layer, NOT this package.


# ---------------------------------------------------------------------------
# Default column mappings (generic fuzzy-match alias list)
# ---------------------------------------------------------------------------

#' Default column name mappings for fuzzy matching
#'
#' @format A named list with elements mapping standard names to character
#'   vectors of acceptable aliases.
#' @export
default_column_mappings <- list(
  patient_id = c(
    "PatientInformation_id", "PatientInformation", "patient_ID",
    "Patient_ID", "PatientID", "Subject_ID", "MRN",
    "medical_record_number", "UHID", "patient_no", "Case_ID"
  ),
  gender = c(
    "Gender", "sex", "Sex", "patient_gender", "Patient_Gender",
    "gender_code", "sex_code"
  ),
  state = c(
    "State", "patient_state", "state_name", "region",
    "State_Name", "state_code"
  ),
  location = c(
    "Location", "city", "hospital_city", "hospital_location",
    "site", "City", "Hospital_Location"
  ),
  DOB = c(
    "date_of_birth", "DateOfBirth", "birth_date", "dob", "DOB",
    "Date_of_Birth", "BirthDate"
  ),
  Age = c("age", "patient_age", "Age_Years", "age_years", "AGE"),
  date_of_admission = c(
    "admission_date", "Date.of.admission",
    "Date_of_admission_in_hospital", "AdmissionDate",
    "hospital_admission_date", "admit_date",
    "Date.of.admission.in.hospital"
  ),
  date_of_culture = c(
    "date_of_event", "Date.of.event", "culture_date",
    "collection_date", "sample_date", "specimen_date",
    "event_date", "culture_collection_date",
    "Date_of_culture", "CultureDate"
  ),
  date_of_final_outcome = c(
    "Date.of.14.day.outcome", "outcome_date",
    "discharge_date", "death_date", "Date_of_outcome",
    "final_outcome_date", "DateOfOutcome",
    "date_of_discharge"
  ),
  final_outcome = c(
    "Final.outcome", "outcome", "patient_outcome",
    "status", "final_status", "discharge_status",
    "Outcome", "Status"
  ),
  organism_name = c(
    "Organism", "organism", "pathogen", "pathogen_name",
    "bacteria", "microorganism", "organism_identified",
    "OrganismName", "isolated_organism"
  ),
  antibiotic_name = c(
    "antibiotic", "drug", "drug_name", "antimicrobial",
    "antimicrobial_name", "antibiotic_tested",
    "drug_tested", "Antibiotic", "AntibioticName"
  ),
  antibiotic_value = c(
    "antibiotic_result", "susceptibility", "resistance",
    "result", "remarks", "antibiotic_remarks",
    "susceptibility_result", "Remarks", "Result",
    "susceptibility_status"
  ),
  specimen_type = c(
    "specimen", "sample_type", "sample", "Sample_type1_name",
    "source", "specimen_source", "sample_source",
    "culture_source", "SpecimenType", "SampleType"
  ),
  diagnosis = c(
    "Diagnosis", "diagnosis_1", "primary_diagnosis",
    "clinical_diagnosis", "Diagnosis_1", "ICD_code",
    "diagnosis_code"
  )
)


# ---------------------------------------------------------------------------
# Moved function: prep_standardize_column_names
# ---------------------------------------------------------------------------

#' Standardize Column Names to Package Convention
#'
#' Convenience wrapper around \code{prep_build_column_map()} +
#' \code{prep_apply_column_map()}. Builds an exact-match map from
#' \code{mapping} aliases, optionally extends it with fuzzy matching, then
#' applies the rename in one call.
#'
#' For more control (e.g. dataset-specific maps, post-rename assertions) use
#' the staged API directly:
#' \code{prep_build_column_map()} -> \code{prep_apply_column_map()} ->
#' \code{prep_assert_standard_names()}.
#'
#' @param data A data frame with raw column names.
#' @param mapping Named list of standard_name -> aliases vectors. Default uses
#'   \code{default_column_mappings}.
#' @param fuzzy_match Logical. Enable fuzzy matching for unresolved standards.
#'   Default TRUE.
#' @param fuzzy_threshold Numeric. Max Jaro-Winkler distance (0-1). Default 0.3.
#' @param interactive Logical. Prompt user to confirm fuzzy matches. Default FALSE.
#'
#' @return List: \code{data} (renamed), \code{mapping_log}, \code{unmapped}.
#' @export
prep_standardize_column_names <- function(data,
                                          mapping         = default_column_mappings,
                                          fuzzy_match     = TRUE,
                                          fuzzy_threshold = 0.3,
                                          interactive     = FALSE) {
  # Build exact-match map via canonical helper
  column_map <- vapply(names(mapping), function(std_name) {
    hit <- intersect(mapping[[std_name]], names(data))
    if (length(hit) > 0L) hit[1L] else NA_character_
  }, character(1))
  column_map <- column_map[!is.na(column_map)]

  mapping_log <- setNames(
    lapply(names(column_map), function(s) list(original = column_map[[s]], method = "exact_match")),
    names(column_map)
  )

  # Optional fuzzy extension
  if (fuzzy_match) {
    already_mapped <- names(data) %in% unname(column_map)
    unmapped_cols  <- names(data)[!already_mapped]
    unresolved     <- setdiff(names(mapping), names(column_map))

    for (std_name in unresolved) {
      if (length(unmapped_cols) == 0L) break
      distances <- stringdist::stringdist(tolower(std_name), tolower(unmapped_cols), method = "jw")
      best_idx  <- which.min(distances)
      if (length(best_idx) > 0 && distances[best_idx] < fuzzy_threshold) {
        accept <- if (interactive) {
          message(sprintf("Fuzzy: '%s' -> '%s' (dist %.2f). Accept? (y/n)",
                          unmapped_cols[best_idx], std_name, distances[best_idx]))
          tolower(readline()) == "y"
        } else {
          message(sprintf("Auto fuzzy: '%s' -> '%s' (dist %.2f)",
                          unmapped_cols[best_idx], std_name, distances[best_idx]))
          TRUE
        }
        if (accept) {
          column_map[std_name]    <- unmapped_cols[best_idx]
          mapping_log[[std_name]] <- list(original = unmapped_cols[best_idx],
                                          method   = "fuzzy_match",
                                          distance = distances[best_idx])
          unmapped_cols <- unmapped_cols[-best_idx]
        }
      }
    }
  }

  # Apply via canonical helper
  data          <- prep_apply_column_map(data, column_map)
  unmapped_final <- setdiff(names(data), names(mapping))

  list(data = data, mapping_log = mapping_log, unmapped = unmapped_final)
}


# ---------------------------------------------------------------------------
# New functions (Layer 2)
# ---------------------------------------------------------------------------

#' Build and Validate a Column Map Against a Dataset
#'
#' Validates a user-supplied column map against the actual columns present in a
#' data frame. Reports which mappings are present, which are absent, and which
#' data columns are not covered.
#'
#' Dataset-specific maps (ICMR stewardship, surveillance, AIIMS, etc.) should
#' be defined in the analysis layer and passed in via \code{column_map}.
#' If no map is supplied, falls back to auto-detection using
#' \code{default_column_mappings} aliases.
#'
#' @param data Data frame with raw column names.
#' @param column_map Named character vector: standard_name = raw_name.
#'   Pass a dataset-specific map from your analysis layer. If NULL,
#'   generic alias detection is attempted.
#' @param custom_map Named character vector. Additional raw -> standard overrides
#'   applied on top of \code{column_map}. Optional.
#'
#' @return Named character vector with only mappings whose raw names exist in data.
#' @export
prep_build_column_map <- function(data,
                                  column_map = NULL,
                                  custom_map = NULL) {
  if (is.null(column_map)) {
    base_map <- vapply(names(default_column_mappings), function(std_name) {
      aliases <- default_column_mappings[[std_name]]
      hit     <- intersect(aliases, names(data))
      if (length(hit) > 0L) hit[1L] else NA_character_
    }, character(1))
    base_map <- base_map[!is.na(base_map)]
  } else {
    if (!is.character(column_map) || is.null(names(column_map)))
      stop("`column_map` must be a named character vector (standard_name = raw_name).")
    base_map <- column_map
  }

  if (!is.null(custom_map)) {
    if (!is.character(custom_map) || is.null(names(custom_map)))
      stop("`custom_map` must be a named character vector.")
    base_map[names(custom_map)] <- custom_map
  }

  present <- base_map %in% names(data)
  if (any(!present))
    message(sprintf("[prep_build_column_map] %d raw column(s) not found in data: %s",
                    sum(!present), paste(base_map[!present], collapse = ", ")))

  base_map[present]
}


#' Apply Column Map to Rename Columns
#'
#' Renames data frame columns using a standard_name = raw_name map.
#' Warns on unmapped columns; does not drop them.
#'
#' @param data Data frame.
#' @param column_map Named character vector from \code{prep_build_column_map()}.
#'
#' @return Data frame with standard column names where mapping existed.
#' @export
prep_apply_column_map <- function(data, column_map) {
  if (is.null(column_map) || length(column_map) == 0L)
    return(data)

  # Invert: raw_name -> standard_name
  inverted <- stats::setNames(names(column_map), column_map)

  # Only rename columns that exist in data
  cols_to_rename <- intersect(names(data), names(inverted))

  if (length(cols_to_rename) > 0L) {
    names(data)[match(cols_to_rename, names(data))] <- inverted[cols_to_rename]
    message(sprintf("[prep_apply_column_map] Renamed %d column(s).", length(cols_to_rename)))
  }

  unmapped <- setdiff(names(data), c(names(column_map), cols_to_rename))
  if (length(unmapped) > 0L)
    message(sprintf("[prep_apply_column_map] %d column(s) not in map (kept as-is): %s",
                    length(unmapped),
                    paste(utils::head(unmapped, 10), collapse = ", ")))

  data
}


#' Assert Standard Names Are Present
#'
#' Confirms that a column mapping succeeded for critical standard columns.
#' Stops on failure when strict = TRUE.
#'
#' @param data Data frame after column mapping.
#' @param required_standard_names Character vector of standard column names that
#'   must be present.
#' @param strict Logical. Stop on failure. Default TRUE.
#'
#' @return Invisibly returns data. Stops or warns on missing columns.
#' @export
prep_assert_standard_names <- function(data,
                                       required_standard_names,
                                       strict = TRUE) {
  missing <- setdiff(required_standard_names, names(data))

  if (length(missing) > 0L) {
    msg <- sprintf(
      "[prep_assert_standard_names] %d required standard column(s) missing after mapping: %s",
      length(missing), paste(missing, collapse = ", ")
    )
    if (strict) stop(msg) else warning(msg)
  } else {
    message(sprintf(
      "[prep_assert_standard_names] All %d required standard columns present.",
      length(required_standard_names)
    ))
  }

  invisible(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 2)
# ---------------------------------------------------------------------------

#' Detect Preprocessing Capabilities
#'
#' Inspects a data frame (after column mapping) and returns a named logical
#' vector indicating which preprocessing steps can be executed given the
#' columns that are present.
#'
#' This is used by \code{run_preprocess()} to decide which layers to run and
#' which to skip gracefully.
#'
#' @param data Data frame after column mapping.
#'
#' @return Named logical vector. TRUE means the step can run; FALSE means
#'   required inputs are absent and the step will be skipped with a warning.
#' @export
detect_preprocessing_capabilities <- function(data) {
  cols <- names(data)

  c(
    parse_dates          = any(c("admission_date", "culture_date", "outcome_date") %in% cols),
    derive_hai           = all(c("admission_date", "culture_date") %in% cols),
    derive_los           = all(c("admission_date", "outcome_date") %in% cols),
    derive_age           = all(c("dob", "culture_date") %in% cols) ||
                            "age" %in% cols,
    standardize_organism = "organism_name" %in% cols,
    standardize_antibiotic = "antibiotic_name" %in% cols,
    standardize_specimen = "sample_type" %in% cols,
    harmonize_ast        = "ast_value_harmonized" %in% cols ||
                            "antibiotic_value" %in% cols,
    create_events        = all(c("patient_id", "culture_date", "organism_name") %in% cols),
    flag_contaminants    = all(c("organism_name", "sample_type") %in% cols),
    flag_polymicrobial   = all(c("patient_id", "organism_name") %in% cols),
    classify_mdr         = "ast_value_harmonized" %in% cols &&
                            "antibiotic_name" %in% cols,
    build_outputs        = "patient_id" %in% cols
  )
}


#' Report Preprocessing Capabilities
#'
#' Prints a human-readable summary of which preprocessing steps can be run
#' given the current data columns.
#'
#' @param data Data frame after column mapping. Or pass the result of
#'   \code{detect_preprocessing_capabilities()} directly.
#'
#' @return Invisibly returns the capabilities vector.
#' @export
prep_report_capabilities <- function(data) {
  caps <- if (is.logical(data)) data else detect_preprocessing_capabilities(data)

  enabled  <- names(caps)[caps]
  disabled <- names(caps)[!caps]

  message("=== Preprocessing Capabilities ===")
  if (length(enabled) > 0) {
    message("Enabled (", length(enabled), "):")
    for (s in enabled) message("  [+] ", s)
  }
  if (length(disabled) > 0) {
    message("Skipped / missing inputs (", length(disabled), "):")
    for (s in disabled) message("  [-] ", s)
  }
  message("===================================")

  invisible(caps)
}


# ---------------------------------------------------------------------------
# Moved function: prep_clean_optional_columns (from prep_clean_and_standardize.R)
# ---------------------------------------------------------------------------

#' Clean Optional Columns
#'
#' Cleans and standardizes all optional columns in a single pass.
#' Handles value normalization, whitespace trimming, and consistency checks.
#'
#' @param data Data frame.
#' @param optional_cols Character vector. Optional column names to groom.
#'   If NULL, grooms all known optional columns. Default NULL.
#'
#' @return Data frame with groomed optional columns.
#' @export
prep_clean_optional_columns <- function(data, optional_cols = NULL) {
  known_optional <- c(
    "infection_type", "device_inserted", "hospital_department",
    "aware_category", "unit_type", "comorbidities", "hospital_location",
    "previous_history", "specimen_category"
  )

  if (is.null(optional_cols)) {
    optional_cols <- intersect(known_optional, names(data))
  } else {
    optional_cols <- intersect(optional_cols, names(data))
  }

  if (length(optional_cols) == 0) {
    message("No optional columns to groom")
    return(data)
  }

  message(sprintf("Grooming %d optional columns: %s",
                  length(optional_cols),
                  paste(optional_cols[1:min(3, length(optional_cols))], collapse = ", ")))

  for (col in optional_cols) {
    if (is.character(data[[col]])) {
      data[[col]] <- trimws(data[[col]])
      data[[col]][data[[col]] == ""] <- NA_character_

      if (!col %in% c("comorbidities", "previous_history"))
        data[[col]] <- tools::toTitleCase(tolower(data[[col]]))

      if (col == "infection_type") {
        data[[col]] <- dplyr::case_when(
          grepl("^comm|^cai",         tolower(data[[col]])) ~ "CAI",
          grepl("^hosp|^hai|^nosoco", tolower(data[[col]])) ~ "HAI",
          TRUE ~ data[[col]]
        )
      }

      if (col == "unit_type") {
        data[[col]] <- dplyr::case_when(
          grepl("icu|intensive|critical", tolower(data[[col]])) ~ "ICU",
          grepl("ward|general",           tolower(data[[col]])) ~ "Ward",
          grepl("er|emergency",           tolower(data[[col]])) ~ "Emergency",
          TRUE ~ data[[col]]
        )
      }
    }
  }

  completeness <- sapply(optional_cols, function(col)
    sum(!is.na(data[[col]])) / nrow(data) * 100)

  message("\nOptional column completeness:")
  print(data.frame(column          = names(completeness),
                   completeness_pct = as.numeric(completeness),
                   stringsAsFactors = FALSE) %>%
          dplyr::arrange(dplyr::desc(completeness_pct)))

  return(data)
}
