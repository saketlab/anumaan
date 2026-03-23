# enrich.R
# Phase 2: Optional variable enrichment and inference

#' Enrich Age
#'
#' Derives Age when missing using DOB and culture date. Uses multi-path
#' logic with confidence scoring.
#'
#' @param data Data frame
#' @param age_col Character. Age column. Default "Age".
#' @param dob_col Character. Date of birth column. Default "DOB".
#' @param date_col Character. Reference date for age calculation.
#'   Default "date_of_culture".
#' @param overwrite Logical. If TRUE, recalculates age even if present.
#'   Default FALSE (only fill missing).
#'
#' @return Data frame with Age column enriched
#' @export
#'
#' @examples
#' \dontrun{
#' data_enriched <- enrich_age(data)
#' }
enrich_age <- function(data,
                       age_col = "Age",
                       dob_col = "DOB",
                       date_col = "date_of_culture",
                       overwrite = FALSE) {
  # Check if age column exists
  has_age <- age_col %in% names(data)
  has_dob <- dob_col %in% names(data)
  has_date <- date_col %in% names(data)

  if (!has_age) {
    data[[age_col]] <- NA_real_
  }

  # Initialize tracking columns if they don't exist
  if (!"age_method" %in% names(data)) {
    data$age_method <- NA_character_
    data$age_confidence <- NA_character_
  }

  # Path A: Calculate from DOB (HIGH CONFIDENCE)
  if (has_dob && has_date) {
    n_before_missing <- sum(is.na(data[[age_col]]))

    data <- data %>%
      dplyr::mutate(
        calculated_age = dplyr::case_when(
          !is.na(!!rlang::sym(dob_col)) & !is.na(!!rlang::sym(date_col)) ~
            as.numeric(difftime(!!rlang::sym(date_col), !!rlang::sym(dob_col), units = "days")) / 365.25,
          TRUE ~ NA_real_
        ),
        !!age_col := dplyr::case_when(
          # Overwrite mode: use calculated
          overwrite & !is.na(calculated_age) ~ calculated_age,
          # Fill missing mode: only fill if missing
          is.na(!!rlang::sym(age_col)) & !is.na(calculated_age) ~ calculated_age,
          # Keep existing
          TRUE ~ !!rlang::sym(age_col)
        ),
        age_method = dplyr::case_when(
          !is.na(calculated_age) & (overwrite | is.na(age_method)) ~ "calculated_from_dob",
          !is.na(age_method) ~ age_method,
          !is.na(!!rlang::sym(age_col)) ~ "provided",
          TRUE ~ NA_character_
        ),
        age_confidence = dplyr::case_when(
          age_method == "calculated_from_dob" ~ "high",
          age_method == "provided" ~ "high",
          TRUE ~ age_confidence
        )
      ) %>%
      dplyr::select(-calculated_age)

    n_after_missing <- sum(is.na(data[[age_col]]))
    n_enriched <- n_before_missing - n_after_missing

    if (n_enriched > 0) {
      message(sprintf(
        "Enriched Age: %d rows filled using DOB calculation",
        n_enriched
      ))
    }
  } else {
    # Mark provided ages
    data <- data %>%
      dplyr::mutate(
        age_method = dplyr::case_when(
          !is.na(!!rlang::sym(age_col)) & is.na(age_method) ~ "provided",
          TRUE ~ age_method
        ),
        age_confidence = dplyr::case_when(
          age_method == "provided" ~ "high",
          TRUE ~ age_confidence
        )
      )
  }

  # Summary
  age_summary <- data %>%
    dplyr::count(age_method, age_confidence) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nAge enrichment summary:")
  print(age_summary)

  n_still_missing <- sum(is.na(data[[age_col]]))
  if (n_still_missing > 0) {
    message(sprintf(
      "⚠ Warning: %d rows still missing Age (%.1f%%)",
      n_still_missing,
      100 * n_still_missing / nrow(data)
    ))
  }

  return(data)
}


#' Enrich Length of Stay
#'
#' Derives Length of Stay (LOS) from admission and outcome dates.
#'
#' @param data Data frame
#' @param los_col Character. LOS column name. Default "Length_of_stay".
#' @param admission_col Character. Admission date. Default "date_of_admission".
#' @param outcome_col Character. Outcome date. Default "date_of_final_outcome".
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with LOS enriched
#' @export
#'
#' @examples
#' \dontrun{
#' data_enriched <- enrich_los(data)
#' }
enrich_los <- function(data,
                       los_col = "Length_of_stay",
                       admission_col = "date_of_admission",
                       outcome_col = "date_of_final_outcome",
                       overwrite = FALSE) {
  # Check columns
  has_los <- los_col %in% names(data)
  has_admission <- admission_col %in% names(data)
  has_outcome <- outcome_col %in% names(data)

  if (!has_los) {
    data[[los_col]] <- NA_real_
  }

  if (!has_admission || !has_outcome) {
    message(sprintf(
      "⚠ Cannot calculate LOS: missing '%s' or '%s'",
      admission_col, outcome_col
    ))
    return(data)
  }

  n_before_missing <- sum(is.na(data[[los_col]]))

  # Calculate LOS
  data <- data %>%
    dplyr::mutate(
      calculated_los = dplyr::case_when(
        !is.na(!!rlang::sym(admission_col)) & !is.na(!!rlang::sym(outcome_col)) ~
          as.numeric(difftime(!!rlang::sym(outcome_col), !!rlang::sym(admission_col), units = "days")),
        TRUE ~ NA_real_
      ),
      !!los_col := dplyr::case_when(
        overwrite & !is.na(calculated_los) ~ calculated_los,
        is.na(!!rlang::sym(los_col)) & !is.na(calculated_los) ~ calculated_los,
        TRUE ~ !!rlang::sym(los_col)
      ),
      los_method = dplyr::case_when(
        !is.na(calculated_los) ~ "calculated_from_dates",
        !is.na(!!rlang::sym(los_col)) ~ "provided",
        TRUE ~ NA_character_
      ),
      los_confidence = dplyr::case_when(
        los_method %in% c("calculated_from_dates", "provided") ~ "high",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-calculated_los)

  # Flag negative LOS (data issue)
  negative_los <- data %>%
    dplyr::filter(!is.na(!!rlang::sym(los_col)) & !!rlang::sym(los_col) < 0)

  if (nrow(negative_los) > 0) {
    warning(sprintf(
      "⚠ Warning: %d rows have negative LOS (date inconsistency)",
      nrow(negative_los)
    ))
  }

  n_after_missing <- sum(is.na(data[[los_col]]))
  n_enriched <- n_before_missing - n_after_missing

  if (n_enriched > 0) {
    message(sprintf(
      "Enriched LOS: %d rows filled",
      n_enriched
    ))
  }

  # Summary statistics
  los_stats <- data %>%
    dplyr::filter(!is.na(!!rlang::sym(los_col))) %>%
    dplyr::summarise(
      mean_los = mean(!!rlang::sym(los_col), na.rm = TRUE),
      median_los = median(!!rlang::sym(los_col), na.rm = TRUE),
      min_los = min(!!rlang::sym(los_col), na.rm = TRUE),
      max_los = max(!!rlang::sym(los_col), na.rm = TRUE)
    )

  message("\nLOS summary:")
  print(los_stats)

  return(data)
}


#' Enrich Infection Type
#'
#' Infers Community-Acquired (CAI) vs Hospital-Acquired (HAI) infection
#' using admission-culture date gap.
#'
#' @param data Data frame
#' @param infection_type_col Character. Infection type column.
#'   Default "infection_type".
#' @param admission_col Character. Admission date. Default "date_of_admission".
#' @param culture_col Character. Culture date. Default "date_of_culture".
#' @param hai_cutoff Numeric. Days after admission to classify as HAI.
#'   Default 2 (48 hours).
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with infection_type enriched
#' @export
#'
#' @examples
#' \dontrun{
#' # Default 2-day cutoff
#' data_enriched <- enrich_infection_type(data)
#'
#' # 3-day cutoff
#' data_enriched <- enrich_infection_type(data, hai_cutoff = 3)
#' }
enrich_infection_type <- function(data,
                                  infection_type_col = "infection_type",
                                  admission_col = "date_of_admission",
                                  culture_col = "date_of_culture",
                                  hai_cutoff = 2,
                                  overwrite = FALSE) {
  # Check columns
  has_infection_type <- infection_type_col %in% names(data)
  has_admission <- admission_col %in% names(data)
  has_culture <- culture_col %in% names(data)

  if (!has_infection_type) {
    data[[infection_type_col]] <- NA_character_
  }

  if (!has_admission || !has_culture) {
    message(sprintf(
      "⚠ Cannot infer infection type: missing '%s' or '%s'",
      admission_col, culture_col
    ))
    return(data)
  }

  n_before_missing <- sum(is.na(data[[infection_type_col]]))

  message(sprintf(
    "Inferring infection type using %d-day HAI cutoff...",
    hai_cutoff
  ))

  # Infer infection type
  data <- data %>%
    dplyr::mutate(
      days_to_culture = dplyr::case_when(
        !is.na(!!rlang::sym(admission_col)) & !is.na(!!rlang::sym(culture_col)) ~
          as.numeric(difftime(!!rlang::sym(culture_col), !!rlang::sym(admission_col), units = "days")),
        TRUE ~ NA_real_
      ),
      inferred_type = dplyr::case_when(
        !is.na(days_to_culture) & days_to_culture >= hai_cutoff ~ "HAI",
        !is.na(days_to_culture) & days_to_culture < hai_cutoff ~ "CAI",
        TRUE ~ NA_character_
      ),
      !!infection_type_col := dplyr::case_when(
        overwrite & !is.na(inferred_type) ~ inferred_type,
        is.na(!!rlang::sym(infection_type_col)) & !is.na(inferred_type) ~ inferred_type,
        TRUE ~ !!rlang::sym(infection_type_col)
      ),
      infection_type_method = dplyr::case_when(
        !is.na(inferred_type) ~ sprintf("inferred_%dday_cutoff", hai_cutoff),
        !is.na(!!rlang::sym(infection_type_col)) ~ "provided",
        TRUE ~ NA_character_
      ),
      infection_type_confidence = dplyr::case_when(
        infection_type_method == "provided" ~ "high",
        !is.na(inferred_type) ~ "medium",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-inferred_type, -days_to_culture)

  n_after_missing <- sum(is.na(data[[infection_type_col]]))
  n_enriched <- n_before_missing - n_after_missing

  if (n_enriched > 0) {
    message(sprintf(
      "Enriched infection_type: %d rows filled",
      n_enriched
    ))
  }

  # Summary
  type_summary <- data %>%
    dplyr::count(!!rlang::sym(infection_type_col), infection_type_method) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nInfection type distribution:")
  print(type_summary)

  return(data)
}


#' Enrich Hospital Department
#'
#' Infers hospital department from contextual data (specimen type,
#' diagnosis, patient demographics).
#'
#' @param data Data frame
#' @param department_col Character. Department column. Default "hospital_department".
#' @param specimen_col Character. Specimen type column. Default "specimen_type".
#' @param diagnosis_col Character. Diagnosis column. Default "diagnosis_1".
#' @param age_col Character. Age column. Default "Age".
#' @param overwrite Logical. Recalculate even if present. Default FALSE.
#'
#' @return Data frame with hospital_department enriched
#' @export
#'
#' @examples
#' \dontrun{
#' data_enriched <- enrich_hospital_department(data)
#' }
enrich_hospital_department <- function(data,
                                       department_col = "hospital_department",
                                       specimen_col = "specimen_type",
                                       diagnosis_col = "diagnosis_1",
                                       age_col = "Age",
                                       overwrite = FALSE) {
  # Check columns
  has_department <- department_col %in% names(data)
  has_specimen <- specimen_col %in% names(data)
  has_diagnosis <- diagnosis_col %in% names(data)
  has_age <- age_col %in% names(data)

  if (!has_department) {
    data[[department_col]] <- NA_character_
  }

  # Need at least one contextual variable
  if (!has_specimen && !has_diagnosis && !has_age) {
    message("⚠ Cannot infer department: no contextual data available")
    return(data)
  }

  n_before_missing <- sum(is.na(data[[department_col]]))

  message("Inferring hospital department from contextual data...")

  # Create dummy columns if they don't exist (for safe case_when evaluation)
  if (!has_age) data[[age_col]] <- NA_real_
  if (!has_specimen) data[[specimen_col]] <- NA_character_
  if (!has_diagnosis) data[[diagnosis_col]] <- NA_character_

  # Inference rules (heuristic-based)
  data <- data %>%
    dplyr::mutate(
      inferred_dept = dplyr::case_when(
        # Pediatrics: Age-based
        !is.na(!!rlang::sym(age_col)) & !!rlang::sym(age_col) < 18 ~ "Pediatrics",

        # ICU: specimen-based (invasive)
        !is.na(!!rlang::sym(specimen_col)) &
          grepl("blood|csf|broncho|balf", tolower(!!rlang::sym(specimen_col))) ~ "ICU",

        # Surgery: specimen-based (sterile sites)
        !is.na(!!rlang::sym(specimen_col)) &
          grepl("peritoneal|abscess|tissue|wound", tolower(!!rlang::sym(specimen_col))) ~ "Surgery",

        # Obstetrics: diagnosis-based
        !is.na(!!rlang::sym(diagnosis_col)) &
          grepl("pregnancy|obstetric|maternal", tolower(!!rlang::sym(diagnosis_col))) ~ "Obstetrics",

        # Default
        TRUE ~ NA_character_
      ),
      !!department_col := dplyr::case_when(
        overwrite & !is.na(inferred_dept) ~ inferred_dept,
        is.na(!!rlang::sym(department_col)) & !is.na(inferred_dept) ~ inferred_dept,
        TRUE ~ !!rlang::sym(department_col)
      ),
      department_method = dplyr::case_when(
        !is.na(inferred_dept) ~ "inferred_heuristic",
        !is.na(!!rlang::sym(department_col)) ~ "provided",
        TRUE ~ NA_character_
      ),
      department_confidence = dplyr::case_when(
        department_method == "provided" ~ "high",
        department_method == "inferred_heuristic" ~ "low",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-inferred_dept)

  n_after_missing <- sum(is.na(data[[department_col]]))
  n_enriched <- n_before_missing - n_after_missing

  if (n_enriched > 0) {
    message(sprintf(
      "Enriched department: %d rows filled (LOW confidence - heuristic)",
      n_enriched
    ))
  }

  # Summary
  dept_summary <- data %>%
    dplyr::count(!!rlang::sym(department_col), department_confidence) %>%
    dplyr::arrange(dplyr::desc(n))

  message("\nDepartment distribution:")
  print(dept_summary)

  return(data)
}


#' Groom Optional Columns
#'
#' Cleans and standardizes all optional columns in a single pass.
#' Handles value normalization, whitespace trimming, and consistency checks.
#'
#' @param data Data frame
#' @param optional_cols Character vector. Optional column names to groom.
#'   If NULL, grooms all known optional columns. Default NULL.
#'
#' @return Data frame with groomed optional columns
#' @export
#'
#' @examples
#' \dontrun{
#' data_groomed <- groom_optional_columns(data)
#'
#' # Groom specific columns
#' data_groomed <- groom_optional_columns(
#'   data,
#'   optional_cols = c("hospital_department", "unit_type", "comorbidities")
#' )
#' }
groom_optional_columns <- function(data,
                                   optional_cols = NULL) {
  # Default optional columns
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

  message(sprintf(
    "Grooming %d optional columns: %s",
    length(optional_cols),
    paste(optional_cols[1:min(3, length(optional_cols))], collapse = ", ")
  ))

  # Groom each column
  for (col in optional_cols) {
    if (is.character(data[[col]])) {
      # Trim whitespace
      data[[col]] <- trimws(data[[col]])

      # Replace empty strings with NA
      data[[col]][data[[col]] == ""] <- NA_character_

      # Title case for consistency (except specific columns)
      if (!col %in% c("comorbidities", "previous_history")) {
        data[[col]] <- tools::toTitleCase(tolower(data[[col]]))
      }

      # Column-specific standardization
      if (col == "infection_type") {
        data[[col]] <- dplyr::case_when(
          grepl("^comm|^cai", tolower(data[[col]])) ~ "CAI",
          grepl("^hosp|^hai|^nosoco", tolower(data[[col]])) ~ "HAI",
          TRUE ~ data[[col]]
        )
      }

      if (col == "unit_type") {
        data[[col]] <- dplyr::case_when(
          grepl("icu|intensive|critical", tolower(data[[col]])) ~ "ICU",
          grepl("ward|general", tolower(data[[col]])) ~ "Ward",
          grepl("er|emergency", tolower(data[[col]])) ~ "Emergency",
          TRUE ~ data[[col]]
        )
      }
    }
  }

  # Summary of completeness
  completeness <- sapply(optional_cols, function(col) {
    sum(!is.na(data[[col]])) / nrow(data) * 100
  })

  completeness_df <- data.frame(
    column = names(completeness),
    completeness_pct = as.numeric(completeness),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(completeness_pct))

  message("\nOptional column completeness:")
  print(completeness_df)

  return(data)
}
