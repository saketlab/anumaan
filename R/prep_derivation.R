# prep_derivation.R
# Layer 5: Derivation - age, LOS, department, and helper functions
#
# Functions:
#   - find_extdata_file (internal helper)
#   - get_age_bins
#   - prep_assign_age_bins
#   - prep_fill_age
#   - prep_derive_los_from_dates  (canonical LOS derivation; merged prep_fill_los)
#   - prep_infer_department
#   - prep_derive_dob_from_components


#' Find Package Data File
#'
#' Locates data files in \code{inst/extdata}. Works both when the package is
#' installed and during development.
#'
#' @param filename Character. File name (e.g., "organisms.csv").
#' @return Character. Full path to the file, or empty string if not found.
#' @keywords internal
#' @noRd
find_extdata_file <- function(filename) {
  file_path <- system.file("extdata", filename, package = "anumaan")
  if (file_path != "" && file.exists(file_path)) return(file_path)

  pkg_root <- tryCatch(
    normalizePath(system.file(package = "anumaan")),
    error = function(e) NULL
  )
  if (!is.null(pkg_root)) {
    extdata_path <- file.path(pkg_root, "inst", "extdata", filename)
    if (file.exists(extdata_path)) return(normalizePath(extdata_path))
  }

  possible_roots <- c(
    ".", "..", "../..", "../../..",
    "../../../..", "anumaan", "../anumaan",
    "../../anumaan", "../../../anumaan"
  )
  for (root in possible_roots) {
    extdata_path <- file.path(root, "inst", "extdata", filename)
    if (file.exists(extdata_path)) return(normalizePath(extdata_path))
  }

  return("")
}


#' Get Age Bin Labels
#'
#' Returns predefined age bin label sets for use with \code{prep_assign_age_bins()}.
#'
#' @param type Character. One of \code{"GBD_standard"} (5-year bins),
#'   \code{"pediatric"}, \code{"geriatric"}, or \code{"neonatal"}.
#'   The \code{"neonatal"} preset uses fractional-year labels corresponding to:
#'   \code{<0.02} (0-7 days), \code{0.02-0.08} (7-28 days),
#'   \code{0.08-0.25} (28-90 days), \code{0.25-1} (3 months-1 year),
#'   then standard pediatric bands through adulthood.
#'
#' @return Character vector of age bin labels.
#' @export
get_age_bins <- function(type = "GBD_standard") {
  switch(type,
    GBD_standard = c(
      "<1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30",
      "30-35", "35-40", "40-45", "45-50", "50-55", "55-60",
      "60-65", "65-70", "70-75", "75-80", "80-85", "85+"
    ),
    pediatric = c(
      "0-0.08", "0.08-1", "1-2", "2-5", "5-10", "10-15", "15-18", "18+"
    ),
    geriatric = c(
      "0-50", "50-60", "60-65", "65-70", "70-75", "75-80",
      "80-85", "85-90", "90+"
    ),
    neonatal = c(
      "<0.02", "0.02-0.08", "0.08-0.25", "0.25-1",
      "1-5", "5-10", "10-15", "15-18", "18+"
    ),
    stop("Unknown age bin type. Use 'GBD_standard', 'pediatric', 'geriatric', or 'neonatal'")
  )
}


#' Assign Age Bins
#'
#' Categorizes age into bins for stratification using GBD standard or custom
#' break points. Supports ages stored as a single column in years, months, or
#' days, as well as compound datasets where years, months, and days are in
#' separate columns (e.g. neonatal records).
#'
#' When \code{age_months_col} or \code{age_days_col} are supplied they are
#' treated as the fractional remainder on top of \code{age_col} (e.g.
#' age_years = 0, age_months = 3, age_days = 5 -> 0.27 years). Missing values
#' in the component columns are treated as zero so that a partially-recorded
#' age is still binnable.
#'
#' @param data Data frame with an age column.
#' @param age_col Character. Primary age column. Default \code{"Age"}.
#' @param bins Character or character vector. Preset name (\code{"GBD_standard"},
#'   \code{"pediatric"}, \code{"geriatric"}, \code{"neonatal"}) or a custom
#'   vector of bin labels. Default \code{"GBD_standard"}.
#' @param age_unit Character. Unit of \code{age_col}: \code{"years"} (default),
#'   \code{"months"}, or \code{"days"}. The value is converted to decimal years
#'   before binning.
#' @param age_months_col Character or NULL. Optional column holding the months
#'   component of age (0-11). Added to \code{age_col} after unit conversion.
#' @param age_days_col Character or NULL. Optional column holding the days
#'   component of age (0-30). Added to \code{age_col} after unit conversion.
#' @param negative_age_strategy Character. How to handle negative ages in
#'   \code{age_col}. \code{"fallback"} (default) tries alternate columns /
#'   DOB logic and falls back to 0 years when unresolved; \code{"na"} sets
#'   negative ages to \code{NA}.
#' @param fallback_years_col Character or NULL. Optional years column used when
#'   \code{age_col} is negative. Default \code{"year"}.
#' @param fallback_months_col Character or NULL. Optional months column used
#'   when \code{age_col} is negative. Default \code{"months"}.
#' @param fallback_days_col Character or NULL. Optional days column used when
#'   \code{age_col} is negative. Default \code{"age_days"}.
#' @param fallback_dob_col Character or NULL. Optional DOB column used with
#'   \code{fallback_admission_col} to derive age when \code{age_col} is
#'   negative. Default \code{"dob"}.
#' @param fallback_admission_col Character or NULL. Optional admission/reference
#'   date column used with \code{fallback_dob_col}. Default
#'   \code{"admission_date"}.
#'
#' @return Data frame with \code{Age_bin} factor column added.
#' @export
prep_assign_age_bins <- function(data,
                                  age_col        = "Age",
                                  bins           = "GBD_standard",
                                  age_unit       = "years",
                                  age_months_col = NULL,
                                  age_days_col   = NULL,
                                  negative_age_strategy = c("fallback", "na"),
                                  fallback_years_col = "year",
                                  fallback_months_col = "months",
                                  fallback_days_col = "age_days",
                                  fallback_dob_col = "dob",
                                  fallback_admission_col = "admission_date") {
  if (!age_col %in% names(data)) stop(sprintf("Age column '%s' not found", age_col))

  age_unit <- match.arg(age_unit, c("years", "months", "days"))
  negative_age_strategy <- match.arg(negative_age_strategy)

  # Convert primary column to decimal years
  age_years <- switch(age_unit,
    years  = as.numeric(data[[age_col]]),
    months = as.numeric(data[[age_col]]) / 12,
    days   = as.numeric(data[[age_col]]) / 365.25
  )

  # Add months component if provided (treat NA as 0)
  if (!is.null(age_months_col)) {
    if (!age_months_col %in% names(data))
      stop(sprintf("age_months_col '%s' not found in data", age_months_col))
    m <- as.numeric(data[[age_months_col]])
    m[is.na(m)] <- 0
    age_years <- age_years + m / 12
  }

  # Add days component if provided (treat NA as 0)
  if (!is.null(age_days_col)) {
    if (!age_days_col %in% names(data))
      stop(sprintf("age_days_col '%s' not found in data", age_days_col))
    d <- as.numeric(data[[age_days_col]])
    d[is.na(d)] <- 0
    age_years <- age_years + d / 365.25
  }

  # Optional recovery for negative ages to avoid data loss in noisy datasets.
  neg_idx <- which(!is.na(age_years) & age_years < 0)
  if (length(neg_idx) > 0L) {
    if (negative_age_strategy == "na") {
      age_years[neg_idx] <- NA_real_
    } else {
      has_col <- function(col) !is.null(col) && col %in% names(data)

      years_vec <- if (has_col(fallback_years_col)) as.numeric(data[[fallback_years_col]]) else rep(NA_real_, nrow(data))
      months_vec <- if (has_col(fallback_months_col)) as.numeric(data[[fallback_months_col]]) else rep(NA_real_, nrow(data))
      days_vec <- if (has_col(fallback_days_col)) as.numeric(data[[fallback_days_col]]) else rep(NA_real_, nrow(data))

      dob_vec <- if (has_col(fallback_dob_col)) suppressWarnings(as.Date(data[[fallback_dob_col]])) else rep(as.Date(NA), nrow(data))
      adm_vec <- if (has_col(fallback_admission_col)) suppressWarnings(as.Date(data[[fallback_admission_col]])) else rep(as.Date(NA), nrow(data))

      recovered <- 0L
      for (i in neg_idx) {
        replacement <- NA_real_

        # Priority 1: explicit years/months/days columns.
        if (!is.na(years_vec[i]) && years_vec[i] >= 0) {
          replacement <- years_vec[i]
          if (!is.na(months_vec[i]) && months_vec[i] > 0) replacement <- replacement + months_vec[i] / 12
          if (!is.na(days_vec[i]) && days_vec[i] > 0) replacement <- replacement + days_vec[i] / 365.25
        } else if (!is.na(months_vec[i]) && months_vec[i] > 0) {
          replacement <- months_vec[i] / 12
          if (!is.na(days_vec[i]) && days_vec[i] > 0) replacement <- replacement + days_vec[i] / 365.25
        } else if (!is.na(days_vec[i]) && days_vec[i] > 0) {
          replacement <- days_vec[i] / 365.25
        } else if (!is.na(years_vec[i]) && years_vec[i] == 0 &&
                   (is.na(months_vec[i]) || months_vec[i] == 0) &&
                   (is.na(days_vec[i]) || days_vec[i] == 0)) {
          replacement <- 0
        }

        # Priority 2: derive from DOB and admission/reference date.
        if (is.na(replacement) && !is.na(dob_vec[i]) && !is.na(adm_vec[i])) {
          candidate <- as.numeric(difftime(adm_vec[i], dob_vec[i], units = "days")) / 365.25
          if (!is.na(candidate) && candidate >= 0) replacement <- candidate
        }

        # Priority 3: final fallback to 0 years (bin into youngest group).
        if (is.na(replacement)) replacement <- 0

        if (!is.na(replacement)) {
          age_years[i] <- replacement
          recovered <- recovered + 1L
        }
      }

      message(sprintf(
        "Recovered %d/%d negative age value(s) using fallback strategy.",
        recovered, length(neg_idx)
      ))
    }
  }

  if (is.character(bins) && length(bins) == 1) {
    bin_labels <- get_age_bins(bins)
  } else {
    bin_labels <- bins
  }

  breaks <- parse_age_bin_labels(bin_labels)

  data$Age_bin <- cut(
    age_years,
    breaks         = breaks$breaks,
    labels         = breaks$labels,
    right          = FALSE,
    include.lowest = TRUE
  )

  n_binned   <- sum(!is.na(data$Age_bin))
  n_unbinned <- sum(is.na(data$Age_bin) & !is.na(age_years))
  message(sprintf("Assigned age bins: %d binned, %d unbinned", n_binned, n_unbinned))
  return(data)
}


#' Fill Missing Age Values
#'
#' Fills missing Age by computing from DOB and a reference date. Tracks
#' derivation method and confidence for each record.
#'
#' @param data Data frame.
#' @param age_col Character. Age column. Default "Age".
#' @param dob_col Character. Date of birth column. Default "DOB".
#' @param date_col Character. Reference date column. Default "date_of_culture".
#' @param overwrite Logical. Recalculate even when Age is present. Default FALSE.
#'
#' @return Data frame with Age enriched plus \code{age_method} and
#'   \code{age_confidence} columns.
#' @export
prep_fill_age <- function(data,
                           age_col   = "Age",
                           dob_col   = "DOB",
                           date_col  = "date_of_culture",
                           overwrite = FALSE) {
  has_age  <- age_col  %in% names(data)
  has_dob  <- dob_col  %in% names(data)
  has_date <- date_col %in% names(data)

  if (!has_age) data[[age_col]] <- NA_real_

  if (!"age_method" %in% names(data)) {
    data$age_method     <- NA_character_
    data$age_confidence <- NA_character_
  }

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
          overwrite & !is.na(calculated_age) ~ calculated_age,
          is.na(!!rlang::sym(age_col)) & !is.na(calculated_age) ~ calculated_age,
          TRUE ~ !!rlang::sym(age_col)
        ),
        age_method = dplyr::case_when(
          !is.na(calculated_age) & (overwrite | is.na(age_method)) ~ "calculated_from_dob",
          !is.na(age_method)                                        ~ age_method,
          !is.na(!!rlang::sym(age_col))                            ~ "provided",
          TRUE                                                       ~ NA_character_
        ),
        age_confidence = dplyr::case_when(
          age_method %in% c("calculated_from_dob", "provided") ~ "high",
          TRUE ~ age_confidence
        )
      ) %>%
      dplyr::select(-calculated_age)

    n_enriched <- n_before_missing - sum(is.na(data[[age_col]]))
    if (n_enriched > 0)
      message(sprintf("Enriched Age: %d rows filled using DOB calculation", n_enriched))
  } else {
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

  message("\nAge enrichment summary:")
  print(dplyr::count(data, age_method, age_confidence) %>% dplyr::arrange(dplyr::desc(n)))

  n_still_missing <- sum(is.na(data[[age_col]]))
  if (n_still_missing > 0)
    message(sprintf("[!] %d rows still missing Age (%.1f%%)",
                    n_still_missing, 100 * n_still_missing / nrow(data)))

  return(data)
}


#' Infer Hospital Department
#'
#' Heuristically infers hospital department from specimen type, diagnosis,
#' and patient age. Low-confidence inference; used only to fill missing values.
#'
#' @param data Data frame.
#' @param department_col Character. Department column. Default "hospital_department".
#' @param specimen_col Character. Specimen type column. Default "specimen_type".
#' @param diagnosis_col Character. Diagnosis column. Default "diagnosis_1".
#' @param age_col Character. Age column. Default "Age".
#' @param overwrite Logical. Overwrite existing values. Default FALSE.
#'
#' @return Data frame with \code{hospital_department} enriched.
#' @export
prep_infer_department <- function(data,
                                   department_col = "hospital_department",
                                   specimen_col   = "specimen_type",
                                   diagnosis_col  = "diagnosis_1",
                                   age_col        = "Age",
                                   overwrite      = FALSE) {
  has_dept      <- department_col %in% names(data)
  has_specimen  <- specimen_col   %in% names(data)
  has_diagnosis <- diagnosis_col  %in% names(data)
  has_age       <- age_col        %in% names(data)

  if (!has_dept)      data[[department_col]] <- NA_character_
  if (!has_age)       data[[age_col]]        <- NA_real_
  if (!has_specimen)  data[[specimen_col]]   <- NA_character_
  if (!has_diagnosis) data[[diagnosis_col]]  <- NA_character_

  if (!has_specimen && !has_diagnosis && !has_age) {
    message("[!] Cannot infer department: no contextual data available")
    return(data)
  }

  n_before_missing <- sum(is.na(data[[department_col]]))
  message("Inferring hospital department from contextual data...")

  data <- data %>%
    dplyr::mutate(
      inferred_dept = dplyr::case_when(
        !is.na(!!rlang::sym(age_col)) & !!rlang::sym(age_col) < 18 ~ "Pediatrics",
        !is.na(!!rlang::sym(specimen_col)) &
          grepl("blood|csf|broncho|balf", tolower(!!rlang::sym(specimen_col))) ~ "ICU",
        !is.na(!!rlang::sym(specimen_col)) &
          grepl("peritoneal|abscess|tissue|wound", tolower(!!rlang::sym(specimen_col))) ~ "Surgery",
        !is.na(!!rlang::sym(diagnosis_col)) &
          grepl("pregnancy|obstetric|maternal", tolower(!!rlang::sym(diagnosis_col))) ~ "Obstetrics",
        TRUE ~ NA_character_
      ),
      !!department_col := dplyr::case_when(
        overwrite & !is.na(inferred_dept)                                     ~ inferred_dept,
        is.na(!!rlang::sym(department_col)) & !is.na(inferred_dept)           ~ inferred_dept,
        TRUE                                                                    ~ !!rlang::sym(department_col)
      ),
      department_method = dplyr::case_when(
        !is.na(inferred_dept)                    ~ "inferred_heuristic",
        !is.na(!!rlang::sym(department_col))     ~ "provided",
        TRUE                                      ~ NA_character_
      ),
      department_confidence = dplyr::case_when(
        department_method == "provided"           ~ "high",
        department_method == "inferred_heuristic" ~ "low",
        TRUE                                      ~ NA_character_
      )
    ) %>%
    dplyr::select(-inferred_dept)

  n_enriched <- n_before_missing - sum(is.na(data[[department_col]]))
  if (n_enriched > 0)
    message(sprintf("Enriched department: %d rows filled (LOW confidence - heuristic)", n_enriched))

  message("\nDepartment distribution:")
  print(dplyr::arrange(dplyr::count(data, !!rlang::sym(department_col), department_confidence),
                       dplyr::desc(n)))

  return(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 5)
# ---------------------------------------------------------------------------

#' Derive Length of Stay from Date Columns
#'
#' Computes LOS when missing (or always when \code{overwrite = TRUE}).
#' Adds \code{los_method} and \code{los_confidence} audit columns.
#' For AIIMS data, uses \code{unit_admission_date} + \code{unit_duration_days}
#' as a fallback when primary date columns are absent.
#'
#' Priority order:
#' \enumerate{
#'   \item \code{outcome_col} - \code{admission_col}
#'   \item \code{unit_admission_date} + \code{unit_duration_days} (AIIMS fallback)
#'   \item Already-present \code{los_col} value retained as-is
#' }
#'
#' @param data Data frame.
#' @param admission_col Character. Admission date column. Default "admission_date".
#' @param outcome_col Character. Outcome/discharge date column. Default "outcome_date".
#' @param los_col Character. Target LOS column to populate. Default "los_days".
#' @param unit_admission_col Character. AIIMS unit admission date fallback.
#'   Default "unit_admission_date".
#' @param unit_duration_col Character. AIIMS unit duration fallback.
#'   Default "unit_duration_days".
#' @param overwrite Logical. Recalculate even when \code{los_col} is already present.
#'   Default FALSE.
#'
#' @return Data frame with \code{los_col}, \code{los_method}, and
#'   \code{los_confidence} populated where possible.
#' @export
prep_derive_los_from_dates <- function(data,
                                        admission_col      = "admission_date",
                                        outcome_col        = "outcome_date",
                                        los_col            = "los_days",
                                        unit_admission_col = "unit_admission_date",
                                        unit_duration_col  = "unit_duration_days",
                                        overwrite          = FALSE) {
  if (!los_col %in% names(data)) data[[los_col]] <- NA_real_
  if (!"los_method"     %in% names(data)) data$los_method     <- NA_character_
  if (!"los_confidence" %in% names(data)) data$los_confidence <- NA_character_

  # Mark already-provided values
  provided <- !is.na(data[[los_col]])
  data$los_method[provided & is.na(data$los_method)]         <- "provided"
  data$los_confidence[provided & is.na(data$los_confidence)] <- "high"

  rows_missing <- if (overwrite) rep(TRUE, nrow(data)) else is.na(data[[los_col]])
  n_missing    <- sum(rows_missing)

  if (n_missing == 0L) {
    message(sprintf("[prep_derive_los_from_dates] %s already complete.", los_col))
    return(data)
  }

  # Primary: outcome_col - admission_col
  if (all(c(admission_col, outcome_col) %in% names(data))) {
    can_calc <- rows_missing & !is.na(data[[admission_col]]) & !is.na(data[[outcome_col]])
    if (any(can_calc)) {
      data[[los_col]][can_calc]     <- as.numeric(
        difftime(data[[outcome_col]][can_calc], data[[admission_col]][can_calc], units = "days"))
      data$los_method[can_calc]     <- "calculated_from_dates"
      data$los_confidence[can_calc] <- "high"
      message(sprintf("[prep_derive_los_from_dates] %d rows filled from %s - %s.",
                      sum(can_calc), outcome_col, admission_col))
      rows_missing <- if (overwrite) rep(FALSE, nrow(data)) else is.na(data[[los_col]])
    }
  }

  # Fallback (AIIMS): unit_duration_days
  if (any(rows_missing) &&
      all(c(unit_admission_col, unit_duration_col) %in% names(data))) {
    can_calc <- rows_missing & !is.na(data[[unit_admission_col]]) &
      !is.na(data[[unit_duration_col]])
    if (any(can_calc)) {
      data[[los_col]][can_calc]     <- as.numeric(data[[unit_duration_col]][can_calc])
      data$los_method[can_calc]     <- "aiims_unit_duration"
      data$los_confidence[can_calc] <- "medium"
      message(sprintf("[prep_derive_los_from_dates] %d rows filled from AIIMS unit duration.",
                      sum(can_calc)))
    }
  }

  n_negative <- sum(!is.na(data[[los_col]]) & data[[los_col]] < 0, na.rm = TRUE)
  if (n_negative > 0)
    warning(sprintf("[prep_derive_los_from_dates] %d row(s) have negative LOS.", n_negative))

  n_filled <- n_missing - sum(is.na(data[[los_col]]))
  message(sprintf("[prep_derive_los_from_dates] %d filled; %d remain missing.",
                  n_filled, sum(is.na(data[[los_col]]))))

  los_non_missing <- dplyr::filter(data, !is.na(!!rlang::sym(los_col)))
  message("\nLOS summary:")
  if (nrow(los_non_missing) == 0L) {
    print(data.frame(
      mean_los = NA_real_,
      median_los = NA_real_,
      min_los = NA_real_,
      max_los = NA_real_
    ))
  } else {
    print(dplyr::summarise(
      los_non_missing,
      mean_los   = mean(!!rlang::sym(los_col),   na.rm = TRUE),
      median_los = median(!!rlang::sym(los_col), na.rm = TRUE),
      min_los    = min(!!rlang::sym(los_col),    na.rm = TRUE),
      max_los    = max(!!rlang::sym(los_col),    na.rm = TRUE)
    ))
  }

  return(data)
}


#' Derive Date of Birth from Age Components
#'
#' Stewardship-specific helper: DOB is sometimes stored as separate
#' year/month/day components or as a decimal age in years. This function
#' reconstructs an approximate DOB or derives age directly when DOB is absent.
#'
#' Strategies (applied in order):
#' \enumerate{
#'   \item If \code{dob_year_col}, \code{dob_month_col}, and \code{dob_day_col}
#'     are present, assemble into a Date.
#'   \item If \code{age_years_col} is present and a reference date exists,
#'     estimate DOB as \code{reference_date - age_years * 365.25}.
#' }
#'
#' @param data Data frame.
#' @param dob_year_col Character. Year component column. Default "dob_year".
#' @param dob_month_col Character. Month component column. Default "dob_month".
#' @param dob_day_col Character. Day component column. Default "dob_day".
#' @param age_years_col Character. Decimal age in years. Default "age_years".
#' @param reference_date_col Character. Reference date for age-based DOB
#'   estimation. Default "admission_date".
#' @param dob_output_col Character. Output DOB column. Default "dob".
#'
#' @return Data frame with \code{dob} column added/populated.
#' @export
prep_derive_dob_from_components <- function(data,
                                             dob_year_col       = "dob_year",
                                             dob_month_col      = "dob_month",
                                             dob_day_col        = "dob_day",
                                             age_years_col      = "age_years",
                                             reference_date_col = "admission_date",
                                             dob_output_col     = "dob") {
  if (!dob_output_col %in% names(data)) data[[dob_output_col]] <- as.Date(NA_character_)

  rows_missing <- is.na(data[[dob_output_col]])
  n_missing    <- sum(rows_missing)
  if (n_missing == 0L) {
    message(sprintf("[prep_derive_dob_from_components] %s already complete.", dob_output_col))
    return(data)
  }

  # Strategy 1: year + month + day components
  component_cols <- c(dob_year_col, dob_month_col, dob_day_col)
  if (all(component_cols %in% names(data))) {
    can_build <- rows_missing &
      !is.na(data[[dob_year_col]]) &
      !is.na(data[[dob_month_col]]) &
      !is.na(data[[dob_day_col]])

    if (any(can_build)) {
      dob_str <- sprintf("%04d-%02d-%02d",
                         as.integer(data[[dob_year_col]][can_build]),
                         as.integer(data[[dob_month_col]][can_build]),
                         as.integer(data[[dob_day_col]][can_build]))
      data[[dob_output_col]][can_build] <- suppressWarnings(as.Date(dob_str))
      n_built <- sum(!is.na(data[[dob_output_col]][can_build]))
      message(sprintf("[prep_derive_dob_from_components] %d DOBs assembled from year/month/day.", n_built))
      rows_missing <- is.na(data[[dob_output_col]])
    }
  }

  # Strategy 2: age in years + reference date
  if (any(rows_missing) &&
      age_years_col %in% names(data) &&
      reference_date_col %in% names(data)) {
    can_estimate <- rows_missing &
      !is.na(data[[age_years_col]]) &
      !is.na(data[[reference_date_col]])

    if (any(can_estimate)) {
      est_dob <- data[[reference_date_col]][can_estimate] -
        round(as.numeric(data[[age_years_col]][can_estimate]) * 365.25)
      data[[dob_output_col]][can_estimate] <- est_dob
      message(sprintf("[prep_derive_dob_from_components] %d DOBs estimated from age + reference date.",
                      sum(can_estimate)))
    }
  }

  n_filled <- n_missing - sum(is.na(data[[dob_output_col]]))
  message(sprintf("[prep_derive_dob_from_components] %d filled; %d remain missing.", n_filled,
                  sum(is.na(data[[dob_output_col]]))))
  return(data)
}
