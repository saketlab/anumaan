# prep_ast.R
# Layer 4d: AST value harmonization
#
# Functions:
#   - prep_clean_ast_values       (canonical: extracts S/I/R from any raw string)
#   - prep_recode_intermediate_ast
#   - prep_harmonize_ast          (pipeline wrapper: clean → recode)
#   - prep_check_organism_ast_consistency
#   - prep_flag_invalid_ast
#   - prep_detect_duplicate_ast
#   - prep_resolve_duplicate_ast
#
# Note: prep_standardize_ast_values removed — exact-match subset of prep_clean_ast_values.


#' Clean Antibiotic Susceptibility Values
#'
#' Extracts clean S/I/R values from messy antibiotic result columns.
#' Handles common patterns like "S (HIGH LEVEL)", "R   Escherichia coli", etc.
#'
#' @param data Data frame with antibiotic susceptibility data.
#' @param value_col Character. Column containing susceptibility values.
#'   Default "antibiotic_value".
#' @param strict Logical. If TRUE, only accept S/I/R values. If FALSE, attempt
#'   to parse from messy strings. Default FALSE.
#'
#' @return Data frame with cleaned \code{antibiotic_value} column.
#' @export
prep_clean_ast_values <- function(data,
                                  value_col = "antibiotic_value",
                                  strict    = FALSE) {
  if (!value_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", value_col))
  }

  n_before         <- nrow(data)
  n_unique_before  <- dplyr::n_distinct(data[[value_col]], na.rm = TRUE)

  message(sprintf("Cleaning antibiotic values: %d unique values found", n_unique_before))

  data$temp_value_clean <- trimws(as.character(data[[value_col]]))

  if (strict) {
    data[[value_col]] <- ifelse(
      data$temp_value_clean %in% c("S", "I", "R"),
      data$temp_value_clean,
      NA_character_
    )
  } else {
    data[[value_col]] <- sapply(data$temp_value_clean, function(val) {
      if (is.na(val) || val == "") return(NA_character_)

      val_upper  <- toupper(val)
      first_char <- substr(val_upper, 1, 1)
      if (first_char %in% c("S", "I", "R")) return(first_char)

      if (grepl("^R[\\s(]|^R$", val_upper, perl = TRUE)) return("R")
      if (grepl("^S[\\s(]|^S$", val_upper, perl = TRUE)) return("S")
      if (grepl("^I[\\s(]|^I$", val_upper, perl = TRUE)) return("I")

      if (grepl("resistant",   val_upper)) return("R")
      if (grepl("susceptible|sensitive", val_upper)) return("S")
      if (grepl("intermediate", val_upper)) return("I")

      NA_character_
    }, USE.NAMES = FALSE)
  }

  data$temp_value_clean <- NULL

  n_unique_after <- dplyr::n_distinct(data[[value_col]], na.rm = TRUE)
  n_na           <- sum(is.na(data[[value_col]]))

  message(sprintf("Cleaned: %d unique values -> %d (S/I/R)", n_unique_before, n_unique_after))
  if (n_na > 0)
    message(sprintf("[!] %d values could not be parsed (%.1f%%)", n_na, 100 * n_na / n_before))

  message("\nValue distribution:")
  print(table(data[[value_col]], useNA = "ifany"))

  return(data)
}


#' Recode Intermediate (I) Susceptibility Values
#'
#' Converts "I" (Intermediate) values to either "S" or "R" based on antibiotic.
#' For Colistin: I -> S (clinical guideline).
#' For all other antibiotics: I -> R (conservative surveillance approach).
#'
#' @param data Data frame with antibiotic susceptibility data.
#' @param antibiotic_col Character. Column with antibiotic names.
#'   Default "antibiotic_name".
#' @param value_col Character. Column with S/I/R values.
#'   Default "antibiotic_value".
#' @param colistin_to_s Logical. Convert Colistin I to S. Default TRUE.
#' @param others_to_r Logical. Convert other antibiotics' I to R. Default TRUE.
#'
#' @return Data frame with recoded intermediate values.
#' @export
prep_recode_intermediate_ast <- function(data,
                                         antibiotic_col = "antibiotic_name",
                                         value_col      = "antibiotic_value",
                                         colistin_to_s  = TRUE,
                                         others_to_r    = TRUE) {
  if (!antibiotic_col %in% names(data)) stop(sprintf("Column '%s' not found", antibiotic_col))
  if (!value_col      %in% names(data)) stop(sprintf("Column '%s' not found", value_col))

  n_intermediate_before <- sum(data[[value_col]] == "I", na.rm = TRUE)

  if (n_intermediate_before == 0) {
    message("No intermediate (I) values found. Returning data unchanged.")
    return(data)
  }

  data$temp_value      <- data[[value_col]]
  data$temp_antibiotic <- tolower(trimws(data[[antibiotic_col]]))

  if (colistin_to_s) {
    idx <- which(data$temp_value == "I" &
                   grepl("colistin", data$temp_antibiotic, ignore.case = TRUE))
    if (length(idx) > 0) {
      data$temp_value[idx] <- "S"
      message(sprintf("Recoded %d Colistin I -> S", length(idx)))
    }
  }

  if (others_to_r) {
    idx <- which(data$temp_value == "I" &
                   !grepl("colistin", data$temp_antibiotic, ignore.case = TRUE))
    if (length(idx) > 0) {
      data$temp_value[idx] <- "R"
      message(sprintf("Recoded %d non-Colistin I -> R", length(idx)))
    }
  }

  data[[value_col]]    <- data$temp_value
  data$temp_value      <- NULL
  data$temp_antibiotic <- NULL

  n_intermediate_after <- sum(data[[value_col]] == "I", na.rm = TRUE)
  message(sprintf("Intermediate values: %d -> %d", n_intermediate_before, n_intermediate_after))

  return(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 4d)
# ---------------------------------------------------------------------------

#' Harmonize AST Values
#'
#' Wrapper that applies the full AST cleaning pipeline in sequence:
#' \enumerate{
#'   \item \code{prep_clean_ast_values()} — extract clean S/I/R from raw strings
#'     (handles exact matches, regex patterns, and keyword scanning in one pass)
#'   \item \code{prep_recode_intermediate_ast()} — resolve I values
#' }
#' Retains the original \code{ast_col} column as \code{ast_value_raw} and adds
#' \code{ast_value_harmonized} and \code{intermediate_recoded}.
#'
#' @param data Data frame with AST data.
#' @param antibiotic_col Character. Standardized antibiotic name column.
#'   Default "antibiotic_name_std".
#' @param ast_col Character. Raw AST value column. Default "ast_value_raw".
#' @param colistin_to_s Logical. Recode Colistin I -> S. Default TRUE.
#' @param others_to_r Logical. Recode non-Colistin I -> R. Default TRUE.
#'
#' @return Data frame with \code{ast_value_harmonized} and
#'   \code{intermediate_recoded} columns added.
#' @export
prep_harmonize_ast <- function(data,
                                antibiotic_col = "antibiotic_name_std",
                                ast_col        = "ast_value_raw",
                                colistin_to_s  = TRUE,
                                others_to_r    = TRUE) {
  if (!ast_col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping AST harmonization.", ast_col))
    data$ast_value_harmonized <- NA_character_
    data$intermediate_recoded <- FALSE
    return(data)
  }

  # Preserve raw values
  if (!"ast_value_raw" %in% names(data))
    data$ast_value_raw <- data[[ast_col]]

  work_col <- "ast_value_harmonized"

  # Step 1: clean messy strings -> S/I/R (handles exact, regex, and keyword patterns)
  data_work <- data
  data_work <- prep_clean_ast_values(data_work, value_col = ast_col, strict = FALSE)
  data[[work_col]] <- data_work[[ast_col]]

  # Track which rows had I before recoding
  had_intermediate <- !is.na(data[[work_col]]) & data[[work_col]] == "I"

  # Step 2: recode I values
  abx_col_available <- antibiotic_col %in% names(data)

  if (abx_col_available) {
    data_work2 <- data
    names(data_work2)[names(data_work2) == work_col] <- "ast_value_recode_tmp"
    data_work2 <- prep_recode_intermediate_ast(
      data_work2,
      antibiotic_col = antibiotic_col,
      value_col      = "ast_value_recode_tmp",
      colistin_to_s  = colistin_to_s,
      others_to_r    = others_to_r
    )
    data[[work_col]] <- data_work2$ast_value_recode_tmp
  } else {
    message("[prep_harmonize_ast] antibiotic column not found; skipping I-recoding step.")
  }

  data$intermediate_recoded <- had_intermediate & !is.na(data[[work_col]]) &
    data[[work_col]] != "I"

  n_harm    <- sum(!is.na(data[[work_col]]))
  n_recoded <- sum(data$intermediate_recoded, na.rm = TRUE)
  message(sprintf("[prep_harmonize_ast] %d/%d rows harmonized; %d intermediate values recoded.",
                  n_harm, nrow(data), n_recoded))

  return(data)
}


#' Check Organism–AST Consistency
#'
#' Flags impossible/implausible AST results where a drug is inappropriate
#' for the organism type (e.g., a Gram-positive organism tested against a
#' Gram-negative-only drug).
#'
#' Uses \code{inst/extdata/Organism_AST_drugs.csv} which maps organism groups
#' to expected antibiotic panels.
#'
#' @param data Data frame with organism and antibiotic columns.
#' @param organism_col Character. Normalized organism column.
#'   Default "organism_normalized".
#' @param antibiotic_col Character. Normalized antibiotic column.
#'   Default "antibiotic_normalized".
#'
#' @return Data frame with \code{is_ast_inconsistent} logical column added.
#' @export
prep_check_organism_ast_consistency <- function(data,
                                                 organism_col   = "organism_normalized",
                                                 antibiotic_col = "antibiotic_normalized") {
  data$is_ast_inconsistent <- FALSE

  missing_cols <- setdiff(c(organism_col, antibiotic_col), names(data))
  if (length(missing_cols) > 0) {
    warning(sprintf("[prep_check_organism_ast_consistency] Column(s) not found: %s. Skipping.",
                    paste(missing_cols, collapse = ", ")))
    data$is_ast_inconsistent <- NA
    return(data)
  }

  ref_path <- find_extdata_file("Organism_AST_drugs.csv")
  if (ref_path == "" || !file.exists(ref_path)) {
    warning("Organism_AST_drugs.csv not found. AST consistency check skipped.")
    data$is_ast_inconsistent <- NA
    return(data)
  }

  ast_ref <- readr::read_csv(ref_path, show_col_types = FALSE)

  # Expected columns: organism_group, antibiotic_name, expected (TRUE/FALSE or "yes"/"no")
  required_cols <- c("organism_group", "antibiotic_name")
  if (!all(required_cols %in% names(ast_ref))) {
    warning(sprintf("Organism_AST_drugs.csv must have columns: %s. Skipping.",
                    paste(required_cols, collapse = ", ")))
    data$is_ast_inconsistent <- NA
    return(data)
  }

  # Build exclusion table: organism_group + antibiotic combos that are NOT expected
  excluded <- ast_ref[!isTRUE(ast_ref$expected) & !ast_ref$expected %in% c("yes", "TRUE", TRUE), ]

  n_flagged <- 0L
  if ("organism_group" %in% names(data) && nrow(excluded) > 0) {
    org_lower <- tolower(trimws(data$organism_group))
    abx_lower <- tolower(trimws(data[[antibiotic_col]]))

    for (i in seq_len(nrow(excluded))) {
      grp_pat <- tolower(trimws(excluded$organism_group[i]))
      abx_pat <- tolower(trimws(excluded$antibiotic_name[i]))
      flag    <- !is.na(org_lower) & !is.na(abx_lower) &
        org_lower == grp_pat & abx_lower == abx_pat
      data$is_ast_inconsistent[flag] <- TRUE
      n_flagged <- n_flagged + sum(flag, na.rm = TRUE)
    }
  }

  message(sprintf("[prep_check_organism_ast_consistency] %d implausible organism-AST pairs flagged.",
                  n_flagged))
  return(data)
}


#' Flag Invalid AST Values
#'
#' Adds a logical column \code{is_ast_invalid} that is TRUE when the harmonized
#' AST value is not in \code{c("S", "I", "R")} and not NA.
#'
#' @param data Data frame.
#' @param col Character. Harmonized AST column. Default "ast_value_harmonized".
#'
#' @return Data frame with \code{is_ast_invalid} column added.
#' @export
prep_flag_invalid_ast <- function(data, col = "ast_value_harmonized") {
  if (!col %in% names(data)) {
    warning(sprintf("Column '%s' not found. All rows flagged as invalid.", col))
    data$is_ast_invalid <- NA
    return(data)
  }

  data$is_ast_invalid <- !is.na(data[[col]]) & !data[[col]] %in% c("S", "I", "R")

  n_invalid <- sum(data$is_ast_invalid, na.rm = TRUE)
  if (n_invalid > 0) {
    message(sprintf("[prep_flag_invalid_ast] %d row(s) with invalid AST value.", n_invalid))
    top_vals <- utils::head(sort(table(data[[col]][data$is_ast_invalid]), decreasing = TRUE), 5)
    message("Top invalid values:")
    print(top_vals)
  } else {
    message("[prep_flag_invalid_ast] All non-NA AST values are valid (S/I/R).")
  }

  return(data)
}


#' Detect Duplicate AST Results
#'
#' Identifies rows where the same patient, organism, antibiotic, and date
#' have conflicting AST results (e.g., one row shows R and another shows S).
#'
#' @param data Data frame with AST data in long format.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param date_col Character. Culture date column. Default "culture_date".
#' @param ast_col Character. Harmonized AST value column.
#'   Default "ast_value_harmonized".
#'
#' @return A data frame of conflict records (with a \code{conflict_group} column),
#'   or an empty data frame if no conflicts exist.
#' @export
prep_detect_duplicate_ast <- function(data,
                                       patient_col    = "patient_id",
                                       organism_col   = "organism_normalized",
                                       antibiotic_col = "antibiotic_normalized",
                                       date_col       = "culture_date",
                                       ast_col        = "ast_value_harmonized") {
  key_cols <- c(patient_col, organism_col, antibiotic_col, date_col, ast_col)
  missing  <- setdiff(key_cols, names(data))
  if (length(missing) > 0) {
    warning(sprintf("[prep_detect_duplicate_ast] Column(s) not found: %s.",
                    paste(missing, collapse = ", ")))
    return(data.frame())
  }

  group_cols <- c(patient_col, organism_col, antibiotic_col, date_col)

  conflict_groups <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_ast_values = dplyr::n_distinct(.data[[ast_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_ast_values > 1)

  if (nrow(conflict_groups) == 0) {
    message("[prep_detect_duplicate_ast] No duplicate AST conflicts detected.")
    return(data.frame())
  }

  conflicts <- dplyr::inner_join(data, conflict_groups, by = group_cols) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))

  conflicts$conflict_group <- as.integer(dplyr::group_indices(
    conflicts, dplyr::across(dplyr::all_of(group_cols))
  ))

  message(sprintf("[prep_detect_duplicate_ast] %d conflict group(s) found (%d rows).",
                  nrow(conflict_groups), nrow(conflicts)))
  return(conflicts)
}


#' Resolve Duplicate AST Results
#'
#' Resolves conflicting AST records for the same patient + organism +
#' antibiotic + date combination.
#'
#' Available strategies:
#' \describe{
#'   \item{resistant_wins}{Keep the R result; for R vs S, R takes precedence.}
#'   \item{susceptible_wins}{Keep the S result.}
#'   \item{first}{Keep the first row in the data for each conflict group.}
#' }
#'
#' @param data Data frame with AST data in long format.
#' @param strategy Character. One of \code{"resistant_wins"} (default),
#'   \code{"susceptible_wins"}, or \code{"first"}.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param date_col Character. Culture date column. Default "culture_date".
#' @param ast_col Character. AST value column. Default "ast_value_harmonized".
#'
#' @return Data frame with conflicts resolved (one row per unique key combination).
#' @export
prep_resolve_duplicate_ast <- function(data,
                                        strategy       = "resistant_wins",
                                        patient_col    = "patient_id",
                                        organism_col   = "organism_normalized",
                                        antibiotic_col = "antibiotic_normalized",
                                        date_col       = "culture_date",
                                        ast_col        = "ast_value_harmonized") {
  strategy <- match.arg(strategy, c("resistant_wins", "susceptible_wins", "first"))

  key_cols    <- c(patient_col, organism_col, antibiotic_col, date_col)
  all_needed  <- c(key_cols, ast_col)
  missing     <- setdiff(all_needed, names(data))
  if (length(missing) > 0) {
    warning(sprintf("[prep_resolve_duplicate_ast] Column(s) not found: %s. Returning data unchanged.",
                    paste(missing, collapse = ", ")))
    return(data)
  }

  n_before <- nrow(data)

  if (strategy == "resistant_wins") {
    priority <- c("R" = 1L, "I" = 2L, "S" = 3L)
    data$.ast_priority <- priority[data[[ast_col]]]
    data$.ast_priority[is.na(data$.ast_priority)] <- 99L

    data <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
      dplyr::slice_min(order_by = .data$.ast_priority, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    data$.ast_priority <- NULL

  } else if (strategy == "susceptible_wins") {
    priority <- c("S" = 1L, "I" = 2L, "R" = 3L)
    data$.ast_priority <- priority[data[[ast_col]]]
    data$.ast_priority[is.na(data$.ast_priority)] <- 99L

    data <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
      dplyr::slice_min(order_by = .data$.ast_priority, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    data$.ast_priority <- NULL

  } else {
    # first: keep first row per key
    data <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }

  n_after   <- nrow(data)
  n_removed <- n_before - n_after
  message(sprintf("[prep_resolve_duplicate_ast] Strategy '%s': removed %d duplicate row(s).",
                  strategy, n_removed))

  return(data)
}
