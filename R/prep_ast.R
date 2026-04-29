# prep_ast.R
# AST data reshaping, value harmonization, and quality control.
#
# Functions:
#   - prep_pivot_ast_wide_to_long  (wide -> long reshape)
#   - prep_create_wide_ast_matrix  (long -> wide reshape)
#   - prep_clean_ast_values        (extracts S/I/R from any raw string)
#   - prep_recode_intermediate_ast
#   - prep_harmonize_ast           (pipeline wrapper: clean -> recode)
#   - prep_check_organism_ast_consistency
#   - prep_flag_invalid_ast
#   - prep_deduplicate_ast        (mode="detect" -> flag + QC summary; mode="remove" -> flag + resolve)
#
# Note: prep_standardize_ast_values removed - exact-match subset of prep_clean_ast_values.


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
#'   \item \code{prep_clean_ast_values()} - extract clean S/I/R from raw strings
#'     (handles exact matches, regex patterns, and keyword scanning in one pass)
#'   \item \code{prep_recode_intermediate_ast()} - resolve I values
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


#' Check Organism-AST Consistency
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
  required_cols <- c("organism_group", "antibiotic_name")

  # The bundled CSV is a semi-structured report with a title row and grouped
  # organism names; normalize it into a simple expected-combination table.
  if (!all(required_cols %in% names(ast_ref))) {
    raw_ref <- readr::read_csv(ref_path, skip = 1, show_col_types = FALSE)
    names(raw_ref) <- gsub("[^a-z0-9]+", "_", tolower(names(raw_ref)))
    group_col <- intersect(names(raw_ref), c("organism_or_organism_group", "organism_group"))
    agents_col <- intersect(names(raw_ref), c("antimicrobial_agents", "antibiotic_name"))

    if (length(group_col) == 1L && length(agents_col) == 1L) {
      ast_ref <- raw_ref[, c(group_col, agents_col), drop = FALSE]
      names(ast_ref) <- required_cols
      ast_ref$organism_group <- trimws(as.character(ast_ref$organism_group))
      ast_ref$organism_group[ast_ref$organism_group == ""] <- NA_character_
      ast_ref$organism_group <- tidyr::fill(ast_ref, organism_group)$organism_group
      ast_ref <- tidyr::separate_longer_delim(
        ast_ref,
        cols = "antibiotic_name",
        delim = ","
      )
      ast_ref$antibiotic_name <- tolower(trimws(ast_ref$antibiotic_name))
      ast_ref$organism_group <- tolower(trimws(ast_ref$organism_group))
      ast_ref <- ast_ref[!is.na(ast_ref$organism_group) & !is.na(ast_ref$antibiotic_name) & nzchar(ast_ref$antibiotic_name), , drop = FALSE]
      ast_ref$expected <- TRUE
    } else {
      warning(sprintf("Organism_AST_drugs.csv must have columns: %s. Skipping.",
                      paste(required_cols, collapse = ", ")))
      data$is_ast_inconsistent <- NA
      return(data)
    }
  }

  n_flagged <- 0L
  if ("organism_group" %in% names(data) && nrow(ast_ref) > 0) {
    org_lower <- tolower(trimws(as.character(data$organism_group)))
    abx_lower <- tolower(trimws(as.character(data[[antibiotic_col]])))
    expected_keys <- unique(paste(ast_ref$organism_group, ast_ref$antibiotic_name, sep = "||"))
    row_keys <- paste(org_lower, abx_lower, sep = "||")
    can_check <- !is.na(org_lower) & nzchar(org_lower) & !is.na(abx_lower) & nzchar(abx_lower)
    data$is_ast_inconsistent[can_check] <- !row_keys[can_check] %in% expected_keys
    n_flagged <- sum(data$is_ast_inconsistent, na.rm = TRUE)
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


#' Handle Duplicate AST Results
#'
#' Detects and optionally resolves conflicting AST records for the same
#' patient + organism + antibiotic + date combination.
#'
#' \describe{
#'   \item{\code{"detect"}}{Flags conflicting rows with \code{is_ast_duplicate = TRUE}
#'     and prints a QC summary of all conflict groups. Returns the full data frame
#'     with the flag column so you can inspect or filter before deciding how to resolve.}
#'   \item{\code{"remove"}}{Runs the detect step first (flag + QC summary), then
#'     applies \code{strategy} to keep one row per key combination and drops the
#'     flag column from the returned data.}
#' }
#'
#' @param data Data frame with AST data in long format.
#' @param mode Character. \code{"detect"} (flag only) or \code{"remove"} (flag then resolve).
#'   Default \code{"detect"}.
#' @param strategy Character. Resolution strategy used only when \code{mode = "remove"}.
#'   One of \code{"resistant_wins"} (default), \code{"susceptible_wins"}, or \code{"first"}.
#' @param patient_col Character. Patient ID column. Default \code{"patient_id"}.
#' @param organism_col Character. Organism column. Default \code{"organism_normalized"}.
#' @param antibiotic_col Character. Antibiotic column. Default \code{"antibiotic_normalized"}.
#' @param date_col Character. Culture date column. Default \code{"culture_date"}.
#' @param ast_col Character. Harmonized AST value column. Default \code{"ast_value_harmonized"}.
#'
#' @return
#' \itemize{
#'   \item \code{mode = "detect"}: original data with \code{is_ast_duplicate} logical column added.
#'   \item \code{mode = "remove"}: data with conflicts resolved (one row per key) and no flag column.
#' }
#' @export
prep_deduplicate_ast <- function(data,
                                  mode           = c("detect", "remove"),
                                  strategy       = c("resistant_wins", "susceptible_wins", "first"),
                                  patient_col    = "patient_id",
                                  organism_col   = "organism_normalized",
                                  antibiotic_col = "antibiotic_normalized",
                                  date_col       = "culture_date",
                                  ast_col        = "ast_value_harmonized") {
  mode     <- match.arg(mode)
  strategy <- match.arg(strategy)

  key_cols   <- c(patient_col, organism_col, antibiotic_col, date_col)
  all_needed <- c(key_cols, ast_col)
  missing    <- setdiff(all_needed, names(data))
  if (length(missing) > 0) {
    warning(sprintf("[prep_deduplicate_ast] Column(s) not found: %s. Returning data unchanged.",
                    paste(missing, collapse = ", ")))
    return(data)
  }

  # --- Step 1: detect conflicts and flag rows ---
  conflict_summary <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::summarise(
      n_ast_values       = dplyr::n_distinct(.data[[ast_col]], na.rm = TRUE),
      conflicting_values = paste(sort(unique(.data[[ast_col]])), collapse = " vs "),
      n_rows             = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_ast_values > 1)

  data <- dplyr::left_join(
    data,
    conflict_summary %>% dplyr::select(dplyr::all_of(key_cols), n_ast_values),
    by = key_cols
  ) %>%
    dplyr::mutate(is_ast_duplicate = !is.na(.data$n_ast_values)) %>%
    dplyr::select(-n_ast_values)

  n_groups <- nrow(conflict_summary)
  n_flagged <- sum(data$is_ast_duplicate, na.rm = TRUE)

  if (n_groups == 0) {
    message("[prep_deduplicate_ast] No duplicate AST conflicts detected.")
  } else {
    message(sprintf("[prep_deduplicate_ast] %d conflict group(s) found, %d rows flagged.",
                    n_groups, n_flagged))
    message("\nConflict summary (top 10):")
    print(utils::head(conflict_summary[, c(key_cols, "conflicting_values", "n_rows")], 10))
  }

  if (mode == "detect") return(data)

  # --- Step 2 (remove only): resolve conflicts via strategy ---
  n_before <- nrow(data)

  if (strategy %in% c("resistant_wins", "susceptible_wins")) {
    priority <- if (strategy == "resistant_wins") c("R" = 1L, "I" = 2L, "S" = 3L)
                else                              c("S" = 1L, "I" = 2L, "R" = 3L)
    data$.ast_priority <- priority[data[[ast_col]]]
    data$.ast_priority[is.na(data$.ast_priority)] <- 99L

    data <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
      dplyr::slice_min(order_by = .data$.ast_priority, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    data$.ast_priority <- NULL

  } else {
    data <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }

  data$is_ast_duplicate <- NULL

  message(sprintf("[prep_deduplicate_ast] Strategy '%s': removed %d duplicate row(s).",
                  strategy, n_before - nrow(data)))
  return(data)
}


# ---------------------------------------------------------------------------
# Reshape helpers
# ---------------------------------------------------------------------------

#' Convert Wide Format to Long Format
#'
#' Converts wide format data (where each antibiotic is a column) to long format
#' (one row per organism-antibiotic combination). This is the first step before
#' normalization and analysis.
#'
#' @param data Data frame in wide format (antibiotics as columns)
#' @param antibiotic_cols Character vector. Names of antibiotic columns to pivot.
#'   If NULL, will auto-detect based on pattern. Default NULL.
#' @param pattern Character. Regex pattern to identify antibiotic columns if
#'   antibiotic_cols not provided. Default NULL (no auto-detection).
#' @param id_cols Character vector. Columns to keep as identifiers (not pivoted).
#'   Default c("patient_id", "event_id", "organism_name", "date_of_culture").
#' @param antibiotic_name_col Character. Name for the new column containing
#'   antibiotic names. Default "antibiotic_name".
#' @param antibiotic_value_col Character. Name for the new column containing
#'   susceptibility results. Default "antibiotic_value".
#' @param remove_missing Logical. Remove rows where antibiotic_value is NA, empty,
#'   or "-". Default TRUE.
#' @param create_event_id Logical. Create event_id column if it doesn't exist
#'   (uses row numbers). Default FALSE.
#'
#' @return Data frame in long format
#' @export
#'
#' @examples
#' \dontrun{
#' # Specify antibiotic columns explicitly
#' long_data <- prep_pivot_ast_wide_to_long(
#'   data = raw_data,
#'   antibiotic_cols = c("AMIKACIN", "GENTAMICIN", "CIPROFLOXACIN")
#' )
#'
#' # Auto-detect columns by pattern (columns 12-53)
#' long_data <- prep_pivot_ast_wide_to_long(
#'   data = raw_data,
#'   antibiotic_cols = names(raw_data)[12:53]
#' )
#'
#' # Auto-detect uppercase antibiotic names
#' long_data <- prep_pivot_ast_wide_to_long(
#'   data = raw_data,
#'   pattern = "^[A-Z]+$"
#' )
#' }
prep_pivot_ast_wide_to_long <- function(data,
                                        antibiotic_cols      = NULL,
                                        pattern              = NULL,
                                        id_cols              = c("patient_id", "event_id",
                                                                 "organism_name", "date_of_culture"),
                                        antibiotic_name_col  = "antibiotic_name",
                                        antibiotic_value_col = "antibiotic_value",
                                        remove_missing       = TRUE,
                                        create_event_id      = FALSE) {
  n_before <- nrow(data)

  if (create_event_id && !"event_id" %in% names(data)) {
    message("Creating event_id column from row numbers...")
    data$event_id <- seq_len(nrow(data))
  }

  if (is.null(antibiotic_cols)) {
    if (!is.null(pattern)) {
      antibiotic_cols <- names(data)[grepl(pattern, names(data))]
      message(sprintf("Auto-detected %d antibiotic columns using pattern '%s'",
                      length(antibiotic_cols), pattern))
    } else {
      stop("Either antibiotic_cols or pattern must be provided")
    }
  }

  if (length(antibiotic_cols) == 0)
    stop("No antibiotic columns found")

  id_cols <- intersect(id_cols, names(data))

  missing_cols <- setdiff(antibiotic_cols, names(data))
  if (length(missing_cols) > 0) {
    warning(sprintf("Some antibiotic columns not found: %s",
                    paste(missing_cols, collapse = ", ")))
    antibiotic_cols <- intersect(antibiotic_cols, names(data))
  }

  message(sprintf("Pivoting %d antibiotic columns to long format...", length(antibiotic_cols)))

  long_data <- data %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(antibiotic_cols),
      names_to  = antibiotic_name_col,
      values_to = antibiotic_value_col
    )

  message(sprintf("Pivoted: %d rows -> %d rows (wide -> long)", n_before, nrow(long_data)))

  if (remove_missing) {
    n_before_filter <- nrow(long_data)
    long_data <- long_data %>%
      dplyr::filter(
        !is.na(!!rlang::sym(antibiotic_value_col)),
        !!rlang::sym(antibiotic_value_col) != "",
        !!rlang::sym(antibiotic_value_col) != "-"
      )
    n_removed <- n_before_filter - nrow(long_data)
    if (n_removed > 0)
      message(sprintf("Removed %d rows with missing/empty values", n_removed))
  }

  n_unique_abx    <- dplyr::n_distinct(long_data[[antibiotic_name_col]])
  n_unique_events <- dplyr::n_distinct(long_data[["event_id"]], na.rm = TRUE)
  message(sprintf("Final long format: %d rows x %d columns (%d antibiotics, %d events)",
                  nrow(long_data), ncol(long_data), n_unique_abx, n_unique_events))

  return(long_data)
}


#' Create Wide Format AST Matrix
#'
#' Converts long format (one row per organism-antibiotic) to wide format
#' (one row per event with antibiotic columns). Useful for analysis and
#' machine learning applications.
#'
#' @param data Data frame in long format
#' @param event_col Character. Event ID column. Default "event_id".
#' @param antibiotic_col Character. Antibiotic/class column to pivot.
#'   Default "antibiotic_normalized".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param prefix Character. Prefix for pivoted columns. Default "abx_".
#' @param keep_cols Character vector. Additional columns to keep from original data.
#'   Default c("patient_id", "organism_normalized", "date_of_culture").
#'
#' @return Wide format data frame (one row per event)
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic wide format
#' wide_data <- prep_create_wide_ast_matrix(data)
#'
#' # Class-level wide format
#' wide_data <- prep_create_wide_ast_matrix(
#'   data,
#'   antibiotic_col = "antibiotic_class",
#'   prefix = "class_"
#' )
#' }
prep_create_wide_ast_matrix <- function(data,
                                        event_col          = "event_id",
                                        antibiotic_col     = "antibiotic_normalized",
                                        susceptibility_col = "antibiotic_value",
                                        prefix             = "abx_",
                                        keep_cols          = c("patient_id", "organism_normalized",
                                                               "date_of_culture")) {
  if (!event_col %in% names(data))
    stop(sprintf("Column '%s' not found", event_col))
  if (!antibiotic_col %in% names(data))
    stop(sprintf("Column '%s' not found", antibiotic_col))
  if (!susceptibility_col %in% names(data))
    stop(sprintf("Column '%s' not found", susceptibility_col))

  keep_cols <- unique(c(event_col, intersect(keep_cols, names(data))))
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  message(sprintf("Creating wide format: pivoting '%s' column...", antibiotic_col))

  wide_data <- data %>%
    dplyr::select(dplyr::all_of(keep_cols),
                  !!rlang::sym(antibiotic_col),
                  !!rlang::sym(susceptibility_col)) %>%
    dplyr::group_by(!!rlang::sym(event_col), !!rlang::sym(antibiotic_col)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      id_cols     = dplyr::all_of(keep_cols),
      names_from  = !!rlang::sym(antibiotic_col),
      values_from = !!rlang::sym(susceptibility_col),
      names_prefix = prefix
    )

  names(wide_data) <- gsub("[^A-Za-z0-9_]", "_", names(wide_data))
  names(wide_data) <- gsub("_{2,}", "_", names(wide_data))

  n_antibiotics <- ncol(wide_data) - length(keep_cols)
  message(sprintf("Created wide format: %d events x %d antibiotics",
                  nrow(wide_data), n_antibiotics))

  if (nrow(wide_data) != n_events_before)
    warning(sprintf("[!] Row count mismatch: %d events in input, %d rows in wide format",
                    n_events_before, nrow(wide_data)))

  return(wide_data)
}
