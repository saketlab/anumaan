# prep_dates.R
# Layer 3a: Date parsing + coercion
#
# Functions:
#   - prep_parse_date_column  (handles Excel serial, Unix ms, ISO/DMY/MDY, swapped digits)
#   - prep_coerce_dates       (auto-detects and converts all date-like columns)
#   - prep_validate_date_logic (6-check cross-column consistency validator)


#' Detect and decode encrypted or non-standard date columns
#'
#' Inspects a column for date-like content and automatically selects the
#' correct parsing strategy (Excel serial, Unix timestamp ms, ISO/DMY/MDY
#' strings, reversed YYYYDDMM).
#'
#' @param x Vector (character, numeric, or Date-like) to parse.
#' @param col_name Character. Column name used in messages.
#' @param table_label Character. Table label used in messages.
#'
#' @return A Date vector of the same length as \code{x}.
#' @export
prep_parse_date_column <- function(x, col_name = "date", table_label = "table") {
  n <- length(x)

  if (inherits(x, "Date"))    return(x)
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.Date(x))

  out   <- rep(as.Date(NA_character_), n)
  x_chr <- trimws(as.character(x))
  placeholders <- c("", "NULL", "null", "NA", "N/A", "None", "none", "nan", "NaN")
  x_chr[x_chr %in% placeholders] <- NA_character_

  idx_present <- which(!is.na(x_chr))
  if (length(idx_present) == 0L) {
    message(sprintf("[%s] '%s': all values missing - returning NA Date vector.",
                    table_label, col_name))
    return(out)
  }

  # ---- Numeric probe -------------------------------------------------------
  x_num      <- suppressWarnings(as.numeric(x_chr))
  is_numeric <- !is.na(x_num)

  # Excel serial dates: 10 000-99 999
  excel_idx <- idx_present[is_numeric[idx_present] & x_num[idx_present] > 10000 &
                              x_num[idx_present] < 100000]
  if (length(excel_idx) > 0L) {
    out[excel_idx] <- as.Date(x_num[excel_idx], origin = "1899-12-30")
    message(sprintf("[%s] '%s': %d value(s) decoded as Excel serial date.",
                    table_label, col_name, length(excel_idx)))
  }

  # Unix timestamps in milliseconds: > 1e10
  unix_idx <- idx_present[is_numeric[idx_present] & x_num[idx_present] > 1e10]
  if (length(unix_idx) > 0L) {
    out[unix_idx] <- as.Date(as.POSIXct(x_num[unix_idx] / 1000, origin = "1970-01-01"))
    message(sprintf("[%s] '%s': %d value(s) decoded as Unix timestamp (ms).",
                    table_label, col_name, length(unix_idx)))
  }

  # ---- String parsing for remaining indices --------------------------------
  # Keep unresolved numeric-looking strings here as well, because patterns like
  # YYYYDDMM (e.g. 20201301) are numeric but still need the swap fallback.
  remaining <- idx_present[is.na(out[idx_present])]

  if (length(remaining) > 0L) {
    parsed <- suppressWarnings(
      lubridate::parse_date_time(
        x_chr[remaining],
        orders = c(
          "Y-m-d", "d-m-Y", "m-d-Y",
          "Y/m/d", "d/m/Y", "m/d/Y",
          "Ymd",   "dmY",   "mdY",
          "Y-m-d H:M:S", "d-m-Y H:M:S", "m-d-Y H:M:S"
        ),
        quiet = TRUE
      )
    )
    parsed_date <- as.Date(parsed)
    succeeded   <- !is.na(parsed_date)
    out[remaining[succeeded]] <- parsed_date[succeeded]

    # Reversed YYYYDDMM / DDMMYYYY attempt on still-failing values
    still_bad <- remaining[!succeeded]
    if (length(still_bad) > 0L) {
      candidates  <- x_chr[still_bad]
      digits_only <- gsub("[^0-9]", "", candidates)
      swapped     <- ifelse(
        nchar(digits_only) == 8L,
        paste0(substr(digits_only, 1, 4), "-",
               substr(digits_only, 7, 8), "-",
               substr(digits_only, 5, 6)),
        NA_character_
      )
      swap_parsed <- suppressWarnings(as.Date(swapped, format = "%Y-%m-%d"))
      swap_ok     <- !is.na(swap_parsed)
      if (any(swap_ok)) {
        out[still_bad[swap_ok]] <- swap_parsed[swap_ok]
        message(sprintf("[%s] '%s': %d value(s) decoded via day/month swap.",
                        table_label, col_name, sum(swap_ok)))
        still_bad <- still_bad[!swap_ok]
      }
    }

    if (length(still_bad) > 0L) {
      encrypted_like <- nchar(x_chr[still_bad]) > 12 &
        !grepl("[/\\-]", x_chr[still_bad])
      n_encrypted <- sum(encrypted_like, na.rm = TRUE)
      n_unparsed  <- length(still_bad)

      if (n_encrypted > 0L) {
        warning(sprintf(
          "[%s] '%s': %d value(s) appear encrypted/undecodable and will be set to NA. Sample: %s",
          table_label, col_name, n_encrypted,
          paste(utils::head(x_chr[still_bad][encrypted_like], 3), collapse = " | ")
        ))
      } else if (n_unparsed > 0L) {
        warning(sprintf(
          "[%s] '%s': %d value(s) could not be parsed as dates and will be set to NA. Sample: %s",
          table_label, col_name, n_unparsed,
          paste(utils::head(x_chr[still_bad], 3), collapse = " | ")
        ))
      }
    }
  }

  n_decoded <- sum(!is.na(out[idx_present]))
  n_failed  <- length(idx_present) - n_decoded
  message(sprintf(
    "[%s] '%s': %d / %d non-missing value(s) successfully parsed as Date (%d failed).",
    table_label, col_name, n_decoded, length(idx_present), n_failed
  ))

  out
}


#' Detect and convert all date-like columns in a table
#'
#' Auto-detects columns with date/time-like names (or uses a supplied list),
#' then calls \code{prep_parse_date_column()} on each.
#'
#' @param data Data frame.
#' @param cols Character vector of column names to convert. When \code{NULL},
#'   columns whose names match a date/time pattern are detected automatically.
#' @param table_label Character. Label used in messages.
#'
#' @return Data frame with date columns converted to \code{Date}.
#' @export
prep_coerce_dates <- function(data, cols = NULL, table_label = "table") {
  if (is.null(cols)) {
    cols <- grep(
      "(^date_|_date$|date$|_date_|^fever_date|^date_HAI|_treat$|date_treat)",
      names(data),
      value = TRUE, ignore.case = TRUE
    )
  }
  cols <- intersect(cols, names(data))
  if (length(cols) == 0L) return(data)
  for (col in cols) {
    data[[col]] <- prep_parse_date_column(data[[col]], col_name = col, table_label = table_label)
  }
  data
}


# ---------------------------------------------------------------------------
# New functions (Layer 3a)
# ---------------------------------------------------------------------------

#' Validate Date Logic
#'
#' Checks that date sequences and derived values are internally consistent.
#' Consolidates all cross-column date/age checks in one place:
#'   1. admission_date <= culture_date <= outcome_date
#'   2. dob <= admission_date
#'   3. age in [0, 120]
#'   4. age vs DOB-computed age within 2-year tolerance
#'   5. "Died" rows must have an outcome date
#'
#' Reports violations as warnings; does not drop or modify rows.
#'
#' @param data Data frame with standard column names.
#' @param admission_col Character. Admission date column. Default "admission_date".
#' @param culture_col Character. Culture/event date column. Default "culture_date".
#' @param outcome_col Character. Outcome date column. Default "outcome_date".
#' @param dob_col Character. Date of birth column. Default "dob".
#' @param age_col Character. Age column. Default "age".
#' @param outcome_value_col Character. Outcome value column for death check.
#'   Default "final_outcome".
#'
#' @return Invisibly returns data. Prints violation summary.
#' @export
prep_validate_date_logic <- function(data,
                                     admission_col     = "admission_date",
                                     culture_col       = "culture_date",
                                     outcome_col       = "outcome_date",
                                     dob_col           = "dob",
                                     age_col           = "age",
                                     outcome_value_col = "final_outcome") {
  violations <- list()

  # Check 1: admission_date <= culture_date
  if (all(c(admission_col, culture_col) %in% names(data))) {
    bad <- which(!is.na(data[[admission_col]]) & !is.na(data[[culture_col]]) &
                   data[[admission_col]] > data[[culture_col]])
    if (length(bad) > 0L) {
      violations$adm_gt_culture <- bad
      warning(sprintf("[prep_validate_date_logic] %d row(s) have admission_date > culture_date.", length(bad)))
    }
  }

  # Check 2: culture_date <= outcome_date
  if (all(c(culture_col, outcome_col) %in% names(data))) {
    bad <- which(!is.na(data[[culture_col]]) & !is.na(data[[outcome_col]]) &
                   data[[culture_col]] > data[[outcome_col]])
    if (length(bad) > 0L) {
      violations$culture_gt_outcome <- bad
      warning(sprintf("[prep_validate_date_logic] %d row(s) have culture_date > outcome_date.", length(bad)))
    }
  }

  # Check 3: dob <= admission_date
  if (all(c(dob_col, admission_col) %in% names(data))) {
    bad <- which(!is.na(data[[dob_col]]) & !is.na(data[[admission_col]]) &
                   data[[dob_col]] > data[[admission_col]])
    if (length(bad) > 0L) {
      violations$dob_gt_admission <- bad
      warning(sprintf("[prep_validate_date_logic] %d row(s) have dob > admission_date.", length(bad)))
    }
  }

  # Check 4: age in [0, 120]
  if (age_col %in% names(data)) {
    bad <- which(!is.na(data[[age_col]]) & (data[[age_col]] < 0 | data[[age_col]] > 120))
    if (length(bad) > 0L) {
      violations$age_out_of_range <- bad
      warning(sprintf("[prep_validate_date_logic] %d row(s) have age outside [0, 120].", length(bad)))
    }
  }

  # Check 5: age vs DOB-computed age (2-year tolerance)
  if (all(c(age_col, dob_col, culture_col) %in% names(data))) {
    mask <- !is.na(data[[age_col]]) & !is.na(data[[dob_col]]) & !is.na(data[[culture_col]])
    if (any(mask)) {
      computed_age <- as.numeric(difftime(data[[culture_col]][mask],
                                          data[[dob_col]][mask], units = "days")) / 365.25
      bad_idx <- which(mask)[abs(data[[age_col]][mask] - computed_age) > 2]
      if (length(bad_idx) > 0L) {
        violations$age_dob_mismatch <- bad_idx
        warning(sprintf("[prep_validate_date_logic] %d row(s) have age vs DOB mismatch > 2 years.", length(bad_idx)))
      }
    }
  }

  # Check 6: "Died" rows must have an outcome date
  if (all(c(outcome_value_col, outcome_col) %in% names(data))) {
    bad <- which(!is.na(data[[outcome_value_col]]) &
                   data[[outcome_value_col]] == "Died" &
                   is.na(data[[outcome_col]]))
    if (length(bad) > 0L) {
      violations$died_no_outcome_date <- bad
      warning(sprintf("[prep_validate_date_logic] %d row(s) have outcome='Died' but no outcome date.", length(bad)))
    }
  }

  if (length(violations) == 0L) {
    message("[prep_validate_date_logic] All date logic checks passed.")
  } else {
    message(sprintf("[prep_validate_date_logic] %d check(s) with violations. See warnings above.",
                    length(violations)))
  }

  invisible(data)
}

