# prep_events.R
# Layer 6a: Event creation, deduplication, and readmission classification
#
# Functions:
#   - prep_create_event_ids
#   - prep_deduplicate_events  (event-aware + generic mode; merged remove_duplicate_rows)
#   - prep_flag_readmission
#   - prep_classify_readmission


#' Create Event IDs from Patient-Level Data
#'
#' Converts patient-level isolate data to event-level data by grouping repeat
#' cultures from the same infection episode.
#'
#' Event classification rules (in order of precedence):
#'
#' | Scenario                                            | Result          |
#' |-----------------------------------------------------|-----------------|
#' | Different body sites (any date)                     | Separate events |
#' | Same site, same day, same organism, same ABG        | One event       |
#' | Same site, same day, same organism, diff ABG        | Separate events |
#' | Same site, same organism, within gap_days, same ABG | Same event      |
#' | Same site, same organism, within gap_days, ABG changed | New event    |
#' | Same site, same organism, after > gap_days          | New event       |
#' | Same site, same day, different organisms            | Separate events |
#'
#' Event IDs are numbered GLOBALLY per patient in chronological order.
#'
#' @param data Data frame (long format -- one row per antibiotic test).
#' @param patient_col  Patient ID column. Default "patient_id".
#' @param date_col     Culture date column. Default "date_of_culture".
#' @param organism_col Organism column. Default "organism_normalized".
#' @param specimen_col Specimen type column. Default "specimen_type".
#' @param antibiotic_col Antibiotic name column. Default "antibiotic_name".
#' @param value_col    Susceptibility result column (S/I/R). Default "antibiotic_value".
#' @param gap_days     Days threshold: gap > gap_days triggers a new event. Default 14.
#'
#' @return Original data frame with \code{event_id} column added.
#' @export
prep_create_event_ids <- function(data,
                                   patient_col    = "patient_id",
                                   date_col       = "date_of_culture",
                                   organism_col   = "organism_normalized",
                                   specimen_col   = "specimen_type",
                                   antibiotic_col = "antibiotic_name",
                                   value_col      = "antibiotic_value",
                                   gap_days       = 14) {
  missing_cols <- setdiff(c(patient_col, date_col, organism_col), names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

  use_specimen   <- specimen_col   %in% names(data)
  use_antibiogram <- antibiotic_col %in% names(data) && value_col %in% names(data)

  if (!use_specimen)
    message(sprintf("Column '%s' not found -- events created without specimen matching.", specimen_col))
  if (!use_antibiogram)
    message("Antibiotic columns not found -- events created without antibiogram comparison.")

  message(sprintf("Creating events (gap threshold: >%d days) ...", gap_days))

  data$.pt  <- as.character(data[[patient_col]])
  data$.dt  <- as.Date(data[[date_col]])
  data$.org <- trimws(tolower(as.character(data[[organism_col]])))
  data$.spc <- if (use_specimen) trimws(tolower(as.character(data[[specimen_col]]))) else "unknown"

  data$.pt[is.na(data$.pt)]   <- "__NA__"
  data$.org[is.na(data$.org)] <- "__NA__"
  data$.spc[is.na(data$.spc)] <- "__NA__"

  if (use_antibiogram) {
    data$.abx <- trimws(as.character(data[[antibiotic_col]]))
    data$.val <- toupper(trimws(as.character(data[[value_col]])))

    df_abg <- data %>%
      dplyr::filter(.val %in% c("S", "I", "R")) %>%
      dplyr::distinct(.pt, .dt, .spc, .org, .abx, .val) %>%
      dplyr::group_by(.pt, .dt, .spc, .org) %>%
      dplyr::summarise(
        .abg_key = paste(sort(paste0(.abx, ":", .val)), collapse = "|"),
        .groups  = "drop"
      )
  }

  df_isolates <- data %>%
    dplyr::distinct(.pt, .dt, .spc, .org) %>%
    dplyr::arrange(.pt, .dt, .org, .spc)

  if (use_antibiogram) {
    df_isolates <- df_isolates %>%
      dplyr::left_join(df_abg, by = c(".pt", ".dt", ".spc", ".org"))
  } else {
    df_isolates$.abg_key <- NA_character_
  }

  df_isolates_grouped <- df_isolates %>%
    dplyr::group_by(.pt, .org, .spc) %>%
    dplyr::arrange(.dt, .by_group = TRUE) %>%
    dplyr::mutate(
      .lag_abg        = dplyr::lag(.abg_key),
      .ep_start       = dplyr::first(.dt),
      .gap_from_start = as.numeric(difftime(.dt, .ep_start, units = "days")),
      .new_ep = dplyr::case_when(
        dplyr::row_number() == 1                                         ~ TRUE,
        .gap_from_start > gap_days                                       ~ TRUE,
        !is.na(.abg_key) & !is.na(.lag_abg) & .abg_key != .lag_abg      ~ TRUE,
        TRUE                                                              ~ FALSE
      ),
      .chain_seq = cumsum(.new_ep),
      .prov_key  = paste(.pt, .org, .spc, .chain_seq, sep = "||")
    ) %>%
    dplyr::ungroup()

  event_meta <- df_isolates_grouped %>%
    dplyr::group_by(.pt, .org, .spc, .prov_key) %>%
    dplyr::summarise(
      .ep_start = {
        d <- .dt[!is.na(.dt)]
        if (length(d) == 0L) as.Date(NA_character_) else min(d)
      },
      .groups = "drop"
    ) %>%
    dplyr::arrange(.pt, .ep_start, .org, .spc)

  patient_event_index <- event_meta %>%
    dplyr::group_by(.pt) %>%
    dplyr::mutate(
      .ep_idx  = dplyr::row_number(),
      event_id = paste0(
        .pt, "_", .spc, "_",
        format(.ep_start, "%Y%m%d"), "_",
        .org, "_",
        sprintf("%03d", .ep_idx)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(.prov_key, event_id)

  df_isolates_with_id <- df_isolates_grouped %>%
    dplyr::left_join(patient_event_index, by = ".prov_key") %>%
    dplyr::select(.pt, .dt, .spc, .org, event_id)

  if ("event_id" %in% names(data)) data <- data %>% dplyr::select(-event_id)

  data <- data %>%
    dplyr::left_join(df_isolates_with_id, by = c(".pt", ".dt", ".spc", ".org"))

  n_patients  <- dplyr::n_distinct(data[[patient_col]], na.rm = TRUE)
  n_events    <- dplyr::n_distinct(data[["event_id"]], na.rm = TRUE)
  n_unmatched <- sum(is.na(data[["event_id"]]))

  message(sprintf("Done: %d patients -> %d events (%.2f per patient)",
                  n_patients, n_events, n_events / max(n_patients, 1)))
  if (n_unmatched > 0)
    message(sprintf("Warning: %d rows have no event_id (join miss).", n_unmatched))

  tmp_cols <- intersect(names(data), c(".pt", ".dt", ".spc", ".org", ".abx", ".val"))
  if (length(tmp_cols) > 0)
    data <- data %>% dplyr::select(-dplyr::all_of(tmp_cols))

  return(data)
}


#' Deduplicate Events
#'
#' Removes duplicate rows within groups. Two modes:
#'
#' \describe{
#'   \item{Event-aware (default)}{Groups by \code{event_col} + \code{organism_col}
#'     + \code{antibiotic_col} and keeps first/last antibiotic test per group.
#'     Requires \code{event_col} to be present.}
#'   \item{Generic (key_cols supplied)}{Groups by \code{key_cols} and drops
#'     duplicates. When \code{keep = "none"}, both copies of every duplicate
#'     are removed. Replaces the former \code{remove_duplicate_rows()} helper.}
#' }
#'
#' @param data Data frame.
#' @param event_col Character. Event ID column (event-aware mode). Default "event_id".
#' @param organism_col Character. Organism column (event-aware mode). Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column (event-aware mode). Default "antibiotic_normalized".
#' @param key_cols Character vector. When supplied, switches to generic mode and
#'   uses these columns for duplicate detection. NULL uses event-aware mode. Default NULL.
#' @param keep Character. "first", "last", or "none" (drop all duplicates). Default "first".
#'
#' @return Deduplicated data frame.
#' @export
prep_deduplicate_events <- function(data,
                                     event_col      = "event_id",
                                     organism_col   = "organism_normalized",
                                     antibiotic_col = "antibiotic_normalized",
                                     key_cols       = NULL,
                                     keep           = "first") {
  if (!keep %in% c("first", "last", "none", "all"))
    stop("keep must be 'first', 'last', 'none', or 'all'")

  n_before <- nrow(data)

  # --- Generic mode (key_cols supplied) ---
  if (!is.null(key_cols)) {
    missing_cols <- setdiff(key_cols, names(data))
    if (length(missing_cols) > 0)
      stop(sprintf("key_cols not found in data: %s", paste(missing_cols, collapse = ", ")))

    df_key <- data[, key_cols, drop = FALSE]
    if (keep == "none") {
      is_dup <- duplicated(df_key) | duplicated(df_key, fromLast = TRUE)
    } else {
      is_dup <- duplicated(df_key, fromLast = (keep == "last"))
    }
    n_dup <- sum(is_dup)
    if (n_dup > 0)
      message(sprintf("[prep_deduplicate_events] Found %d duplicate rows (%.1f%%).",
                      n_dup, 100 * n_dup / n_before))
    data <- data[!is_dup, ]
    message(sprintf("Deduplication: %d rows removed (%d -> %d, %.1f%% retained)",
                    n_before - nrow(data), n_before, nrow(data),
                    100 * nrow(data) / n_before))
    return(data)
  }

  # --- Event-aware mode ---
  missing_cols <- setdiff(c(event_col, organism_col, antibiotic_col), names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

  data <- data %>%
    dplyr::group_by(
      !!rlang::sym(event_col),
      !!rlang::sym(organism_col),
      !!rlang::sym(antibiotic_col)
    ) %>%
    dplyr::mutate(
      .dup_rank = dplyr::row_number(),
      .is_dup   = dplyr::n() > 1
    ) %>%
    dplyr::ungroup()

  n_duplicates <- sum(data$.is_dup & data$.dup_rank > 1)

  if (keep == "first") {
    data <- data %>% dplyr::filter(.dup_rank == 1)
  } else if (keep == "last") {
    data <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::filter(.dup_rank == max(.dup_rank)) %>%
      dplyr::ungroup()
  } else if (keep == "none") {
    data <- data %>% dplyr::filter(!.is_dup)
  }
  # keep == "all": retain everything

  data <- data %>% dplyr::select(-.dup_rank, -.is_dup)

  n_after   <- nrow(data)
  n_removed <- n_before - n_after
  message(sprintf("Deduplication: %d rows removed (%d -> %d, %.1f%% retained)",
                  n_removed, n_before, n_after, 100 * n_after / n_before))
  if (n_duplicates > 0)
    message(sprintf("  %d duplicate antibiotic tests removed", n_duplicates))

  return(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 6a)
# ---------------------------------------------------------------------------

#' Flag and Classify Readmissions
#'
#' For each patient, classifies admissions into:
#' \describe{
#'   \item{index}{First admission for this patient.}
#'   \item{linked_readmission}{Readmission within \code{gap_linked_days} -- treated
#'     as the same episode. Excluded from HAI incidence counts.}
#'   \item{new_readmission}{Readmission between \code{gap_linked_days} and
#'     \code{gap_new_days}.}
#'   \item{late_readmission}{Readmission after \code{gap_new_days} -- fully new
#'     event.}
#' }
#'
#' @param data Data frame.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param admission_col Character. Admission date column.
#'   Default "admission_date".
#' @param gap_linked_days Numeric. Gap threshold for linked readmissions.
#'   Default 30.
#' @param gap_new_days Numeric. Gap threshold for new vs late readmissions.
#'   Default 90.
#'
#' @return Data frame with \code{readmission_class} column added.
#' @export
prep_flag_readmission <- function(data,
                                   patient_col     = "patient_id",
                                   admission_col   = "admission_date",
                                   gap_linked_days = 30,
                                   gap_new_days    = 90) {
  missing_cols <- setdiff(c(patient_col, admission_col), names(data))
  if (length(missing_cols) > 0) {
    warning(sprintf("[prep_flag_readmission] Column(s) not found: %s. Skipping.",
                    paste(missing_cols, collapse = ", ")))
    data$readmission_class <- NA_character_
    return(data)
  }

  # Work on unique admission dates per patient
  data_sorted <- data %>%
    dplyr::arrange(!!rlang::sym(patient_col), !!rlang::sym(admission_col))

  data_sorted <- data_sorted %>%
    dplyr::group_by(!!rlang::sym(patient_col)) %>%
    dplyr::mutate(
      .adm_date   = as.Date(!!rlang::sym(admission_col)),
      .prev_adm   = dplyr::lag(.adm_date),
      .gap_days   = as.numeric(difftime(.adm_date, .prev_adm, units = "days")),
      readmission_class = dplyr::case_when(
        dplyr::row_number() == 1        ~ "index",
        .gap_days <= gap_linked_days    ~ "linked_readmission",
        .gap_days <= gap_new_days       ~ "new_readmission",
        TRUE                             ~ "late_readmission"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.adm_date, -.prev_adm, -.gap_days)

  dist <- table(data_sorted$readmission_class, useNA = "ifany")
  message("[prep_flag_readmission] Readmission classification:")
  print(dist)

  return(data_sorted)
}


#' Classify Readmission Type
#'
#' Applies readmission classification rules to a data frame that already has a
#' \code{readmission_class} column, re-standardizing the values.
#'
#' Useful when data arrives with non-standard labels that need mapping to the
#' controlled vocabulary used by \code{prep_flag_readmission()}.
#'
#' @param data Data frame.
#' @param readmission_col Character. Readmission column to standardize.
#'   Default "readmission_class".
#'
#' @return Data frame with standardized readmission classification.
#' @export
prep_classify_readmission <- function(data, readmission_col = "readmission_class") {
  if (!readmission_col %in% names(data)) {
    warning(sprintf("[prep_classify_readmission] Column '%s' not found. Skipping.", readmission_col))
    return(data)
  }

  val_up <- toupper(trimws(as.character(data[[readmission_col]])))

  data[[readmission_col]] <- dplyr::case_when(
    val_up %in% c("INDEX", "FIRST", "INITIAL", "PRIMARY")                    ~ "index",
    val_up %in% c("LINKED_READMISSION", "LINKED", "EARLY", "SAME_EPISODE")   ~ "linked_readmission",
    val_up %in% c("NEW_READMISSION", "NEW", "READMIT")                       ~ "new_readmission",
    val_up %in% c("LATE_READMISSION", "LATE")                                ~ "late_readmission",
    TRUE                                                                       ~ NA_character_
  )

  dist <- table(data[[readmission_col]], useNA = "ifany")
  message("[prep_classify_readmission] Readmission distribution after standardization:")
  print(dist)

  return(data)
}
