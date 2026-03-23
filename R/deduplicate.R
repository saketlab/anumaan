# deduplicate.R
# Event creation and deduplication for AMR data

#' Create Event IDs from Patient-Level Data
#'
#' Converts patient-level isolate data to event-level data by grouping repeat
#' cultures from the same infection episode.
#'
#' Event classification rules (in order of precedence):
#'
#' | Scenario                                          | Result         |
#' |---------------------------------------------------|----------------|
#' | Different body sites (any date)                   | Separate events|
#' | Same site, same day, same organism, same ABG      | One event      |
#' | Same site, same day, same organism, diff ABG      | Separate events|
#' | Same site, same organism, within gap_days, same ABG | Same event   |
#' | Same site, same organism, within gap_days, ABG changed | New event  |
#' | Same site, same organism, after > gap_days        | New event      |
#' | Same site, same day, different organisms          | Separate events|
#'
#' Event IDs are numbered GLOBALLY per patient in chronological order, so a
#' patient's events across all organisms/sites form one consistent sequence.
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
#' @return Original data frame with event_id column added.
#' @export
create_event_ids <- function(data,
                             patient_col = "patient_id",
                             date_col = "date_of_culture",
                             organism_col = "organism_normalized",
                             specimen_col = "specimen_type",
                             antibiotic_col = "antibiotic_name",
                             value_col = "antibiotic_value",
                             gap_days = 14) {
  # ---------- validate ----------
  missing_cols <- setdiff(c(patient_col, date_col, organism_col), names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  use_specimen <- specimen_col %in% names(data)
  use_antibiogram <- antibiotic_col %in% names(data) && value_col %in% names(data)

  if (!use_specimen) {
    message(sprintf("Column '%s' not found -- events created without specimen matching.", specimen_col))
  }
  if (!use_antibiogram) {
    message("Antibiotic columns not found -- events created without antibiogram comparison.")
  }

  message(sprintf("Creating events (gap threshold: >%d days) ...", gap_days))

  # ============================================================================
  # STEP 1: Standardise into internal temp columns
  # ============================================================================
  data$.pt <- as.character(data[[patient_col]])
  data$.dt <- as.Date(data[[date_col]])
  data$.org <- trimws(tolower(as.character(data[[organism_col]])))
  data$.spc <- if (use_specimen) trimws(tolower(as.character(data[[specimen_col]]))) else "unknown"

  data$.pt[is.na(data$.pt)] <- "__NA__"
  data$.org[is.na(data$.org)] <- "__NA__"
  data$.spc[is.na(data$.spc)] <- "__NA__"

  # ============================================================================
  # STEP 2: Per-isolate antibiogram fingerprint
  # Group: patient + date + specimen + organism
  # Key  : sorted "antibiotic:value" string (e.g. "amikacin:S|cipro:R")
  # Use distinct on antibiotic+value first to collapse duplicate long-format rows.
  # ============================================================================
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

  # ============================================================================
  # STEP 3: Unique isolate table -- one row per (patient, date, specimen, organism)
  # Join antibiogram key onto it.
  # ============================================================================
  df_isolates <- data %>%
    dplyr::distinct(.pt, .dt, .spc, .org) %>%
    dplyr::arrange(.pt, .dt, .org, .spc)

  if (use_antibiogram) {
    df_isolates <- df_isolates %>%
      dplyr::left_join(df_abg, by = c(".pt", ".dt", ".spc", ".org"))
  } else {
    df_isolates$.abg_key <- NA_character_
  }

  # ============================================================================
  # STEP 4: Assign provisional event keys within each patient + organism + specimen
  #
  # A NEW provisional event starts when ANY of the following is TRUE:
  #   (a) First culture in this organism+specimen chain for this patient
  #   (b) Gap from the previous culture in the chain > gap_days
  #   (c) Antibiogram changed vs the previous culture in the chain
  #
  # provisional_event_key = "patientID||organism||specimen||chain_seq"
  # This key uniquely identifies one infection episode per organism+site.
  # ============================================================================

  df_isolates_grouped <- df_isolates %>%
    dplyr::group_by(.pt, .org, .spc) %>%
    dplyr::arrange(.dt, .by_group = TRUE) %>%
    dplyr::mutate(
      .lag_abg = dplyr::lag(.abg_key),

      # initialize episode start date
      .ep_start = dplyr::first(.dt),
      .gap_from_start = as.numeric(difftime(.dt, .ep_start, units = "days")),
      .new_ep = dplyr::case_when(
        dplyr::row_number() == 1 ~ TRUE,
        .gap_from_start > gap_days ~ TRUE,
        !is.na(.abg_key) & !is.na(.lag_abg) & .abg_key != .lag_abg ~ TRUE,
        TRUE ~ FALSE
      ),
      .chain_seq = cumsum(.new_ep),
      .prov_key = paste(.pt, .org, .spc, .chain_seq, sep = "||")
    ) %>%
    dplyr::ungroup()

  # ============================================================================
  # STEP 5: Map provisional keys -> global sequential event IDs per patient
  #
  # Each unique provisional_event_key gets one entry in the event index,
  # ordered by its episode start date (then organism, specimen for tie-breaking).
  # The sequential index counts ALL events for that patient across all organisms
  # and specimen sites.
  # ============================================================================
  event_meta <- df_isolates_grouped %>%
    dplyr::group_by(.pt, .org, .spc, .prov_key) %>%
    dplyr::summarise(
      .ep_start = {
        d <- .dt[!is.na(.dt)]
        if (length(d) == 0L) as.Date(NA_character_) else min(d)
      },
      .groups = "drop"
    ) %>%
    dplyr::arrange(.pt, .ep_start, .org, .spc) # deterministic tie-breaking

  patient_event_index <- event_meta %>%
    dplyr::group_by(.pt) %>%
    dplyr::mutate(
      .ep_idx = dplyr::row_number(),
      event_id = paste0(
        .pt, "_",
        .spc, "_",
        format(.ep_start, "%Y%m%d"), "_",
        .org, "_",
        sprintf("%03d", .ep_idx)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(.prov_key, event_id)

  # ============================================================================
  # STEP 6: Join event_id back to the isolate table, then to the full long data
  # ============================================================================
  df_isolates_with_id <- df_isolates_grouped %>%
    dplyr::left_join(patient_event_index, by = ".prov_key") %>%
    dplyr::select(.pt, .dt, .spc, .org, event_id)

  # Drop pre-existing event_id to avoid .x/.y suffix collision
  if ("event_id" %in% names(data)) data <- data %>% dplyr::select(-event_id)

  data <- data %>%
    dplyr::left_join(df_isolates_with_id, by = c(".pt", ".dt", ".spc", ".org"))

  # ============================================================================
  # STEP 7: Summary
  # ============================================================================
  n_patients <- dplyr::n_distinct(data[[patient_col]], na.rm = TRUE)
  n_events <- dplyr::n_distinct(data[["event_id"]], na.rm = TRUE)
  message(sprintf(
    "Done: %d patients -> %d events (%.2f per patient)",
    n_patients, n_events, n_events / max(n_patients, 1)
  ))

  n_unmatched <- sum(is.na(data[["event_id"]]))
  if (n_unmatched > 0) {
    message(sprintf("Warning: %d rows have no event_id (join miss).", n_unmatched))
  }

  # Remove all temp columns
  tmp_cols <- names(data)[names(data) %in% c(".pt", ".dt", ".spc", ".org", ".abx", ".val")]
  if (length(tmp_cols) > 0) {
    data <- data %>% dplyr::select(-dplyr::all_of(tmp_cols))
  }

  return(data)
}


#' Deduplicate Events
#'
#' Removes duplicate isolate tests within the same event. Keeps the first
#' occurrence of each organism-antibiotic combination per event.
#'
#' @param data Data frame with event_id column
#' @param event_col Character. Event ID column. Default "event_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param keep Character. Which duplicate to keep: "first", "last", or "all".
#'   Default "first".
#'
#' @return Deduplicated data frame
#' @export
deduplicate_events <- function(data,
                               event_col = "event_id",
                               organism_col = "organism_normalized",
                               antibiotic_col = "antibiotic_normalized",
                               keep = "first") {
  missing_cols <- setdiff(c(event_col, organism_col, antibiotic_col), names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  n_before <- nrow(data)

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
  } else if (keep != "all") {
    stop("keep must be 'first', 'last', or 'all'")
  }

  data <- data %>% dplyr::select(-.dup_rank, -.is_dup)

  n_after <- nrow(data)
  n_removed <- n_before - n_after
  message(sprintf(
    "Deduplication: %d rows removed (%d -> %d, %.1f%% retained)",
    n_removed, n_before, n_after, 100 * n_after / n_before
  ))

  if (n_duplicates > 0) {
    message(sprintf("  %d duplicate antibiotic tests removed", n_duplicates))
  }

  return(data)
}
