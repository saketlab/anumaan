# collapse.R
# Data aggregation and format conversion functions


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
#' long_data <- pivot_wide_to_long(
#'   data = raw_data,
#'   antibiotic_cols = c("AMIKACIN", "GENTAMICIN", "CIPROFLOXACIN")
#' )
#'
#' # Auto-detect columns by pattern (columns 12-53)
#' long_data <- pivot_wide_to_long(
#'   data = raw_data,
#'   antibiotic_cols = names(raw_data)[12:53]
#' )
#'
#' # Auto-detect uppercase antibiotic names
#' long_data <- pivot_wide_to_long(
#'   data = raw_data,
#'   pattern = "^[A-Z]+$"
#' )
#' }
pivot_wide_to_long <- function(data,
                               antibiotic_cols = NULL,
                               pattern = NULL,
                               id_cols = c("patient_id", "event_id", "organism_name", "date_of_culture"),
                               antibiotic_name_col = "antibiotic_name",
                               antibiotic_value_col = "antibiotic_value",
                               remove_missing = TRUE,
                               create_event_id = FALSE) {
  n_before <- nrow(data)

  # Create event_id if needed
  if (create_event_id && !"event_id" %in% names(data)) {
    message("Creating event_id column from row numbers...")
    data$event_id <- seq_len(nrow(data))
  }

  # Auto-detect antibiotic columns if not provided
  if (is.null(antibiotic_cols)) {
    if (!is.null(pattern)) {
      # Use pattern to detect
      antibiotic_cols <- names(data)[grepl(pattern, names(data))]
      message(sprintf(
        "Auto-detected %d antibiotic columns using pattern '%s'",
        length(antibiotic_cols), pattern
      ))
    } else {
      stop("Either antibiotic_cols or pattern must be provided")
    }
  }

  if (length(antibiotic_cols) == 0) {
    stop("No antibiotic columns found")
  }

  # Filter id_cols to only existing columns
  id_cols <- intersect(id_cols, names(data))

  message(sprintf(
    "Pivoting %d antibiotic columns to long format...",
    length(antibiotic_cols)
  ))

  # Check if all antibiotic columns exist
  missing_cols <- setdiff(antibiotic_cols, names(data))
  if (length(missing_cols) > 0) {
    warning(sprintf(
      "Some antibiotic columns not found: %s",
      paste(missing_cols, collapse = ", ")
    ))
    antibiotic_cols <- intersect(antibiotic_cols, names(data))
  }

  # Pivot to long format
  long_data <- data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(antibiotic_cols),
      names_to = antibiotic_name_col,
      values_to = antibiotic_value_col
    )

  n_after <- nrow(long_data)

  message(sprintf(
    "Pivoted: %d rows -> %d rows (wide -> long)",
    n_before, n_after
  ))

  # Remove missing values if requested
  if (remove_missing) {
    n_before_filter <- nrow(long_data)

    long_data <- long_data %>%
      dplyr::filter(
        !is.na(!!rlang::sym(antibiotic_value_col)),
        !!rlang::sym(antibiotic_value_col) != "",
        !!rlang::sym(antibiotic_value_col) != "-"
      )

    n_removed <- n_before_filter - nrow(long_data)

    if (n_removed > 0) {
      message(sprintf(
        "Removed %d rows with missing/empty values",
        n_removed
      ))
    }
  }

  message(sprintf(
    "Final long format: %d rows x %d columns",
    nrow(long_data), ncol(long_data)
  ))

  # Summary
  n_unique_abx <- dplyr::n_distinct(long_data[[antibiotic_name_col]])
  n_unique_events <- dplyr::n_distinct(long_data[["event_id"]], na.rm = TRUE)

  message(sprintf(
    "Contains: %d unique antibiotics across %d events",
    n_unique_abx, n_unique_events
  ))

  return(long_data)
}


#' Collapse to Antibiotic Level (OPTIONAL - Run When YOU Decide)
#'
#' **IMPORTANT**: This function removes duplicate tests. Only run this when
#' you have reviewed your data and decided to collapse duplicates.
#'
#' Aggregates multiple test results for the same organism-antibiotic combination
#' within an event. Uses "any R -> R" logic where resistance in any test
#' results in resistant classification.
#'
#' **When to use**: After you've cleaned and normalized data, if you have
#' multiple tests for the same patient-organism-antibiotic and want one
#' result per combination.
#'
#' **What it removes**: Duplicate rows based on event_id + organism + antibiotic
#'
#' @param data Data frame with susceptibility results
#' @param event_col Character. Event ID column. Default "event_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param susceptibility_col Character. Susceptibility column (S/I/R).
#'   Default "antibiotic_value".
#' @param aggregation_rule Character. Rule for aggregation: "any_R" (default,
#'   any R -> R), "most_resistant" (R > I > S), "most_common" (mode).
#'
#' @return Aggregated data frame (one row per event-organism-antibiotic)
#' @export
#'
#' @examples
#' \dontrun{
#' # Check for duplicates first
#' data %>%
#'   group_by(event_id, organism_normalized, antibiotic_normalized) %>%
#'   filter(n() > 1)
#'
#' # Then decide to collapse using any R -> R
#' collapsed <- collapse_to_antibiotic_level(data)
#'
#' # Or use most common result
#' collapsed <- collapse_to_antibiotic_level(
#'   data,
#'   aggregation_rule = "most_common"
#' )
#' }
collapse_to_antibiotic_level <- function(data,
                                         event_col = "event_id",
                                         organism_col = "organism_normalized",
                                         antibiotic_col = "antibiotic_normalized",
                                         susceptibility_col = "antibiotic_value",
                                         aggregation_rule = "any_R") {
  # Validate columns
  required_cols <- c(event_col, organism_col, antibiotic_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Validate rule
  valid_rules <- c("any_R", "most_resistant", "most_common")
  if (!aggregation_rule %in% valid_rules) {
    stop(sprintf("aggregation_rule must be one of: %s", paste(valid_rules, collapse = ", ")))
  }

  n_before <- nrow(data)

  message(sprintf(
    "Collapsing to antibiotic level using rule: %s",
    aggregation_rule
  ))

  # Identify duplicates
  duplicates <- data %>%
    dplyr::group_by(
      !!rlang::sym(event_col),
      !!rlang::sym(organism_col),
      !!rlang::sym(antibiotic_col)
    ) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  n_duplicates <- dplyr::n_distinct(
    duplicates[[event_col]],
    duplicates[[organism_col]],
    duplicates[[antibiotic_col]]
  )

  if (n_duplicates > 0) {
    message(sprintf(
      "Found %d event-organism-antibiotic combinations with multiple tests",
      n_duplicates
    ))
  }

  # Apply aggregation rule
  if (aggregation_rule == "any_R") {
    # Any R -> R
    collapsed <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::arrange(dplyr::desc(!!rlang::sym(susceptibility_col))) %>%
      dplyr::mutate(
        final_result = dplyr::case_when(
          any(!!rlang::sym(susceptibility_col) == "R", na.rm = TRUE) ~ "R",
          any(!!rlang::sym(susceptibility_col) == "I", na.rm = TRUE) ~ "I",
          any(!!rlang::sym(susceptibility_col) == "S", na.rm = TRUE) ~ "S",
          TRUE ~ NA_character_
        ),
        aggregation_method = "any_R",
        n_tests_aggregated = dplyr::n()
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  } else if (aggregation_rule == "most_resistant") {
    # R > I > S hierarchy
    collapsed <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::mutate(
        resistance_rank = dplyr::case_when(
          !!rlang::sym(susceptibility_col) == "R" ~ 3,
          !!rlang::sym(susceptibility_col) == "I" ~ 2,
          !!rlang::sym(susceptibility_col) == "S" ~ 1,
          TRUE ~ 0
        ),
        n_tests_aggregated = dplyr::n()
      ) %>%
      dplyr::arrange(dplyr::desc(resistance_rank)) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        final_result = !!rlang::sym(susceptibility_col),
        aggregation_method = "most_resistant"
      ) %>%
      dplyr::select(-resistance_rank) %>%
      dplyr::ungroup()
  } else if (aggregation_rule == "most_common") {
    # Mode (most frequent result)
    collapsed <- data %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::mutate(n_tests_aggregated = dplyr::n()) %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col),
        !!rlang::sym(susceptibility_col),
        .add = FALSE
      ) %>%
      dplyr::mutate(result_count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(
        !!rlang::sym(event_col),
        !!rlang::sym(organism_col),
        !!rlang::sym(antibiotic_col)
      ) %>%
      dplyr::arrange(dplyr::desc(result_count)) %>%
      dplyr::slice(1) %>%
      dplyr::mutate(
        final_result = !!rlang::sym(susceptibility_col),
        aggregation_method = "most_common"
      ) %>%
      dplyr::select(-result_count) %>%
      dplyr::ungroup()
  }

  # Update susceptibility column with final result
  collapsed <- collapsed %>%
    dplyr::mutate(!!susceptibility_col := final_result) %>%
    dplyr::select(-final_result)

  n_after <- nrow(collapsed)
  n_removed <- n_before - n_after

  message(sprintf(
    "Collapsed: %d -> %d rows (%d duplicates removed)",
    n_before, n_after, n_removed
  ))

  return(collapsed)
}


#' Collapse to Class Level
#'
#' Aggregates resistance at antibiotic class level instead of individual drugs.
#' Uses "any R in class -> class R" logic.
#'
#' @param data Data frame with antibiotic class information
#' @param event_col Character. Event ID column. Default "event_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param class_col Character. Antibiotic class column. Default "antibiotic_class".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param extra_cols Character vector or NULL. Additional columns to carry
#'   through the aggregation. Default NULL.
#'
#' @return Aggregated data frame (one row per event-organism-class)
#' @export
#'
#' @examples
#' \dontrun{
#' class_level <- collapse_to_class_level(data)
#' }
collapse_to_class_level <- function(data,
                                    event_col = "event_id",
                                    organism_col = "organism_normalized",
                                    class_col = "antibiotic_class",
                                    susceptibility_col = "antibiotic_value",
                                    extra_cols = NULL) {
  # ---------------------------------------------------------------------------
  # Validate required columns
  # ---------------------------------------------------------------------------
  required_cols <- c(event_col, organism_col, class_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Validate optional extra columns
  if (!is.null(extra_cols)) {
    missing_extra <- setdiff(extra_cols, names(data))
    if (length(missing_extra) > 0) {
      stop(sprintf(
        "Missing extra columns requested: %s",
        paste(missing_extra, collapse = ", ")
      ))
    }
  }

  n_before <- nrow(data)
  message("Collapsing to antibiotic class level...")

  # ---------------------------------------------------------------------------
  # Build grouping variables dynamically
  # ---------------------------------------------------------------------------
  group_vars <- c(event_col, organism_col, class_col)

  # ---------------------------------------------------------------------------
  # Core aggregation
  # ---------------------------------------------------------------------------
  collapsed <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      class_resistance = dplyr::case_when(
        any(.data[[susceptibility_col]] == "R", na.rm = TRUE) ~ "R",
        any(.data[[susceptibility_col]] == "I", na.rm = TRUE) ~ "I",
        any(.data[[susceptibility_col]] == "S", na.rm = TRUE) ~ "S",
        TRUE ~ NA_character_
      ),
      n_drugs_in_class = dplyr::n(),
      n_resistant = sum(.data[[susceptibility_col]] == "R", na.rm = TRUE),
      pct_resistant_in_class = 100 * n_resistant / n_drugs_in_class,
      drugs_tested = paste(
        sort(unique(.data[["antibiotic_normalized"]])),
        collapse = "; "
      ),
      .groups = "drop"
    )

  # ---------------------------------------------------------------------------
  # Attach optional extra columns (collapsed safely)
  # ---------------------------------------------------------------------------
  if (!is.null(extra_cols)) {
    extra_summary <- data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(extra_cols),
          ~ paste(sort(unique(.x)), collapse = "; "),
          .names = "{.col}"
        ),
        .groups = "drop"
      )

    collapsed <- collapsed %>%
      dplyr::left_join(extra_summary, by = group_vars)
  }

  # ---------------------------------------------------------------------------
  # Metadata
  # ---------------------------------------------------------------------------
  collapsed <- collapsed %>%
    dplyr::mutate(collapse_method = "class_any_R")

  n_after <- nrow(collapsed)

  message(sprintf(
    "Collapsed: %d rows -> %d class-level rows",
    n_before, n_after
  ))

  # ---------------------------------------------------------------------------
  # Class-level resistance summary (informational)
  # ---------------------------------------------------------------------------
  class_summary <- collapsed %>%
    dplyr::count(.data[[class_col]], class_resistance) %>%
    tidyr::pivot_wider(
      names_from = class_resistance,
      values_from = n,
      values_fill = 0
    )

  message("\nClass-level resistance distribution:")
  print(class_summary)

  return(collapsed)
}


#' Create Wide Format Dataset
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
#' wide_data <- create_wide_format(data)
#'
#' # Class-level wide format
#' wide_data <- create_wide_format(
#'   data,
#'   antibiotic_col = "antibiotic_class",
#'   prefix = "class_"
#' )
#' }
create_wide_format <- function(data,
                               event_col = "event_id",
                               antibiotic_col = "antibiotic_normalized",
                               susceptibility_col = "antibiotic_value",
                               prefix = "abx_",
                               keep_cols = c("patient_id", "organism_normalized", "date_of_culture")) {
  # Validate required columns
  if (!event_col %in% names(data)) {
    stop(sprintf("Column '%s' not found", event_col))
  }
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found", antibiotic_col))
  }
  if (!susceptibility_col %in% names(data)) {
    stop(sprintf("Column '%s' not found", susceptibility_col))
  }

  # Filter keep_cols to only existing columns
  keep_cols <- intersect(keep_cols, names(data))
  keep_cols <- unique(c(event_col, keep_cols))

  message(sprintf(
    "Creating wide format: pivoting '%s' column...",
    antibiotic_col
  ))

  n_before <- nrow(data)
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  # Select relevant columns
  data_subset <- data %>%
    dplyr::select(
      dplyr::all_of(keep_cols),
      !!rlang::sym(antibiotic_col),
      !!rlang::sym(susceptibility_col)
    )

  # Pivot to wide format
  wide_data <- data_subset %>%
    # Handle duplicates: take first value
    dplyr::group_by(
      !!rlang::sym(event_col),
      !!rlang::sym(antibiotic_col)
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # Pivot
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(keep_cols),
      names_from = !!rlang::sym(antibiotic_col),
      values_from = !!rlang::sym(susceptibility_col),
      names_prefix = prefix
    )

  # Clean column names (replace spaces/special chars with underscores)
  names(wide_data) <- gsub("[^A-Za-z0-9_]", "_", names(wide_data))
  names(wide_data) <- gsub("_{2,}", "_", names(wide_data)) # Remove multiple underscores

  n_after <- nrow(wide_data)
  n_antibiotics <- ncol(wide_data) - length(keep_cols)

  message(sprintf(
    "Created wide format: %d events x %d antibiotics",
    n_after, n_antibiotics
  ))

  # Check for data loss
  if (n_after != n_events_before) {
    warning(sprintf(
      "[!] Row count mismatch: %d events in original data, %d rows in wide format",
      n_events_before, n_after
    ))
  }

  return(wide_data)
}


#' Create Resistance Profile
#'
#' Generates a resistance profile string for each event summarizing
#' resistance patterns. Useful for identifying common resistance phenotypes.
#'
#' @param data Data frame with resistance data
#' @param event_col Character. Event ID column. Default "event_id".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param format Character. Output format:
#'   - "resistant_list": List resistant drugs only (default)
#'   - "full_pattern": Full S/R pattern string
#'   - "class_summary": Resistant classes only
#' @param class_col Character. Class column (required for format = "class_summary").
#'   Default "antibiotic_class".
#'
#' @return Data frame with resistance_profile column added
#' @export
#'
#' @examples
#' \dontrun{
#' # List resistant drugs
#' data_with_profile <- create_resistance_profile(data)
#'
#' # Full S/R pattern
#' data_with_profile <- create_resistance_profile(
#'   data,
#'   format = "full_pattern"
#' )
#'
#' # Resistant classes only
#' data_with_profile <- create_resistance_profile(
#'   data,
#'   format = "class_summary"
#' )
#' }
create_resistance_profile <- function(data,
                                      event_col = "event_id",
                                      antibiotic_col = "antibiotic_normalized",
                                      susceptibility_col = "antibiotic_value",
                                      format = "resistant_list",
                                      class_col = "antibiotic_class") {
  # Validate columns
  required_cols <- c(event_col, antibiotic_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (format == "class_summary" && !class_col %in% names(data)) {
    stop(sprintf("Column '%s' required for format = 'class_summary'", class_col))
  }

  # Validate format
  valid_formats <- c("resistant_list", "full_pattern", "class_summary")
  if (!format %in% valid_formats) {
    stop(sprintf("format must be one of: %s", paste(valid_formats, collapse = ", ")))
  }

  message(sprintf("Creating resistance profiles (format: %s)...", format))

  # Format A: List resistant drugs only
  if (format == "resistant_list") {
    profiles <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R") %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(sort(unique(!!rlang::sym(antibiotic_col))), collapse = "; "),
        n_resistant = dplyr::n_distinct(!!rlang::sym(antibiotic_col)),
        .groups = "drop"
      )

    # Add profile to original data
    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(
        resistance_profile = dplyr::coalesce(resistance_profile, "None"),
        n_resistant = dplyr::coalesce(n_resistant, 0L),
        profile_format = "resistant_list"
      )
  }

  # Format B: Full S/R pattern
  else if (format == "full_pattern") {
    profiles <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(!!rlang::sym(antibiotic_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(
          sprintf(
            "%s:%s",
            !!rlang::sym(antibiotic_col),
            !!rlang::sym(susceptibility_col)
          ),
          collapse = "; "
        ),
        n_resistant = sum(!!rlang::sym(susceptibility_col) == "R", na.rm = TRUE),
        n_tested = dplyr::n(),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(profile_format = "full_pattern")
  }

  # Format C: Resistant classes summary
  else if (format == "class_summary") {
    profiles <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R") %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(sort(unique(!!rlang::sym(class_col))), collapse = "; "),
        n_resistant_classes = dplyr::n_distinct(!!rlang::sym(class_col)),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(
        resistance_profile = dplyr::coalesce(resistance_profile, "None"),
        n_resistant_classes = dplyr::coalesce(n_resistant_classes, 0L),
        profile_format = "class_summary"
      )
  }

  # Summary of common profiles
  profile_freq <- data %>%
    dplyr::filter(!duplicated(!!rlang::sym(event_col))) %>%
    dplyr::count(resistance_profile, sort = TRUE) %>%
    utils::head(10)

  message("\nTop 10 resistance profiles:")
  print(profile_freq)

  return(data)
}
