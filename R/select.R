# select.R
# Resistance class selection using beta-lactam hierarchy and RR ranking

#' Select Resistance Class
#'
#' Selects a single resistance class per event using beta-lactam hierarchy
#' and relative risk (RR) values. Prevents double-counting in burden estimation
#' by choosing the most clinically relevant resistant class.
#'
#' Selection logic:
#' 1. Filter to resistant classes only (R)
#' 2. Apply beta-lactam hierarchy (Carbapenems > 4GC > 3GC > ...)
#' 3. Within same hierarchy rank, prioritize by RR value (higher RR first)
#' 4. If tied, select alphabetically for reproducibility
#'
#' @param data Data frame with resistance and RR information
#' @param event_col Character. Event ID column. Default "event_id".
#' @param class_col Character. Antibiotic class column. Default "antibiotic_class".
#' @param susceptibility_col Character. Susceptibility column. Default "antibiotic_value".
#' @param rr_col Character. RR value column. Default "rr_value".
#'   If missing, only hierarchy is used.
#' @param hierarchy Named numeric vector. Custom hierarchy (class name → rank).
#'   If NULL, uses default from get_beta_lactam_hierarchy().
#' @param filter_resistant Logical. If TRUE, only consider resistant (R) classes.
#'   Default TRUE.
#'
#' @return Data frame filtered to one resistance class per event
#' @export
#'
#' @examples
#' \dontrun{
#' # Select single resistance class per event
#' selected <- select_resistance_class(data)
#'
#' # Include susceptible classes too
#' selected <- select_resistance_class(data, filter_resistant = FALSE)
#' }
select_resistance_class <- function(data,
                                    event_col = "event_id",
                                    class_col = "antibiotic_class",
                                    susceptibility_col = "antibiotic_value",
                                    rr_col = "rr_value",
                                    hierarchy = NULL,
                                    filter_resistant = TRUE) {
  # Validate columns
  required_cols <- c(event_col, class_col, susceptibility_col)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Get hierarchy
  if (is.null(hierarchy)) {
    hierarchy <- get_beta_lactam_hierarchy()
  }

  # Check if RR column exists
  use_rr <- rr_col %in% names(data)
  if (!use_rr) {
    message(sprintf(
      "RR column '%s' not found. Using hierarchy only for selection.",
      rr_col
    ))
  }

  n_before <- nrow(data)
  n_events_before <- dplyr::n_distinct(data[[event_col]])

  message(sprintf(
    "Selecting resistance classes using hierarchy%s...",
    ifelse(use_rr, " + RR", "")
  ))

  # Filter to resistant classes if requested
  if (filter_resistant) {
    data <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R")

    message(sprintf(
      "Filtered to resistant classes: %d rows",
      nrow(data)
    ))
  }

  # Apply prioritization
  selected <- prioritize_resistance(
    data = data,
    event_col = event_col,
    class_col = class_col,
    rr_col = if (use_rr) rr_col else NULL,
    hierarchy = hierarchy
  )

  n_after <- nrow(selected)
  n_events_after <- dplyr::n_distinct(selected[[event_col]])

  message(sprintf(
    "Selected: %d rows from %d events (avg %.2f classes/event before → 1.0 after)",
    n_after,
    n_events_after,
    n_before / n_events_before
  ))

  # Show selection summary
  selection_summary <- selected %>%
    dplyr::count(!!rlang::sym(class_col), name = "n_events") %>%
    dplyr::arrange(dplyr::desc(n_events)) %>%
    utils::head(10)

  message("\nTop 10 selected classes:")
  print(selection_summary)

  return(selected)
}


#' Prioritize Resistance
#'
#' Helper function that applies hierarchy + RR ranking to select
#' the most important resistance class per event.
#'
#' @param data Data frame
#' @param event_col Character. Event ID column.
#' @param class_col Character. Class column.
#' @param rr_col Character or NULL. RR column. If NULL, uses hierarchy only.
#' @param hierarchy Named numeric vector. Hierarchy mapping.
#'
#' @return Data frame with one row per event (highest priority class)
#' @export
#'
#' @keywords internal
prioritize_resistance <- function(data,
                                  event_col,
                                  class_col,
                                  rr_col = NULL,
                                  hierarchy) {
  # Map classes to hierarchy ranks
  hierarchy_df <- data.frame(
    class = names(hierarchy),
    hierarchy_rank = as.numeric(hierarchy),
    stringsAsFactors = FALSE
  )
  names(hierarchy_df)[1] <- class_col

  data <- data %>%
    dplyr::left_join(hierarchy_df, by = class_col)

  # Assign default rank for unmapped classes (lowest priority)
  max_rank <- max(hierarchy_df$hierarchy_rank, na.rm = TRUE)
  data <- data %>%
    dplyr::mutate(
      hierarchy_rank = dplyr::coalesce(hierarchy_rank, max_rank + 1)
    )

  # Priority logic
  if (!is.null(rr_col) && rr_col %in% names(data)) {
    # Priority: hierarchy rank (lower = better) → RR (higher = better) → alphabetical
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(
        hierarchy_rank, # Lower rank = higher priority
        dplyr::desc(!!rlang::sym(rr_col)), # Higher RR = higher priority
        !!rlang::sym(class_col) # Alphabetical tie-breaker
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        selection_method = "hierarchy_rr",
        selection_confidence = "high"
      )
  } else {
    # Priority: hierarchy rank only → alphabetical
    selected <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(
        hierarchy_rank,
        !!rlang::sym(class_col)
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        selection_method = "hierarchy_only",
        selection_confidence = "medium"
      )
  }

  # Clean up temporary column
  selected <- selected %>%
    dplyr::select(-hierarchy_rank)

  # Count how many classes were dropped per event
  classes_per_event <- data %>%
    dplyr::group_by(!!rlang::sym(event_col)) %>%
    dplyr::summarise(n_classes = dplyr::n(), .groups = "drop")

  multi_class_events <- classes_per_event %>%
    dplyr::filter(n_classes > 1)

  if (nrow(multi_class_events) > 0) {
    message(sprintf(
      "Applied selection to %d events with multiple resistant classes",
      nrow(multi_class_events)
    ))

    # Show example of selection
    example_event <- multi_class_events[[event_col]][1]
    example_before <- data %>%
      dplyr::filter(!!rlang::sym(event_col) == example_event) %>%
      dplyr::select(
        !!rlang::sym(event_col),
        !!rlang::sym(class_col),
        dplyr::any_of(c(rr_col, "hierarchy_rank"))
      )

    example_after <- selected %>%
      dplyr::filter(!!rlang::sym(event_col) == example_event) %>%
      dplyr::select(
        !!rlang::sym(event_col),
        !!rlang::sym(class_col),
        selection_method
      )

    message(sprintf("\nExample selection for event '%s':", example_event))
    message("Before (all resistant classes):")
    print(example_before)
    message("After (selected class):")
    print(example_after)
  }

  return(selected)
}
