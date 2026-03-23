#' Flag Polymicrobial Infections (No Specimen Type, No Episode ID)
#'
#' Flags polymicrobial status by counting distinct organisms per patient
#' (optionally within facility/syndrome scope). No specimen type and no
#' episode_id are used.
#'
#' @param data          Data frame with patient and organism columns.
#' @param patient_col   Character. Patient ID column. Default \code{"patient_id"}.
#' @param organism_col  Character. Organism column.
#'   Default \code{"organism_normalized"}.
#' @param facility_col  Character or \code{NULL}. Facility column. When supplied,
#'   polymicrobial status is scoped within each facility.
#' @param facility_name Character or \code{NULL}. When supplied with
#'   \code{facility_col}, data are first filtered to that facility.
#' @param syndrome_col  Character or \code{NULL}. Syndrome column. When supplied,
#'   polymicrobial status is scoped within each syndrome.
#' @param syndrome_name Character or \code{NULL}. When supplied with
#'   \code{syndrome_col}, data are first filtered to that syndrome.
#'
#' @return Data frame with \code{n_organisms} and \code{is_polymicrobial} (0/1).
#' @export
flag_polymicrobial <- function(data,
                               patient_col = "patient_id",
                               organism_col = "organism_normalized",
                               facility_col = NULL,
                               facility_name = NULL,
                               syndrome_col = NULL,
                               syndrome_name = NULL) {
  # -- Validate columns --------------------------------------------------------
  required_cols <- c(patient_col, organism_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be supplied when facility_name is specified.")
  }
  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }

  # -- Optional filters --------------------------------------------------------
  if (!is.null(facility_name)) {
    data <- data %>% dplyr::filter(.data[[facility_col]] == facility_name)
    message(sprintf("Filtered to facility: %s (%d rows)", facility_name, nrow(data)))
  }
  if (!is.null(syndrome_name)) {
    data <- data %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
    message(sprintf("Filtered to syndrome: %s (%d rows)", syndrome_name, nrow(data)))
  }

  message("Identifying polymicrobial infections (no specimen type, no episode_id)...")

  # -- Temp columns ------------------------------------------------------------
  data <- data %>%
    dplyr::mutate(
      .temp_patient  = as.character(.data[[patient_col]]),
      .temp_organism = as.character(.data[[organism_col]])
    )

  if (!is.null(facility_col)) data$.temp_facility <- as.character(data[[facility_col]])
  if (!is.null(syndrome_col)) data$.temp_syndrome <- as.character(data[[syndrome_col]])

  data$.temp_patient[is.na(data$.temp_patient)] <- "__NA__"

  # Group scope: patient (+ optional facility/syndrome)
  group_cols <- c(
    ".temp_patient",
    if (!is.null(facility_col)) ".temp_facility",
    if (!is.null(syndrome_col)) ".temp_syndrome"
  )

  # Distinct patient-context x organism
  distinct_cols <- c(group_cols, ".temp_organism")
  df_unique <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(distinct_cols)))

  # Count organisms per patient-context
  poly_counts <- df_unique %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_organisms      = dplyr::n_distinct(.temp_organism),
      is_polymicrobial = as.integer(n_organisms > 1),
      .groups          = "drop"
    )

  # Join back to full data
  data <- data %>%
    dplyr::left_join(poly_counts, by = group_cols)

  # Summary
  n_groups <- poly_counts %>% nrow()
  n_poly <- sum(poly_counts$is_polymicrobial == 1, na.rm = TRUE)
  pct_poly <- if (n_groups > 0) 100 * n_poly / n_groups else 0

  message(sprintf(
    "\nPolymicrobial: %d/%d patient-context groups (%.1f%%)",
    n_poly, n_groups, pct_poly
  ))
  message("\nOrganism count distribution per patient-context:")
  org_dist <- poly_counts %>%
    dplyr::count(n_organisms, name = "n_groups") %>%
    dplyr::arrange(n_organisms)
  print(org_dist)

  # Cleanup temp columns
  data <- data %>%
    dplyr::select(-dplyr::starts_with(".temp_"))

  return(data)
}


#' Compute Polymicrobial Weights
#'
#' Calculates proportional weights for polymicrobial infections.  Monomicrobial
#' patients always receive weight = 1.0.  For polymicrobial patients, weight is
#' computed per organism within each episode using one of three methods.
#'
#' When \code{facility_col} or \code{syndrome_col} are supplied, the
#' monomicrobial reference pool (method \code{"monomicrobial_proportion"}) is
#' computed within each stratum, so the reference distribution is local to each
#' facility / syndrome rather than global.
#'
#' @param data             Data frame with \code{episode_id},
#'   \code{is_polymicrobial} (0/1), and organism columns (output of
#'   \code{flag_polymicrobial()}).
#' @param episode_col      Character. Episode ID column.
#'   Default \code{"episode_id"}.
#' @param organism_col     Character. Organism column.
#'   Default \code{"organism_normalized"}.
#' @param polymicrobial_col Character. Polymicrobial flag column (0/1).
#'   Default \code{"is_polymicrobial"}.
#' @param method           Character. Weighting method:
#'   \code{"monomicrobial_proportion"} (default), \code{"equal"}, or
#'   \code{"manual"}.
#' @param weight_map       Named numeric vector. Custom organism weights when
#'   \code{method = "manual"}.
#' @param facility_col     Character or \code{NULL}. Facility column. When
#'   supplied, monomicrobial proportions are computed per-facility so each
#'   facility's local organism distribution is used as reference.
#' @param facility_name    Character or \code{NULL}. When supplied together
#'   with \code{facility_col}, data are first filtered to that facility.
#' @param syndrome_col     Character or \code{NULL}. Syndrome column. When
#'   supplied, monomicrobial proportions are computed per-syndrome.
#' @param syndrome_name    Character or \code{NULL}. When supplied together
#'   with \code{syndrome_col}, data are first filtered to that syndrome.
#'
#' @return Data frame with \code{polymicrobial_weight} column (range 0-1),
#'   plus \code{weight_method} and \code{weight_confidence} audit columns.
#'   \code{episode_id} is removed (internal use only).
#' @export
#'
#' @examples
#' \dontrun{
#' # Global monomicrobial proportion weights
#' data_weighted <- compute_polymicrobial_weight(data_flagged)
#'
#' # Per-facility reference distribution
#' data_weighted <- compute_polymicrobial_weight(
#'   data_flagged,
#'   facility_col = "center_name"
#' )
#'
#' # Filter to one facility + one syndrome, then weight
#' data_weighted <- compute_polymicrobial_weight(
#'   data_flagged,
#'   facility_col  = "center_name",
#'   facility_name = "PGIMER",
#'   syndrome_col  = "infectious_syndrome",
#'   syndrome_name = "Bloodstream infection"
#' )
#' }
compute_polymicrobial_weight <- function(data,
                                         episode_col = "episode_id",
                                         organism_col = "organism_normalized",
                                         polymicrobial_col = "is_polymicrobial",
                                         method = "monomicrobial_proportion",
                                         weight_map = NULL,
                                         facility_col = NULL,
                                         facility_name = NULL,
                                         syndrome_col = NULL,
                                         syndrome_name = NULL) {
  # -- Validate columns --------------------------------------------------------
  required_cols <- c(episode_col, organism_col, polymicrobial_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  valid_methods <- c("monomicrobial_proportion", "equal", "manual")
  if (!method %in% valid_methods) {
    stop(sprintf("method must be one of: %s", paste(valid_methods, collapse = ", ")))
  }

  if (method == "manual" && is.null(weight_map)) {
    stop("weight_map must be provided when method = 'manual'.")
  }

  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be supplied when facility_name is specified.")
  }
  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }

  # -- Optional filters --------------------------------------------------------
  if (!is.null(facility_name)) {
    data <- data %>% dplyr::filter(.data[[facility_col]] == facility_name)
    message(sprintf("Filtered to facility: %s (%d rows)", facility_name, nrow(data)))
  }
  if (!is.null(syndrome_name)) {
    data <- data %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
    message(sprintf("Filtered to syndrome: %s (%d rows)", syndrome_name, nrow(data)))
  }

  message(sprintf("Computing polymicrobial weights using method: %s", method))

  # Stratification columns for per-facility / per-syndrome reference
  strat_cols <- c(
    if (!is.null(facility_col)) facility_col,
    if (!is.null(syndrome_col)) syndrome_col
  )

  # -- Initialise weight column ------------------------------------------------
  data <- data %>% dplyr::mutate(polymicrobial_weight = 1.0)

  # -- Method A: Monomicrobial proportion -------------------------------------
  if (method == "monomicrobial_proportion") {
    # Compute monomicrobial proportions -- per stratum if strat_cols provided,
    # globally otherwise.
    mono_base <- data %>% dplyr::filter(.data[[polymicrobial_col]] == 0)

    if (length(strat_cols) > 0) {
      mono_proportions <- mono_base %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(strat_cols, organism_col)))) %>%
        dplyr::summarise(n_mono = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(strat_cols))) %>%
        dplyr::mutate(
          total_mono = sum(n_mono),
          proportion = n_mono / total_mono
        ) %>%
        dplyr::ungroup()
    } else {
      mono_proportions <- mono_base %>%
        dplyr::count(!!rlang::sym(organism_col), name = "n_mono") %>%
        dplyr::mutate(
          total_mono = sum(n_mono),
          proportion = n_mono / total_mono
        )
    }

    message(sprintf(
      "Calculated monomicrobial proportions for %d organism%s%s",
      nrow(mono_proportions),
      if (nrow(mono_proportions) != 1) "s" else "",
      if (length(strat_cols) > 0) {
        sprintf(" across %s strata", paste(strat_cols, collapse = " x "))
      } else {
        ""
      }
    ))

    # Join proportions back (by organism + strat_cols if present)
    join_by <- c(strat_cols, organism_col)
    data <- data %>%
      dplyr::left_join(
        mono_proportions %>%
          dplyr::select(dplyr::all_of(c(join_by, "proportion"))),
        by = join_by
      ) %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        polymicrobial_weight = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ 1.0,
          is.na(proportion) ~ 1.0 / dplyr::n(), # equal fallback
          TRUE ~ proportion / sum(proportion, na.rm = TRUE)
        ),
        weight_method = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ "monomicrobial",
          is.na(proportion) ~ "equal_fallback",
          TRUE ~ "monomicrobial_proportion"
        ),
        weight_confidence = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ "high",
          weight_method == "equal_fallback" ~ "low",
          TRUE ~ "high"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-proportion)

    n_fallback <- sum(data$weight_method == "equal_fallback", na.rm = TRUE)
    if (n_fallback > 0) {
      message(sprintf(
        "Warning: %d isolate(s) used equal weighting (no monomicrobial reference in stratum).",
        n_fallback
      ))
    }
  }

  # -- Method B: Equal weighting -----------------------------------------------
  else if (method == "equal") {
    data <- data %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        polymicrobial_weight = 1.0 / dplyr::n(),
        weight_method        = "equal",
        weight_confidence    = "medium"
      ) %>%
      dplyr::ungroup()

    message("Applied equal weighting to all organisms within episodes.")
  }

  # -- Method C: Manual weights ------------------------------------------------
  else if (method == "manual") {
    manual_weights_df <- data.frame(
      organism = names(weight_map),
      manual_weight = as.numeric(weight_map),
      stringsAsFactors = FALSE
    )
    names(manual_weights_df)[1L] <- organism_col

    data <- data %>%
      dplyr::left_join(manual_weights_df, by = organism_col) %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        polymicrobial_weight = dplyr::case_when(
          !is.na(manual_weight) ~ manual_weight / sum(manual_weight, na.rm = TRUE),
          TRUE ~ 1.0 / dplyr::n() # fallback for unmapped organisms
        ),
        weight_method = dplyr::case_when(
          !is.na(manual_weight) ~ "manual",
          TRUE ~ "equal_fallback"
        ),
        weight_confidence = dplyr::case_when(
          weight_method == "manual" ~ "user_defined",
          TRUE ~ "low"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-manual_weight)

    message(sprintf("Applied manual weights for %d organisms.", length(weight_map)))
  }

  # -- Validation: weights should sum to 1.0 per episode ----------------------
  weight_sums <- data %>%
    dplyr::group_by(!!rlang::sym(episode_col)) %>%
    dplyr::summarise(
      total_weight = sum(polymicrobial_weight, na.rm = TRUE),
      .groups      = "drop"
    )

  invalid_sums <- weight_sums %>%
    dplyr::filter(abs(total_weight - 1.0) > 0.01)

  if (nrow(invalid_sums) > 0) {
    warning(sprintf(
      "%d episode(s) have weights not summing to 1.0 (mixed-method fallback).",
      nrow(invalid_sums)
    ))
  }

  # -- Summary -----------------------------------------------------------------
  message("\nWeight statistics (polymicrobial isolates only):")
  weight_summary <- data %>%
    dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
    dplyr::summarise(
      mean_weight   = mean(polymicrobial_weight, na.rm = TRUE),
      median_weight = median(polymicrobial_weight, na.rm = TRUE),
      min_weight    = min(polymicrobial_weight, na.rm = TRUE),
      max_weight    = max(polymicrobial_weight, na.rm = TRUE)
    )
  print(weight_summary)

  # Remove episode_id (internal use only; avoid confusion with event_id)
  if (episode_col %in% names(data)) {
    data <- data %>% dplyr::select(-dplyr::all_of(episode_col))
  }

  return(data)
}
