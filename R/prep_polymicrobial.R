# prep_polymicrobial.R
# Layer 6d: Polymicrobial episode handling
#
# Functions moved and renamed from prep_ast_and_syndrome.R:
#   - flag_polymicrobial        -> prep_flag_polymicrobial
#   - compute_polymicrobial_weight -> prep_compute_poly_weights
#
# New functions: prep_split_poly_episode


#' Flag Polymicrobial Infections
#'
#' Flags polymicrobial status by counting distinct organisms per patient
#' (optionally scoped within facility and/or syndrome).
#'
#' @param data Data frame.
#' @param patient_col Character. Patient ID column. Default "patient_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param facility_col Character or NULL. Facility column for scoped counting.
#' @param facility_name Character or NULL. Filter to this facility before flagging.
#' @param syndrome_col Character or NULL. Syndrome column for scoped counting.
#' @param syndrome_name Character or NULL. Filter to this syndrome before flagging.
#'
#' @return Data frame with \code{n_organisms} and \code{is_polymicrobial} (0/1).
#' @export
prep_flag_polymicrobial <- function(data,
                                     patient_col   = "patient_id",
                                     organism_col  = "organism_normalized",
                                     facility_col  = NULL,
                                     facility_name = NULL,
                                     syndrome_col  = NULL,
                                     syndrome_name = NULL) {
  required_cols <- c(patient_col, organism_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

  if (!is.null(facility_name) && is.null(facility_col))
    stop("facility_col must be supplied when facility_name is specified.")
  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("syndrome_col must be supplied when syndrome_name is specified.")

  if (!is.null(facility_name)) {
    data <- data %>% dplyr::filter(.data[[facility_col]] == facility_name)
    message(sprintf("Filtered to facility: %s (%d rows)", facility_name, nrow(data)))
  }
  if (!is.null(syndrome_name)) {
    data <- data %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
    message(sprintf("Filtered to syndrome: %s (%d rows)", syndrome_name, nrow(data)))
  }

  message("Identifying polymicrobial infections ...")

  data <- data %>%
    dplyr::mutate(
      .temp_patient  = as.character(.data[[patient_col]]),
      .temp_organism = as.character(.data[[organism_col]])
    )

  if (!is.null(facility_col)) data$.temp_facility <- as.character(data[[facility_col]])
  if (!is.null(syndrome_col)) data$.temp_syndrome <- as.character(data[[syndrome_col]])

  data$.temp_patient[is.na(data$.temp_patient)] <- "__NA__"

  group_cols   <- c(".temp_patient",
                    if (!is.null(facility_col)) ".temp_facility",
                    if (!is.null(syndrome_col)) ".temp_syndrome")
  distinct_cols <- c(group_cols, ".temp_organism")

  poly_counts <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(distinct_cols))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n_organisms      = dplyr::n_distinct(.temp_organism),
      is_polymicrobial = as.integer(n_organisms > 1),
      .groups          = "drop"
    )

  data <- data %>% dplyr::left_join(poly_counts, by = group_cols)

  n_groups <- nrow(poly_counts)
  n_poly   <- sum(poly_counts$is_polymicrobial == 1, na.rm = TRUE)
  pct_poly <- if (n_groups > 0) 100 * n_poly / n_groups else 0

  message(sprintf("\nPolymicrobial: %d/%d patient-context groups (%.1f%%)",
                  n_poly, n_groups, pct_poly))
  message("\nOrganism count distribution per patient-context:")
  print(dplyr::arrange(dplyr::count(poly_counts, n_organisms, name = "n_groups"), n_organisms))

  data <- data %>% dplyr::select(-dplyr::starts_with(".temp_"))
  return(data)
}

# Deprecated alias — use prep_flag_polymicrobial()
#' @rdname prep_flag_polymicrobial
#' @export
flag_polymicrobial <- prep_flag_polymicrobial


#' Compute Polymicrobial Weights
#'
#' Calculates proportional weights for polymicrobial infections.
#' Monomicrobial patients always receive weight = 1.0. For polymicrobial
#' patients, weight is computed per organism within each episode using one of
#' three methods.
#'
#' @param data Data frame with \code{episode_id}, \code{is_polymicrobial} (0/1),
#'   and organism columns (output of \code{prep_flag_polymicrobial()}).
#' @param episode_col Character. Episode ID column. Default "episode_id".
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param polymicrobial_col Character. Polymicrobial flag column (0/1).
#'   Default "is_polymicrobial".
#' @param method Character. One of "monomicrobial_proportion" (default),
#'   "equal", or "manual".
#' @param weight_map Named numeric vector. Custom organism weights when
#'   \code{method = "manual"}.
#' @param facility_col Character or NULL. Scope monomicrobial reference pool
#'   per facility.
#' @param facility_name Character or NULL. Filter to this facility.
#' @param syndrome_col Character or NULL. Scope reference pool per syndrome.
#' @param syndrome_name Character or NULL. Filter to this syndrome.
#'
#' @return Data frame with \code{poly_weight}, \code{weight_method}, and
#'   \code{weight_confidence} columns.
#' @export
prep_compute_poly_weights <- function(data,
                                       episode_col       = "episode_id",
                                       organism_col      = "organism_normalized",
                                       polymicrobial_col = "is_polymicrobial",
                                       method            = "monomicrobial_proportion",
                                       weight_map        = NULL,
                                       facility_col      = NULL,
                                       facility_name     = NULL,
                                       syndrome_col      = NULL,
                                       syndrome_name     = NULL) {
  required_cols <- c(episode_col, organism_col, polymicrobial_col)
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

  method <- match.arg(method, c("monomicrobial_proportion", "equal", "manual"))
  if (method == "manual" && is.null(weight_map))
    stop("weight_map must be provided when method = 'manual'.")
  if (!is.null(facility_name) && is.null(facility_col))
    stop("facility_col must be supplied when facility_name is specified.")
  if (!is.null(syndrome_name) && is.null(syndrome_col))
    stop("syndrome_col must be supplied when syndrome_name is specified.")

  if (!is.null(facility_name)) {
    data <- data %>% dplyr::filter(.data[[facility_col]] == facility_name)
    message(sprintf("Filtered to facility: %s (%d rows)", facility_name, nrow(data)))
  }
  if (!is.null(syndrome_name)) {
    data <- data %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
    message(sprintf("Filtered to syndrome: %s (%d rows)", syndrome_name, nrow(data)))
  }

  message(sprintf("Computing polymicrobial weights using method: %s", method))

  strat_cols <- c(if (!is.null(facility_col)) facility_col,
                  if (!is.null(syndrome_col)) syndrome_col)

  data <- data %>% dplyr::mutate(poly_weight = 1.0)

  if (method == "monomicrobial_proportion") {
    mono_base <- data %>% dplyr::filter(.data[[polymicrobial_col]] == 0)

    if (length(strat_cols) > 0) {
      mono_proportions <- mono_base %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(strat_cols, organism_col)))) %>%
        dplyr::summarise(n_mono = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(strat_cols))) %>%
        dplyr::mutate(total_mono = sum(n_mono), proportion = n_mono / total_mono) %>%
        dplyr::ungroup()
    } else {
      mono_proportions <- mono_base %>%
        dplyr::count(!!rlang::sym(organism_col), name = "n_mono") %>%
        dplyr::mutate(total_mono = sum(n_mono), proportion = n_mono / total_mono)
    }

    message(sprintf("Calculated monomicrobial proportions for %d organism(s).",
                    nrow(mono_proportions)))

    join_by <- c(strat_cols, organism_col)
    data <- data %>%
      dplyr::left_join(
        mono_proportions %>% dplyr::select(dplyr::all_of(c(join_by, "proportion"))),
        by = join_by
      ) %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        poly_weight = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ 1.0,
          is.na(proportion)               ~ 1.0 / dplyr::n(),
          TRUE                             ~ proportion / sum(proportion, na.rm = TRUE)
        ),
        weight_method = dplyr::case_when(
          .data[[polymicrobial_col]] == 0 ~ "monomicrobial",
          is.na(proportion)               ~ "equal_fallback",
          TRUE                             ~ "monomicrobial_proportion"
        ),
        weight_confidence = dplyr::case_when(
          .data[[polymicrobial_col]] == 0        ~ "high",
          weight_method == "equal_fallback"       ~ "low",
          TRUE                                    ~ "high"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-proportion)

    n_fallback <- sum(data$weight_method == "equal_fallback", na.rm = TRUE)
    if (n_fallback > 0)
      message(sprintf("Warning: %d isolate(s) used equal weighting (no monomicrobial reference).",
                      n_fallback))

  } else if (method == "equal") {
    data <- data %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        poly_weight       = 1.0 / dplyr::n(),
        weight_method     = "equal",
        weight_confidence = "medium"
      ) %>%
      dplyr::ungroup()
    message("Applied equal weighting to all organisms within episodes.")

  } else if (method == "manual") {
    manual_df <- data.frame(
      org_col       = names(weight_map),
      manual_weight = as.numeric(weight_map),
      stringsAsFactors = FALSE
    )
    names(manual_df)[1L] <- organism_col

    data <- data %>%
      dplyr::left_join(manual_df, by = organism_col) %>%
      dplyr::group_by(!!rlang::sym(episode_col)) %>%
      dplyr::mutate(
        poly_weight = dplyr::case_when(
          !is.na(manual_weight) ~ manual_weight / sum(manual_weight, na.rm = TRUE),
          TRUE                  ~ 1.0 / dplyr::n()
        ),
        weight_method = dplyr::case_when(
          !is.na(manual_weight) ~ "manual",
          TRUE                  ~ "equal_fallback"
        ),
        weight_confidence = dplyr::case_when(
          weight_method == "manual" ~ "user_defined",
          TRUE                      ~ "low"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-manual_weight)

    message(sprintf("Applied manual weights for %d organisms.", length(weight_map)))
  }

  invalid_sums <- data %>%
    dplyr::group_by(!!rlang::sym(episode_col)) %>%
    dplyr::summarise(total_weight = sum(poly_weight, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(abs(total_weight - 1.0) > 0.01)

  if (nrow(invalid_sums) > 0)
    warning(sprintf("%d episode(s) have weights not summing to 1.0.", nrow(invalid_sums)))

  message("\nWeight statistics (polymicrobial isolates only):")
  print(dplyr::summarise(
    dplyr::filter(data, .data[[polymicrobial_col]] == 1),
    mean_weight   = mean(poly_weight,   na.rm = TRUE),
    median_weight = median(poly_weight, na.rm = TRUE),
    min_weight    = min(poly_weight,    na.rm = TRUE),
    max_weight    = max(poly_weight,    na.rm = TRUE)
  ))

  if (episode_col %in% names(data))
    data <- data %>% dplyr::select(-dplyr::all_of(episode_col))

  return(data)
}

# Deprecated alias — use prep_compute_poly_weights()
#' @rdname prep_compute_poly_weights
#' @export
compute_polymicrobial_weight <- prep_compute_poly_weights


# ---------------------------------------------------------------------------
# New functions (Layer 6d)
# ---------------------------------------------------------------------------

#' Split Polymicrobial Episodes by Strategy
#'
#' Applies a chosen strategy for handling polymicrobial episodes:
#' \describe{
#'   \item{fractional}{Keep all rows; \code{poly_weight} is applied downstream
#'     during DALY calculation.}
#'   \item{exclude}{Drop polymicrobial episodes from the main analysis frame.}
#'   \item{separate}{Retain polymicrobial episodes in a separate data frame
#'     stored in \code{attr(result, "poly_data")} for sensitivity analyses.}
#' }
#'
#' @param data Data frame with \code{is_polymicrobial} column (output of
#'   \code{prep_flag_polymicrobial()}).
#' @param strategy Character. "fractional" (default), "exclude", or "separate".
#' @param polymicrobial_col Character. Polymicrobial flag column (0/1).
#'   Default "is_polymicrobial".
#'
#' @return Data frame with polymicrobial rows handled per strategy. When
#'   \code{strategy = "separate"}, the poly subset is in \code{attr(result, "poly_data")}.
#' @export
prep_split_poly_episode <- function(data,
                                     strategy          = c("fractional", "exclude", "separate"),
                                     polymicrobial_col = "is_polymicrobial") {
  strategy <- match.arg(strategy)

  if (!polymicrobial_col %in% names(data)) {
    warning(sprintf("[prep_split_poly_episode] Column '%s' not found. Returning data unchanged.",
                    polymicrobial_col))
    return(data)
  }

  n_total <- nrow(data)
  n_poly  <- sum(data[[polymicrobial_col]] == 1, na.rm = TRUE)

  message(sprintf("[prep_split_poly_episode] Strategy: '%s'. Polymicrobial rows: %d / %d (%.1f%%).",
                  strategy, n_poly, n_total, 100 * n_poly / n_total))

  if (strategy == "fractional") {
    attr(data, "poly_strategy") <- "fractional"
    return(data)
  }

  is_poly <- !is.na(data[[polymicrobial_col]]) & data[[polymicrobial_col]] == 1

  if (strategy == "exclude") {
    result <- data[!is_poly, ]
    message(sprintf("[prep_split_poly_episode] Excluded %d polymicrobial rows; %d rows retained.",
                    sum(is_poly), nrow(result)))
    attr(result, "poly_strategy") <- "exclude"
    return(result)
  }

  if (strategy == "separate") {
    result  <- data[!is_poly, ]
    poly_df <- data[is_poly, ]
    message(sprintf("[prep_split_poly_episode] Main: %d rows; poly_data: %d rows.",
                    nrow(result), nrow(poly_df)))
    attr(result, "poly_strategy") <- "separate"
    attr(result, "poly_data")     <- poly_df
    return(result)
  }
}
