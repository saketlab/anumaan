# burden_yld.R
# YLD (Years Lived with Disability) calculation functions for AMR burden estimation
#
# Implements the GBD AMR methodology for computing deaths by infectious syndrome,
# as the foundation for YLD associated and YLD attributable to AMR.
#
# Three core population-level components:
#   D_J  : Deaths for underlying cause J
#   S_J  : Fraction of deaths related to infection for cause J
#   M_LJ : Infectious syndrome fraction (syndrome L given cause J)
#
# Synthesis:
#   D_L : Deaths by syndrome L = sum_J(D_J * S_J * M_LJ)
#          OR direct count of unique patient deaths by syndrome from facility data
#
# References:
#   Antimicrobial Resistance Collaborators. Lancet. 2022.
#   GBD 2019 Diseases and Injuries Collaborators. Lancet. 2020.


# ── D_J ───────────────────────────────────────────────────────────────────────

#' Calculate Deaths by Underlying Cause (D_J)
#'
#' Computes D_J, the number of deaths for each underlying cause J from
#' population-level vital registration or mortality data. Optionally stratified
#' by grouping variables such as age group, sex, year, or location.
#'
#' This function requires population-level data. If only facility-level data
#' are available, use \code{calculate_syndrome_deaths()} with
#' \code{facility_data} instead, which directly counts deaths by syndrome.
#'
#' @param pop_data Data frame. Population-level vital registration or mortality
#'   data. Required — this function does not accept facility-level data.
#' @param cause_col Character. Column containing the underlying cause of death
#'   (cause J), e.g. an ICD-10 code or cause name. Default \code{"cause_of_death"}.
#' @param deaths_col Character. Column with pre-aggregated death counts. Set to
#'   \code{NULL} if each row represents one individual death record.
#'   Default \code{NULL}.
#' @param groupby_cols Character vector. Additional stratification columns
#'   (e.g., \code{c("Age_bin", "gender", "year")}). Default \code{NULL}.
#'
#' @return Data frame with columns: \code{cause_col}, any \code{groupby_cols},
#'   \code{D_J} (death count), \code{D_J_method} (\code{"population"}),
#'   \code{D_J_confidence} (\code{"high"}).
#' @export
#'
#' @references
#' Antimicrobial Resistance Collaborators. Global burden of bacterial
#' antimicrobial resistance in 2019. Lancet. 2022.
#'
#' @examples
#' \dontrun{
#' # One row per death record
#' d_j <- calculate_deaths_by_cause(
#'   pop_data  = vital_reg,
#'   cause_col = "icd10_cause"
#' )
#'
#' # Pre-aggregated counts, stratified by age and sex
#' d_j <- calculate_deaths_by_cause(
#'   pop_data     = vital_reg,
#'   cause_col    = "icd10_cause",
#'   deaths_col   = "n_deaths",
#'   groupby_cols = c("Age_bin", "gender")
#' )
#' }
calculate_deaths_by_cause <- function(pop_data,
                                      cause_col = "cause_of_death",
                                      deaths_col = NULL,
                                      groupby_cols = NULL) {
  if (missing(pop_data) || is.null(pop_data)) {
    stop(
      "pop_data is required for calculate_deaths_by_cause(). ",
      "If only facility-level data are available, use calculate_syndrome_deaths() ",
      "with facility_data instead."
    )
  }

  if (!cause_col %in% names(pop_data)) {
    stop(sprintf("Cause column '%s' not found in pop_data.", cause_col))
  }
  if (!is.null(deaths_col) && !deaths_col %in% names(pop_data)) {
    stop(sprintf("Deaths column '%s' not found in pop_data.", deaths_col))
  }
  if (!is.null(groupby_cols)) {
    missing_grp <- setdiff(groupby_cols, names(pop_data))
    if (length(missing_grp) > 0) {
      stop(sprintf(
        "Grouping column(s) not found in pop_data: %s",
        paste(missing_grp, collapse = ", ")
      ))
    }
  }

  message("Computing D_J from population-level data...")

  group_vars <- unique(c(groupby_cols, cause_col))

  if (!is.null(deaths_col)) {
    result <- pop_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        D_J = sum(!!rlang::sym(deaths_col), na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    result <- pop_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(D_J = dplyr::n(), .groups = "drop")
  }

  result <- result %>%
    dplyr::mutate(
      D_J_method     = "population",
      D_J_confidence = "high"
    )

  message(sprintf(
    "D_J computed: %d total deaths across %d cause(s).",
    sum(result$D_J, na.rm = TRUE),
    dplyr::n_distinct(result[[cause_col]])
  ))

  return(result)
}


# ── S_J ───────────────────────────────────────────────────────────────────────

#' Calculate Infection Fraction of Deaths by Cause (S_J)
#'
#' Computes S_J, the fraction of deaths for underlying cause J that were
#' related to infection. Values range from 0 (no infection involvement) to 1
#' (all deaths for cause J involved infection).
#'
#' This function requires population-level data where each death record carries
#' a binary flag indicating whether infection was involved. If only facility
#' data are available, use \code{calculate_syndrome_deaths()} with
#' \code{facility_data} instead.
#'
#' @param pop_data Data frame. Population-level mortality data. Required.
#' @param cause_col Character. Underlying cause of death column (cause J).
#'   Default \code{"cause_of_death"}.
#' @param infection_flag_col Character. Binary column (TRUE/FALSE or 1/0) in
#'   \code{pop_data} indicating whether the death involved infection.
#'   Default \code{"is_infection_death"}.
#' @param groupby_cols Character vector. Stratification columns present in
#'   \code{pop_data}. Default \code{NULL}.
#'
#' @return Data frame with columns: \code{cause_col}, any \code{groupby_cols},
#'   \code{D_J} (total deaths), \code{infection_deaths} (infection-related
#'   death count), \code{S_J} (infection fraction 0-1), \code{S_J_method},
#'   \code{S_J_confidence}.
#' @export
#'
#' @references
#' Antimicrobial Resistance Collaborators. Global burden of bacterial
#' antimicrobial resistance in 2019. Lancet. 2022.
#'
#' @examples
#' \dontrun{
#' s_j <- calculate_infection_fraction(
#'   pop_data           = vital_reg,
#'   cause_col          = "icd10_cause",
#'   infection_flag_col = "is_infectious",
#'   groupby_cols       = c("Age_bin", "gender")
#' )
#' }
calculate_infection_fraction <- function(pop_data,
                                         cause_col = "cause_of_death",
                                         infection_flag_col = "is_infection_death",
                                         groupby_cols = NULL) {
  if (missing(pop_data) || is.null(pop_data)) {
    stop(
      "pop_data is required for calculate_infection_fraction(). ",
      "If only facility-level data are available, use calculate_syndrome_deaths() ",
      "with facility_data instead."
    )
  }

  required_cols <- c(cause_col, infection_flag_col)
  missing_cols <- setdiff(required_cols, names(pop_data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in pop_data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(groupby_cols)) {
    missing_grp <- setdiff(groupby_cols, names(pop_data))
    if (length(missing_grp) > 0) {
      stop(sprintf(
        "Grouping column(s) not found in pop_data: %s",
        paste(missing_grp, collapse = ", ")
      ))
    }
  }

  message("Computing S_J from population-level data...")

  group_vars <- unique(c(groupby_cols, cause_col))

  result <- pop_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      D_J = dplyr::n(),
      infection_deaths = sum(
        as.numeric(!!rlang::sym(infection_flag_col)),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      S_J            = dplyr::if_else(D_J > 0, infection_deaths / D_J, 0),
      S_J_method     = "population",
      S_J_confidence = "high"
    )

  message(sprintf(
    "S_J computed: mean infection fraction = %.3f across %d cause(s).",
    mean(result$S_J, na.rm = TRUE),
    dplyr::n_distinct(result[[cause_col]])
  ))

  return(result)
}


# ── M_LJ ──────────────────────────────────────────────────────────────────────

#' Calculate Infectious Syndrome Fraction (M_LJ)
#'
#' Computes M_LJ, the fraction of infection-related deaths for underlying cause
#' J that are attributed to infectious syndrome L. This distributes infection
#' deaths across clinical syndrome categories (e.g., bloodstream infection,
#' pneumonia, urinary tract infection).
#'
#' This function requires population-level data. If only facility data are
#' available, use \code{calculate_syndrome_deaths()} with \code{facility_data}
#' instead.
#'
#' @param pop_data Data frame. Population-level mortality data with cause and
#'   syndrome columns. Required.
#' @param cause_col Character. Underlying cause of death column (cause J).
#'   Default \code{"cause_of_death"}.
#' @param syndrome_col Character. Infectious syndrome column (syndrome L),
#'   e.g., \code{"infectious_syndrome"}. Default \code{"syndrome"}.
#' @param infection_flag_col Character. Binary column (TRUE/FALSE or 1/0)
#'   indicating infection involvement in \code{pop_data}.
#'   Default \code{"is_infection_death"}.
#' @param groupby_cols Character vector. Additional stratification columns.
#'   Default \code{NULL}.
#'
#' @return Data frame with columns: \code{cause_col}, \code{syndrome_col},
#'   any \code{groupby_cols}, \code{infection_deaths_LJ} (deaths for syndrome L
#'   given cause J), \code{infection_deaths_J} (total infection deaths for
#'   cause J), \code{M_LJ} (syndrome fraction 0-1), \code{M_LJ_method},
#'   \code{M_LJ_confidence}.
#' @export
#'
#' @references
#' Antimicrobial Resistance Collaborators. Global burden of bacterial
#' antimicrobial resistance in 2019. Lancet. 2022.
#'
#' @examples
#' \dontrun{
#' m_lj <- calculate_syndrome_fraction(
#'   pop_data           = vital_reg,
#'   cause_col          = "icd10_cause",
#'   syndrome_col       = "infectious_syndrome",
#'   infection_flag_col = "is_infectious",
#'   groupby_cols       = c("Age_bin")
#' )
#' }
calculate_syndrome_fraction <- function(pop_data,
                                        cause_col = "cause_of_death",
                                        syndrome_col = "syndrome",
                                        infection_flag_col = "is_infection_death",
                                        groupby_cols = NULL) {
  if (missing(pop_data) || is.null(pop_data)) {
    stop(
      "pop_data is required for calculate_syndrome_fraction(). ",
      "If only facility-level data are available, use calculate_syndrome_deaths() ",
      "with facility_data instead."
    )
  }

  required_cols <- c(cause_col, syndrome_col, infection_flag_col)
  missing_cols <- setdiff(required_cols, names(pop_data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Column(s) not found in pop_data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(groupby_cols)) {
    missing_grp <- setdiff(groupby_cols, names(pop_data))
    if (length(missing_grp) > 0) {
      stop(sprintf(
        "Grouping column(s) not found in pop_data: %s",
        paste(missing_grp, collapse = ", ")
      ))
    }
  }

  message("Computing M_LJ from population-level data...")

  group_vars_lj <- unique(c(groupby_cols, cause_col, syndrome_col))
  group_vars_j <- unique(c(groupby_cols, cause_col))

  # Restrict to infection-related deaths only
  infection_pop <- pop_data %>%
    dplyr::filter(as.numeric(!!rlang::sym(infection_flag_col)) == 1)

  if (nrow(infection_pop) == 0) {
    warning(sprintf(
      "No infection deaths found — all values of '%s' are 0 or FALSE.",
      infection_flag_col
    ))
  }

  # Deaths per cause J + syndrome L
  deaths_lj <- infection_pop %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_lj))) %>%
    dplyr::summarise(infection_deaths_LJ = dplyr::n(), .groups = "drop")

  # Total infection deaths per cause J (denominator for M_LJ)
  deaths_j <- deaths_lj %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_j))) %>%
    dplyr::summarise(
      infection_deaths_J = sum(infection_deaths_LJ),
      .groups = "drop"
    )

  result <- deaths_lj %>%
    dplyr::left_join(deaths_j, by = group_vars_j) %>%
    dplyr::mutate(
      M_LJ = dplyr::if_else(
        infection_deaths_J > 0,
        infection_deaths_LJ / infection_deaths_J,
        0
      ),
      M_LJ_method = "population",
      M_LJ_confidence = "high"
    )

  message(sprintf(
    "M_LJ computed: %d cause-syndrome combination(s) across %d cause(s).",
    nrow(result),
    dplyr::n_distinct(result[[cause_col]])
  ))

  return(result)
}


# ── Deaths by Syndrome (D_L) ──────────────────────────────────────────────────

#' Calculate Deaths by Infectious Syndrome (D_L)
#'
#' Computes D_L, the number of deaths due to each infectious syndrome L.
#'
#' Three operating modes, evaluated in priority order:
#' \enumerate{
#'   \item \strong{Pre-computed components} (HIGH confidence): pass the outputs
#'     of \code{calculate_deaths_by_cause()}, \code{calculate_infection_fraction()},
#'     and \code{calculate_syndrome_fraction()} via \code{d_j}, \code{s_j}, and
#'     \code{m_lj}. Computes: D_L = sum_J(D_J * S_J * M_LJ).
#'   \item \strong{Raw population data} (HIGH confidence): pass \code{pop_data}
#'     and all three components are computed internally before combining.
#'     Computes: D_L = sum_J(D_J * S_J * M_LJ).
#'   \item \strong{Facility fallback} (LOW confidence): when no population
#'     inputs are available, counts the number of \strong{unique patients} who
#'     died for each syndrome directly from \code{facility_data}.
#' }
#'
#' Use \code{syndrome} to restrict results to a single syndrome of interest.
#' Use \code{facility_col} and \code{facility_name} to restrict the facility
#' fallback to a specific site.
#'
#' @param d_j Data frame. Output of \code{calculate_deaths_by_cause()}.
#'   Supply together with \code{s_j} and \code{m_lj} for Mode 1.
#'   Default \code{NULL}.
#' @param s_j Data frame. Output of \code{calculate_infection_fraction()}.
#'   Default \code{NULL}.
#' @param m_lj Data frame. Output of \code{calculate_syndrome_fraction()}.
#'   Default \code{NULL}.
#' @param pop_data Data frame. Raw population-level data for Mode 2.
#'   Default \code{NULL}.
#' @param facility_data Data frame. Facility-level data for Mode 3 fallback.
#'   Default \code{NULL}.
#' @param cause_col Character. Underlying cause column used in Modes 1 and 2.
#'   Default \code{"cause_of_death"}.
#' @param syndrome_col Character. Column containing syndrome labels in both
#'   population and facility data. Default \code{"syndrome"}.
#' @param syndrome Character. Optional. Name of a specific syndrome to filter
#'   results to (e.g., \code{"Bloodstream infection"}). \code{NULL} returns
#'   all syndromes. Default \code{NULL}.
#' @param deaths_col Character. Pre-aggregated deaths column in \code{pop_data}
#'   (Mode 2 only). \code{NULL} means each row is one death record.
#'   Default \code{NULL}.
#' @param infection_flag_col Character. Binary infection flag column in
#'   \code{pop_data} (Modes 1 and 2). Default \code{"is_infection_death"}.
#' @param outcome_col Character. Outcome column in \code{facility_data}
#'   (Mode 3). Default \code{"final_outcome"}.
#' @param death_value Character. Value in \code{outcome_col} that represents
#'   death. Accepts any string, e.g., \code{"Died"} or \code{"Death"}.
#'   Default \code{"Died"}.
#' @param patient_col Character. Unique patient identifier column in
#'   \code{facility_data}. Deaths are counted as distinct patients, not rows.
#'   Default \code{"patient_id"}.
#' @param facility_col Character. Column containing facility names in
#'   \code{facility_data}. Required when \code{facility_name} is specified.
#'   Default \code{NULL}.
#' @param facility_name Character. Name of a specific facility to restrict the
#'   analysis to. \code{NULL} uses all facilities combined. Default \code{NULL}.
#' @param groupby_cols Character vector. Additional stratification columns
#'   (e.g., \code{c("Age_bin", "gender")}). Default \code{NULL}.
#'
#' @return Data frame with columns: \code{syndrome_col}, any \code{groupby_cols},
#'   \code{D_L} (deaths by syndrome), \code{D_L_method}, \code{D_L_confidence}.
#' @export
#'
#' @references
#' Antimicrobial Resistance Collaborators. Global burden of bacterial
#' antimicrobial resistance in 2019. Lancet. 2022.
#'
#' @examples
#' \dontrun{
#' # Mode 1: pass pre-computed components, all syndromes
#' d_j <- calculate_deaths_by_cause(pop_data = vr, cause_col = "icd10")
#' s_j <- calculate_infection_fraction(pop_data = vr, cause_col = "icd10")
#' m_lj <- calculate_syndrome_fraction(
#'   pop_data = vr, cause_col = "icd10", syndrome_col = "syndrome"
#' )
#' d_l <- calculate_syndrome_deaths(
#'   d_j = d_j, s_j = s_j, m_lj = m_lj,
#'   cause_col = "icd10", syndrome_col = "syndrome"
#' )
#'
#' # Mode 2: raw population data, filter to one syndrome
#' d_l <- calculate_syndrome_deaths(
#'   pop_data     = vital_reg,
#'   cause_col    = "icd10_cause",
#'   syndrome_col = "infectious_syndrome",
#'   syndrome     = "Bloodstream infection"
#' )
#'
#' # Mode 3: facility fallback, one facility, all syndromes
#' d_l <- calculate_syndrome_deaths(
#'   facility_data = amr_data,
#'   syndrome_col  = "specimen_normalized",
#'   outcome_col   = "final_outcome",
#'   death_value   = "Died",
#'   patient_col   = "patient_id",
#'   facility_col  = "location",
#'   facility_name = "Mumbai"
#' )
#'
#' # Mode 3: facility fallback, all facilities, one syndrome, stratified by age
#' d_l <- calculate_syndrome_deaths(
#'   facility_data = amr_data,
#'   syndrome_col  = "specimen_normalized",
#'   syndrome      = "Blood",
#'   death_value   = "Died",
#'   groupby_cols  = c("Age_bin")
#' )
#' }
calculate_syndrome_deaths <- function(d_j = NULL,
                                      s_j = NULL,
                                      m_lj = NULL,
                                      pop_data = NULL,
                                      facility_data = NULL,
                                      cause_col = "cause_of_death",
                                      syndrome_col = "syndrome",
                                      syndrome = NULL,
                                      deaths_col = NULL,
                                      infection_flag_col = "is_infection_death",
                                      outcome_col = "final_outcome",
                                      death_value = "Died",
                                      patient_col = "patient_id",
                                      facility_col = NULL,
                                      facility_name = NULL,
                                      groupby_cols = NULL) {
  has_components <- !is.null(d_j) && !is.null(s_j) && !is.null(m_lj)
  has_pop_data <- !is.null(pop_data)
  has_facility <- !is.null(facility_data)

  if (!has_components && !has_pop_data && !has_facility) {
    stop(
      "No data source provided. Supply one of:\n",
      "  (1) pre-computed d_j + s_j + m_lj\n",
      "  (2) pop_data\n",
      "  (3) facility_data"
    )
  }

  # ── Mode 1: Pre-computed components ──────────────────────────────────────────
  if (has_components) {
    message("Mode 1: Computing D_L from pre-computed D_J, S_J, M_LJ components...")

    group_vars_j <- unique(c(groupby_cols, cause_col))
    group_vars_lj <- unique(c(groupby_cols, cause_col, syndrome_col))
    group_vars_l <- unique(c(groupby_cols, syndrome_col))

    # Validate required columns in each component
    if (!cause_col %in% names(d_j)) {
      stop(sprintf("Cause column '%s' not found in d_j.", cause_col))
    }
    if (!"D_J" %in% names(d_j)) stop("Column 'D_J' not found in d_j.")
    if (!"S_J" %in% names(s_j)) stop("Column 'S_J' not found in s_j.")
    if (!"M_LJ" %in% names(m_lj)) stop("Column 'M_LJ' not found in m_lj.")
    if (!syndrome_col %in% names(m_lj)) {
      stop(sprintf("Syndrome column '%s' not found in m_lj.", syndrome_col))
    }

    # Join D_J and S_J on cause J strata, compute D_J * S_J
    dj_sj <- d_j %>%
      dplyr::select(dplyr::all_of(c(group_vars_j, "D_J"))) %>%
      dplyr::left_join(
        s_j %>% dplyr::select(dplyr::all_of(c(group_vars_j, "S_J"))),
        by = group_vars_j
      ) %>%
      dplyr::mutate(DS_J = D_J * S_J)

    # Join M_LJ, compute D_J * S_J * M_LJ per cause J + syndrome L cell
    combined <- m_lj %>%
      dplyr::select(dplyr::all_of(c(group_vars_lj, "M_LJ"))) %>%
      dplyr::left_join(dj_sj, by = group_vars_j) %>%
      dplyr::mutate(D_L_contribution = DS_J * M_LJ)

    # Sum across all causes J to get D_L per syndrome L
    result <- combined %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_l))) %>%
      dplyr::summarise(
        D_L = sum(D_L_contribution, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        D_L_method     = "population_components",
        D_L_confidence = "high"
      )

    # Syndrome filter
    if (!is.null(syndrome)) {
      n_before <- nrow(result)
      result <- result %>%
        dplyr::filter(!!rlang::sym(syndrome_col) == syndrome)
      if (nrow(result) == 0) {
        warning(sprintf(
          "Syndrome '%s' not found in results. Check values in '%s'.",
          syndrome, syndrome_col
        ))
      } else {
        message(sprintf(
          "Filtered to syndrome '%s' (%d of %d row(s) retained).",
          syndrome, nrow(result), n_before
        ))
      }
    }

    message(sprintf(
      "D_L computed: %.1f total deaths across %d syndrome(s).",
      sum(result$D_L, na.rm = TRUE),
      dplyr::n_distinct(result[[syndrome_col]])
    ))

    return(result)
  }

  # ── Mode 2: Raw population data — compute all components internally ───────────
  if (has_pop_data) {
    message("Mode 2: Computing D_L from raw population data (all components computed internally)...")

    d_j_int <- calculate_deaths_by_cause(
      pop_data     = pop_data,
      cause_col    = cause_col,
      deaths_col   = deaths_col,
      groupby_cols = groupby_cols
    )
    s_j_int <- calculate_infection_fraction(
      pop_data           = pop_data,
      cause_col          = cause_col,
      infection_flag_col = infection_flag_col,
      groupby_cols       = groupby_cols
    )
    m_lj_int <- calculate_syndrome_fraction(
      pop_data           = pop_data,
      cause_col          = cause_col,
      syndrome_col       = syndrome_col,
      infection_flag_col = infection_flag_col,
      groupby_cols       = groupby_cols
    )

    # Recurse into Mode 1 with the computed components
    return(calculate_syndrome_deaths(
      d_j          = d_j_int,
      s_j          = s_j_int,
      m_lj         = m_lj_int,
      cause_col    = cause_col,
      syndrome_col = syndrome_col,
      syndrome     = syndrome,
      groupby_cols = groupby_cols
    ))
  }

  # ── Mode 3: Facility fallback — count unique patients who died by syndrome ────
  message("Mode 3: No population data — counting deaths by syndrome from facility data (LOW confidence)...")

  # Validate required columns
  if (!syndrome_col %in% names(facility_data)) {
    stop(sprintf("Syndrome column '%s' not found in facility_data.", syndrome_col))
  }
  if (!outcome_col %in% names(facility_data)) {
    stop(sprintf("Outcome column '%s' not found in facility_data.", outcome_col))
  }
  if (!patient_col %in% names(facility_data)) {
    stop(sprintf(
      "Patient ID column '%s' not found in facility_data.", patient_col
    ))
  }

  # Validate facility filter arguments
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop(
      "facility_col must be specified when facility_name is provided. ",
      "Supply the name of the column that contains facility identifiers."
    )
  }
  if (!is.null(facility_col) && !facility_col %in% names(facility_data)) {
    stop(sprintf(
      "Facility column '%s' not found in facility_data.", facility_col
    ))
  }

  # Validate groupby_cols against facility_data
  if (!is.null(groupby_cols)) {
    missing_grp <- setdiff(groupby_cols, names(facility_data))
    if (length(missing_grp) > 0) {
      warning(sprintf(
        "Grouping column(s) not found in facility_data: %s. Ignoring.",
        paste(missing_grp, collapse = ", ")
      ))
      groupby_cols <- intersect(groupby_cols, names(facility_data))
    }
  }

  working_data <- facility_data

  # Apply facility filter
  if (!is.null(facility_col) && !is.null(facility_name)) {
    n_before <- nrow(working_data)
    working_data <- working_data %>%
      dplyr::filter(!!rlang::sym(facility_col) == facility_name)
    n_after <- nrow(working_data)
    if (n_after == 0) {
      stop(sprintf(
        "No records found for facility '%s' in column '%s'. Check facility_name.",
        facility_name, facility_col
      ))
    }
    message(sprintf(
      "Filtered to facility '%s': %d of %d record(s) retained.",
      facility_name, n_after, n_before
    ))
  } else if (!is.null(facility_col) && is.null(facility_name)) {
    n_facilities <- dplyr::n_distinct(facility_data[[facility_col]])
    message(sprintf(
      "Using all %d facilities combined. Set facility_name to restrict to one.",
      n_facilities
    ))
  }

  # Apply syndrome filter
  if (!is.null(syndrome)) {
    n_before <- nrow(working_data)
    working_data <- working_data %>%
      dplyr::filter(!!rlang::sym(syndrome_col) == syndrome)
    n_after <- nrow(working_data)
    if (n_after == 0) {
      warning(sprintf(
        "No records found for syndrome '%s' in column '%s'. Check syndrome value.",
        syndrome, syndrome_col
      ))
    } else {
      message(sprintf(
        "Filtered to syndrome '%s': %d of %d record(s) retained.",
        syndrome, n_after, n_before
      ))
    }
  }

  # Filter to deaths only
  working_data <- working_data %>%
    dplyr::filter(!!rlang::sym(outcome_col) == death_value)

  if (nrow(working_data) == 0) {
    warning(sprintf(
      "No records found where '%s' == '%s'. Check death_value.",
      outcome_col, death_value
    ))
  }

  # Count distinct patients per syndrome (and any groupby strata)
  group_vars_l <- unique(c(groupby_cols, syndrome_col))

  result <- working_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars_l))) %>%
    dplyr::summarise(
      D_L = dplyr::n_distinct(!!rlang::sym(patient_col)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      D_L_method     = "facility_direct_count",
      D_L_confidence = "low"
    )

  message(sprintf(
    "D_L (facility fallback): %d unique patient deaths across %d syndrome(s).",
    sum(result$D_L, na.rm = TRUE),
    dplyr::n_distinct(result[[syndrome_col]])
  ))

  warning(
    "D_L computed from facility data (unique patient counts). ",
    "Results reflect in-hospital deaths only and do not use the ",
    "GBD D_J * S_J * M_LJ decomposition."
  )

  return(result)
}


# ── R'Kd : Non-fatal prevalence of resistance ────────────────────────────────

#' Calculate non-fatal prevalence of resistance (R'_{k,d})
#'
#' Computes isolate-wise prevalence of resistance for each
#' pathogen–drug combination, with optional facility stratification.
#' Antibiotics are collapsed to classes, retaining the maximum
#' resistance prevalence within each class.
#'
#' @param ast_data Data frame containing AST results.
#' @param isolate_col Character. Unique isolate ID column.
#' @param pathogen_col Character. Pathogen column (K).
#' @param antibiotic_col Character. Antibiotic name column (d).
#' @param ast_result_col Character. AST result column ("R", "S", "I").
#' @param facility_col Character or NULL. Facility column if present.
#' @param antibiotic_class_map Data frame with columns:
#'   antibiotic, drug_class.
#'
#' @return Data frame with columns:
#'   pathogen, drug_class, R_kd_prime, N_tested, N_resistant
#'   (+ facility if provided)
#'
#' @export
calculate_Rkd_prime <- function(ast_data,
                                isolate_col = "isolate_id",
                                pathogen_col = "pathogen",
                                antibiotic_col = "antibiotic",
                                ast_result_col = "ast_result",
                                # Antibiotic class handling
                                drug_class_col = NULL,
                                antibiotic_class_map = NULL,
                                # Facility handling
                                facility_col = NULL,
                                facility_name = NULL,
                                # Pathogen filter
                                pathogen_name = NULL) {
  ## ── sanity checks ─────────────────────────────────────────────
  required_cols <- c(
    isolate_col, pathogen_col,
    antibiotic_col, ast_result_col
  )
  missing_cols <- setdiff(required_cols, names(ast_data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s): %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  ## ── facility filter (explicit) ────────────────────────────────
  if (!is.null(facility_name)) {
    if (is.null(facility_col)) {
      stop("facility_col must be provided if facility_name is specified.")
    }
    ast_data <- ast_data %>%
      dplyr::filter(.data[[facility_col]] == facility_name)
  }

  ## ── pathogen filter ───────────────────────────────────────────
  if (!is.null(pathogen_name)) {
    ast_data <- ast_data %>%
      dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(ast_data) == 0) {
      stop(sprintf(
        "No records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  ## ── antibiotic class resolution ───────────────────────────────
  if (!is.null(drug_class_col)) {
    if (!drug_class_col %in% names(ast_data)) {
      stop("Specified drug_class_col not found in data.")
    }

    working <- ast_data %>%
      dplyr::mutate(drug_class = .data[[drug_class_col]])
  } else {
    if (is.null(antibiotic_class_map)) {
      stop("Provide either drug_class_col or antibiotic_class_map.")
    }

    if (!all(c("antibiotic", "drug_class") %in%
      names(antibiotic_class_map))) {
      stop("antibiotic_class_map must have columns: antibiotic, drug_class")
    }

    working <- ast_data %>%
      dplyr::left_join(
        antibiotic_class_map,
        by = c(antibiotic_col = "antibiotic")
      )
  }

  ## ── resistance flag ───────────────────────────────────────────
  working <- working %>%
    dplyr::mutate(
      resistant_flag = dplyr::if_else(
        .data[[ast_result_col]] == "R", 1L, 0L
      )
    )

  ## ── isolate-wise collapse ─────────────────────────────────────
  isolate_level <- working %>%
    dplyr::group_by(
      .data[[isolate_col]],
      .data[[pathogen_col]],
      drug_class
    ) %>%
    dplyr::summarise(
      resistant = max(resistant_flag, na.rm = TRUE),
      tested    = 1L,
      .groups   = "drop"
    )

  ## ── pathogen × drug-class aggregation ─────────────────────────
  result <- isolate_level %>%
    dplyr::group_by(
      .data[[pathogen_col]],
      drug_class
    ) %>%
    dplyr::summarise(
      N_tested    = sum(tested),
      N_resistant = sum(resistant),
      R_kd_prime  = N_resistant / N_tested,
      .groups     = "drop"
    )

  return(result)
}


# ── P'LK : Non-fatal pathogen distribution ────────────────────────────────────

#' Calculate non-fatal pathogen distribution (P'_{Lk})
#'
#' Computes the non-fatal pathogen distribution for a given infectious syndrome
#' (L) using facility-level microbiology data. This quantity represents the
#' fractional contribution of each pathogen (k) to non-fatal infection cases,
#' and is used in YLD estimation.
#'
#' The unit of analysis is the \strong{patient}. Each patient contributes
#' total weight 1, distributed equally across their valid pathogens.
#' For a patient with \eqn{m_r} valid pathogens, each pathogen receives
#' weight \eqn{1/m_r}.
#'
#' For \strong{polymicrobial patients} (\code{polymicrobial_col == 1}),
#' only pathogens listed in the GLASS reference (\code{glass_ref}) for the
#' given specimen type are retained before weighting. Monomicrobial patients
#' (\code{polymicrobial_col == 0}) are never filtered.
#'
#' The pooled formula across facilities is:
#' \deqn{P'_{LK}^{\text{pooled}} =
#'   \frac{\sum_f N^{NF}_{f,L,K}}{\sum_f N^{NF}_{f,L}}}
#'
#' @param data            Data frame of facility-level microbiology records.
#' @param syndrome_col    Character. Column containing infectious syndrome labels (L).
#' @param syndrome_name   Character. Syndrome to analyse.
#' @param specimen_col    Character. Column containing specimen type.
#' @param specimen_name   Character. Specimen to restrict to (e.g., \code{"Blood"}).
#' @param polymicrobial_col Character. Column flagging polymicrobial patients
#'   (1 = polymicrobial, 0 = monomicrobial).
#' @param patient_col     Character. Unique patient identifier column.
#' @param pathogen_col    Character. Pathogen (organism) column (k).
#' @param outcome_col     Character. Final patient outcome column.
#' @param discharged_value Character. Value indicating non-fatal discharge.
#'   Default \code{"Discharged"}.
#' @param glass_ref       Character vector of valid pathogen names, or a data
#'   frame with columns \code{specimen} and \code{pathogen}. Applied to
#'   polymicrobial patients only. \code{NULL} skips GLASS filtering.
#' @param facility_col    Character or NULL. Facility identifier column. When
#'   provided without \code{facility_name}, returns both facility-level and
#'   pooled P'LK.
#' @param facility_name   Character or NULL. Restricts to a single facility.
#' @param pathogen_name   Character vector or NULL. Filter to specific pathogen(s).
#'
#' @return A list:
#'   \describe{
#'     \item{P_Lk_prime}{Pooled P'LK data frame: \code{pathogen_col},
#'       \code{N_NF_LK}, \code{N_NF_L}, \code{P_Lk_prime}.}
#'     \item{facility_level}{Per-facility P'LK (only when \code{facility_col}
#'       is supplied and \code{facility_name} is NULL).}
#'   }
#'
#' @export
calculate_P_Lk_prime <- function(data,
                                 syndrome_col,
                                 syndrome_name,
                                 specimen_col = NULL,
                                 specimen_name = NULL,
                                 polymicrobial_col,
                                 patient_col,
                                 pathogen_col,
                                 outcome_col,
                                 discharged_value = "Discharged",
                                 glass_ref = NULL,
                                 facility_col = NULL,
                                 facility_name = NULL,
                                 pathogen_name = NULL) {
  # ── Input validation ──────────────────────────────────────────────────────
  if (xor(is.null(specimen_col), is.null(specimen_name))) {
    stop("specimen_col and specimen_name must both be provided or both be NULL.")
  }

  required_cols <- c(
    syndrome_col, specimen_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  # ── Step 1: Filter syndrome + specimen (optional) + non-fatal ─────────────
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[outcome_col]] == discharged_value
    )
  if (!is.null(specimen_col)) {
    df <- df %>%
      dplyr::filter(.data[[specimen_col]] == specimen_name)
  }
  if (nrow(df) == 0) {
    stop("No non-fatal records remain after syndrome/specimen filtering.")
  }

  # ── Step 2: Optional single-facility restriction ───────────────────────────
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
  }

  # ── Step 3: Optional pathogen filter ─────────────────────────────────────
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # ── Step 4: GLASS filter — polymicrobial patients only ───────────────────
  # Monomicrobial patients (polymicrobial_col == 0) are never filtered.
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      valid_pathogens <- glass_ref %>%
        dplyr::filter(.data[["specimen"]] == specimen_name) %>%
        dplyr::pull(.data[["pathogen"]])
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
      dplyr::filter(.data[[pathogen_col]] %in% valid_pathogens)

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial patients (pathogen not on GLASS list for '%s').",
        n_removed, specimen_name
      ))
    }
    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No records remain after GLASS reference filtering.")
    }
  }

  # ── Step 5: Deduplicate to one row per patient × pathogen ────────────────
  # The raw data has one row per antibiotic tested. Collapse to patient-pathogen
  # level first so that each patient-pathogen pair contributes exactly once.
  group_cols <- if (!is.null(facility_col)) c(facility_col, patient_col) else patient_col

  patient_pathogen <- df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_cols, pathogen_col))))

  # ── Step 5b: Fractional weight 1/m_r ─────────────────────────────────────
  # m_r = distinct valid pathogens for patient r (within facility if relevant).
  patient_pathogen <- patient_pathogen %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r    = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r
    ) %>%
    dplyr::ungroup()

  # ── Step 6: N^NF_LK = weighted sum per (facility, pathogen) ──────────────
  agg_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  N_LK <- patient_pathogen %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(N_NF_LK = sum(weight, na.rm = TRUE), .groups = "drop")

  # ── Step 7: N^NF_L = unique non-fatal patients per facility ──────────────
  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  N_L <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
    dplyr::summarise(
      N_NF_L = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  # ── Step 8: Compute P'_LK and return ─────────────────────────────────────
  if (!is.null(facility_col)) {
    facility_level <- dplyr::left_join(N_LK, N_L, by = facility_col) %>%
      dplyr::mutate(P_Lk_prime = N_NF_LK / N_NF_L)

    # Pooled: sum numerators / sum denominators (not average of ratios)
    pooled <- N_LK %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(N_NF_LK = sum(N_NF_LK), .groups = "drop") %>%
      dplyr::mutate(
        N_NF_L     = sum(N_L$N_NF_L),
        P_Lk_prime = N_NF_LK / N_NF_L
      )

    message(sprintf(
      "P'LK computed: %d pathogen(s) across %d facility/facilities; pooled P'LK also returned.",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]])
    ))
    return(list(P_Lk_prime = pooled, facility_level = facility_level))
  } else {
    N_total <- N_L$N_NF_L
    result <- N_LK %>%
      dplyr::mutate(
        N_NF_L     = N_total,
        P_Lk_prime = N_NF_LK / N_NF_L
      )
    message(sprintf(
      "P'LK computed: %d pathogen(s), N^NF_L = %d patient(s).",
      nrow(result), N_total
    ))
    return(list(P_Lk_prime = result))
  }
}


# ── P_LK : Fatal pathogen distribution ────────────────────────────────────────

#' Calculate fatal pathogen distribution (P_{Lk})
#'
#' Computes the fatal pathogen distribution for a given infectious syndrome
#' (L) using facility-level microbiology data. This quantity represents the
#' fractional contribution of each pathogen (k) to \strong{fatal} infection
#' cases, and is used in YLL / mortality burden estimation.
#'
#' The unit of analysis is the \strong{patient}. Each patient contributes
#' total weight 1, distributed equally across their valid pathogens.
#' For a patient with \eqn{m_r} valid pathogens, each pathogen receives
#' weight \eqn{1/m_r}.
#'
#' For \strong{polymicrobial patients} (\code{polymicrobial_col == 1}),
#' only pathogens listed in the GLASS reference (\code{glass_ref}) for the
#' given specimen type are retained before weighting. Monomicrobial patients
#' (\code{polymicrobial_col == 0}) are never filtered.
#'
#' The pooled formula across facilities is:
#' \deqn{P_{LK}^{\text{pooled}} =
#'   \frac{\sum_f N^{F}_{f,L,K}}{\sum_f N^{F}_{f,L}}}
#'
#' @param data            Data frame of facility-level microbiology records.
#' @param syndrome_col    Character. Column containing infectious syndrome labels (L).
#' @param syndrome_name   Character. Syndrome to analyse.
#' @param specimen_col    Character or NULL. Column containing specimen type.
#'   If \code{NULL}, no specimen filter is applied.
#' @param specimen_name   Character or NULL. Specimen to restrict to (e.g.,
#'   \code{"Blood"}). Required when \code{specimen_col} is provided.
#' @param polymicrobial_col Character. Column flagging polymicrobial patients
#'   (1 = polymicrobial, 0 = monomicrobial).
#' @param patient_col     Character. Unique patient identifier column.
#' @param pathogen_col    Character. Pathogen (organism) column (k).
#' @param outcome_col     Character. Final patient outcome column.
#' @param death_value     Character. Value indicating a fatal outcome.
#'   Default \code{"Death"}.
#' @param glass_ref       Character vector of valid pathogen names, or a data
#'   frame with columns \code{specimen} and \code{pathogen}. Applied to
#'   polymicrobial patients only. When \code{specimen_name} is \code{NULL} and
#'   \code{glass_ref} is a data frame, all pathogens in the reference are used
#'   regardless of specimen. \code{NULL} skips GLASS filtering.
#' @param facility_col    Character or NULL. Facility identifier column. When
#'   provided without \code{facility_name}, returns both facility-level and
#'   pooled P_LK.
#' @param facility_name   Character or NULL. Restricts to a single facility.
#' @param pathogen_name   Character vector or NULL. Filter to specific pathogen(s).
#'
#' @return A list:
#'   \describe{
#'     \item{P_Lk_fatal}{Pooled P_LK data frame: \code{pathogen_col},
#'       \code{N_F_LK}, \code{N_F_L}, \code{P_Lk_fatal}.}
#'     \item{facility_level}{Per-facility P_LK (only when \code{facility_col}
#'       is supplied and \code{facility_name} is NULL).}
#'   }
#'
#' @export
calculate_P_Lk_fatal <- function(data,
                                 syndrome_col,
                                 syndrome_name,
                                 specimen_col = NULL,
                                 specimen_name = NULL,
                                 polymicrobial_col,
                                 patient_col,
                                 pathogen_col,
                                 outcome_col,
                                 death_value = "Death",
                                 glass_ref = NULL,
                                 facility_col = NULL,
                                 facility_name = NULL,
                                 pathogen_name = NULL) {
  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(
    syndrome_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  if (!is.null(specimen_col)) required_cols <- c(required_cols, specimen_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }
  if (!is.null(specimen_col) && is.null(specimen_name)) {
    stop("specimen_name must be provided when specimen_col is specified.")
  }

  # ── Step 1: Filter syndrome + (optional specimen) + fatal outcome ─────────
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[outcome_col]] == death_value
    )
  if (!is.null(specimen_col)) {
    df <- df %>% dplyr::filter(.data[[specimen_col]] == specimen_name)
  }
  if (nrow(df) == 0) {
    stop(sprintf(
      "No fatal records remain after syndrome%s filtering (death_value='%s').",
      if (!is.null(specimen_col)) "/specimen" else "",
      death_value
    ))
  }

  # ── Step 2: Optional single-facility restriction ───────────────────────────
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No fatal records found for facility '%s'.", facility_name))
    }
  }

  # ── Step 3: Optional pathogen filter ──────────────────────────────────────
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No fatal records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # ── Step 4: GLASS filter — polymicrobial patients only ────────────────────
  # Monomicrobial patients (polymicrobial_col == 0) are never filtered.
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      if (!is.null(specimen_name)) {
        valid_pathogens <- glass_ref %>%
          dplyr::filter(.data[["specimen"]] == specimen_name) %>%
          dplyr::pull(.data[["pathogen"]])
      } else {
        # No specimen filter: use all pathogens in the reference
        valid_pathogens <- unique(glass_ref[["pathogen"]])
      }
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
      dplyr::filter(.data[[pathogen_col]] %in% valid_pathogens)

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial fatal patients (pathogen not on GLASS list%s).",
        n_removed,
        if (!is.null(specimen_name)) sprintf(" for '%s'", specimen_name) else ""
      ))
    }
    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No fatal records remain after GLASS reference filtering.")
    }
  }

  # ── Step 5: Patient-level fractional weight 1/m_r ─────────────────────────
  # m_r = distinct valid pathogens for fatal patient r.
  group_cols <- if (!is.null(facility_col)) c(facility_col, patient_col) else patient_col

  df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r    = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r
    ) %>%
    dplyr::ungroup()

  # ── Step 6: N^F_LK = weighted sum per (facility, pathogen) ───────────────
  agg_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  N_LK <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(N_F_LK = sum(weight, na.rm = TRUE), .groups = "drop")

  # ── Step 7: N^F_L = unique fatal patients per facility ────────────────────
  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  N_L <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
    dplyr::summarise(
      N_F_L = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  # ── Step 8: Compute P_LK and return ───────────────────────────────────────
  if (!is.null(facility_col)) {
    facility_level <- dplyr::left_join(N_LK, N_L, by = facility_col) %>%
      dplyr::mutate(P_Lk_fatal = N_F_LK / N_F_L)

    # Pooled: sum numerators / sum denominators (not average of ratios)
    pooled <- N_LK %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(N_F_LK = sum(N_F_LK), .groups = "drop") %>%
      dplyr::mutate(
        N_F_L      = sum(N_L$N_F_L),
        P_Lk_fatal = N_F_LK / N_F_L
      )

    message(sprintf(
      "P_LK (fatal) computed: %d pathogen(s) across %d facility/facilities; pooled P_LK also returned.",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]])
    ))
    return(list(P_Lk_fatal = pooled, facility_level = facility_level))
  } else {
    N_total <- N_L$N_F_L
    result <- N_LK %>%
      dplyr::mutate(
        N_F_L      = N_total,
        P_Lk_fatal = N_F_LK / N_F_L
      )
    message(sprintf(
      "P_LK (fatal) computed: %d pathogen(s), N^F_L = %d fatal patient(s).",
      nrow(result), N_total
    ))
    return(list(P_Lk_fatal = result))
  }
}


# ── CFR_LK : Case fatality ratio by syndrome and pathogen ──────────────────────

#' Calculate case fatality ratio by syndrome and pathogen (CFR_{Lk})
#'
#' Computes the case fatality ratio (CFR) for each pathogen (k) within a
#' specified infectious syndrome (L) using facility-level microbiology data.
#'
#' The unit of analysis is the \strong{patient}. Each patient contributes
#' total weight 1, distributed equally across their valid pathogens via a
#' fractional weight \eqn{1/m_r} (where \eqn{m_r} is the number of distinct
#' valid pathogens for patient r). This ensures polymicrobial patients are
#' not double-counted.
#'
#' For \strong{polymicrobial patients} (\code{polymicrobial_col == 1}), only
#' pathogens listed in the GLASS reference (\code{glass_ref}) for the given
#' specimen type are retained before weighting. Monomicrobial patients
#' (\code{polymicrobial_col == 0}) are never filtered.
#'
#' When \code{facility_col} is supplied and \code{facility_name} is NULL,
#' both per-facility and pooled CFR are returned. The pooled CFR sums
#' weighted deaths and totals across facilities before dividing.
#'
#' @param data             Data frame of facility-level microbiology records.
#' @param syndrome_col     Character. Column containing infectious syndrome labels (L).
#' @param syndrome_name    Character. Syndrome to analyse.
#' @param specimen_col     Character. Column containing specimen type.
#' @param specimen_name    Character. Specimen to restrict to (e.g., \code{"Blood"}).
#' @param polymicrobial_col Character. Column flagging polymicrobial patients
#'   (1 = polymicrobial, 0 = monomicrobial).
#' @param patient_col      Character. Unique patient identifier column.
#' @param pathogen_col     Character. Pathogen (organism) column (k).
#' @param outcome_col      Character. Final patient outcome column.
#' @param death_value      Character. Value indicating death. Default \code{"Died"}.
#' @param glass_ref        Character vector of valid pathogen names, or a data
#'   frame with columns \code{specimen} and \code{pathogen}. Applied to
#'   polymicrobial patients only. \code{NULL} skips GLASS filtering.
#' @param facility_col     Character or NULL. Facility identifier column.
#' @param facility_name    Character or NULL. Restricts to a single facility.
#' @param pathogen_name    Character vector or NULL. Filter to specific pathogen(s).
#'
#' @return A list:
#'   \describe{
#'     \item{cfr_table}{Pooled CFR: \code{pathogen_col}, \code{weighted_deaths},
#'       \code{weighted_total}, \code{CFR_LK}.}
#'     \item{facility_level}{Per-facility CFR (only when \code{facility_col} is
#'       supplied and \code{facility_name} is NULL).}
#'   }
#'
#' @export

calculate_cfr_lk <- function(data,
                             syndrome_col,
                             syndrome_name,
                             specimen_col,
                             specimen_name,
                             polymicrobial_col,
                             patient_col,
                             pathogen_col,
                             outcome_col,
                             death_value = "Died",
                             glass_ref = NULL,
                             facility_col = NULL,
                             facility_name = NULL,
                             pathogen_name = NULL) {
  # ── Input validation ──────────────────────────────────────────────────────
  required_cols <- c(
    syndrome_col, specimen_col, polymicrobial_col,
    patient_col, pathogen_col, outcome_col
  )
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  # ── Step 1: Filter syndrome + specimen ───────────────────────────────────
  df <- data %>%
    dplyr::filter(
      .data[[syndrome_col]] == syndrome_name,
      .data[[specimen_col]] == specimen_name
    )
  if (nrow(df) == 0) {
    stop("No records remain after syndrome/specimen filtering.")
  }

  # ── Step 2: Optional single-facility restriction ──────────────────────────
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
  }

  # ── Step 3: Optional pathogen filter ─────────────────────────────────────
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # ── Step 4: GLASS filter — polymicrobial patients only ───────────────────
  if (!is.null(glass_ref)) {
    if (is.data.frame(glass_ref)) {
      valid_pathogens <- glass_ref %>%
        dplyr::filter(.data[["specimen"]] == specimen_name) %>%
        dplyr::pull(.data[["pathogen"]])
    } else {
      valid_pathogens <- glass_ref
    }

    df_mono <- df %>% dplyr::filter(.data[[polymicrobial_col]] == 0)
    df_poly <- df %>%
      dplyr::filter(.data[[polymicrobial_col]] == 1) %>%
      dplyr::filter(.data[[pathogen_col]] %in% valid_pathogens)

    n_removed <- sum(df[[polymicrobial_col]] == 1) - nrow(df_poly)
    if (n_removed > 0) {
      message(sprintf(
        "GLASS filter: removed %d row(s) from polymicrobial patients (pathogen not on GLASS list for '%s').",
        n_removed, specimen_name
      ))
    }
    df <- dplyr::bind_rows(df_mono, df_poly)
    if (nrow(df) == 0) {
      stop("No records remain after GLASS reference filtering.")
    }
  }

  # ── Step 5: Patient-level fractional weight 1/m_r ────────────────────────
  group_cols <- if (!is.null(facility_col)) c(facility_col, patient_col) else patient_col

  df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      m_r = dplyr::n_distinct(.data[[pathogen_col]]),
      weight = 1 / m_r,
      fatal = .data[[outcome_col]] == death_value
    ) %>%
    dplyr::ungroup()

  # ── Step 6: Weighted deaths + totals per (facility, pathogen) ─────────────
  agg_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  cfr_raw <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) %>%
    dplyr::summarise(
      weighted_deaths = sum(weight * fatal, na.rm = TRUE),
      weighted_total = sum(weight, na.rm = TRUE),
      CFR_LK = dplyr::if_else(
        weighted_total > 0,
        weighted_deaths / weighted_total,
        NA_real_
      ),
      .groups = "drop"
    )

  # ── Step 7: Pooled across facilities (if applicable) ──────────────────────
  if (!is.null(facility_col) && is.null(facility_name)) {
    pooled <- cfr_raw %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        weighted_deaths = sum(weighted_deaths),
        weighted_total = sum(weighted_total),
        CFR_LK = dplyr::if_else(
          weighted_total > 0,
          weighted_deaths / weighted_total,
          NA_real_
        ),
        .groups = "drop"
      )

    message(sprintf(
      "CFR_LK computed: %d pathogen(s) across %d facility/facilities; pooled CFR also returned.",
      dplyr::n_distinct(df[[pathogen_col]]),
      dplyr::n_distinct(df[[facility_col]])
    ))
    return(list(cfr_table = pooled, facility_level = cfr_raw))
  }

  message(sprintf(
    "CFR_LK computed: %d pathogen(s).",
    dplyr::n_distinct(df[[pathogen_col]])
  ))
  return(list(cfr_table = cfr_raw))
}


# ── Incident cases by syndrome (direct count) ─────────────────────────────────

#' Count incident cases by syndrome from facility data
#'
#' Counts the number of unique patients per infectious syndrome directly from
#' facility-level data. This is the direct-count approach to incidence —
#' no CFR, no CR_L adjustment, no pathogen weighting. Use this when you
#' want raw facility-reported case counts rather than the formula-derived
#' estimate from \code{calculate_incidence_L()}.
#'
#' Results are returned facility-wise when \code{facility_col} is supplied
#' and \code{facility_name} is NULL. When \code{facility_name} is specified,
#' only that facility is returned. When no facility information is provided,
#' a single pooled count across all records is returned.
#'
#' @param data          Data frame of facility-level records.
#' @param syndrome_col  Character. Column containing infectious syndrome labels.
#' @param syndrome_name Character. Syndrome to count cases for
#'   (e.g., \code{"Bloodstream infections"}).
#' @param patient_col   Character. Unique patient identifier column. Cases
#'   are counted as distinct patients, not rows.
#' @param facility_col  Character or NULL. Facility identifier column.
#'   When provided without \code{facility_name}, counts are broken down
#'   per facility. Default \code{NULL}.
#' @param facility_name Character or NULL. If provided, restricts the count
#'   to that facility only. Default \code{NULL}.
#'
#' @return Data frame with columns:
#'   \code{syndrome_col}, \code{n_cases} (unique patient count),
#'   and \code{facility_col} if supplied.
#' @export
count_incident_cases <- function(data,
                                 syndrome_col,
                                 syndrome_name,
                                 patient_col,
                                 facility_col = NULL,
                                 facility_name = NULL,
                                 pathogen_col = NULL,
                                 pathogen_name = NULL) {
  # ── Input validation ──────────────────────────────────────────────────────
  required_cols <- c(syndrome_col, patient_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }
  if (!is.null(facility_col) && !facility_col %in% names(data)) {
    stop(sprintf("facility_col '%s' not found in data.", facility_col))
  }
  if (!is.null(pathogen_name) && is.null(pathogen_col)) {
    stop("pathogen_col must be provided when pathogen_name is specified.")
  }
  if (!is.null(pathogen_col) && !pathogen_col %in% names(data)) {
    stop(sprintf("pathogen_col '%s' not found in data.", pathogen_col))
  }

  # ── Step 1: Filter to syndrome ────────────────────────────────────────────
  df <- data %>%
    dplyr::filter(.data[[syndrome_col]] == syndrome_name)

  if (nrow(df) == 0) {
    warning(sprintf("No records found for syndrome '%s'.", syndrome_name))
    return(data.frame())
  }

  # ── Step 1b: Optional pathogen filter ────────────────────────────────────
  if (!is.null(pathogen_name)) {
    df <- df %>% dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(df) == 0) {
      stop(sprintf(
        "No records found for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
    message(sprintf(
      "Pathogen filter applied: %s",
      paste(pathogen_name, collapse = ", ")
    ))
  }

  # ── Step 2: Optional single-facility restriction ───────────────────────────
  if (!is.null(facility_name)) {
    n_before <- nrow(df)
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
    message(sprintf(
      "Restricted to facility '%s': %d of %d record(s) retained.",
      facility_name, nrow(df), n_before
    ))
  }

  # ── Step 3: Count unique patients ─────────────────────────────────────────
  # Group by facility if facility_col is supplied (and no specific facility
  # was requested — i.e., return one row per facility).
  group_vars <- if (!is.null(facility_col)) {
    c(facility_col, syndrome_col)
  } else {
    syndrome_col
  }

  result <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      n_cases = dplyr::n_distinct(.data[[patient_col]]),
      .groups = "drop"
    )

  message(sprintf(
    "Incident cases for '%s': %d total unique patient(s)%s.",
    syndrome_name,
    sum(result$n_cases),
    if (!is.null(facility_col)) sprintf(" across %d facility/facilities", nrow(result)) else ""
  ))

  return(result)
}


# ── Top N pathogens ───────────────────────────────────────────────────────────

#' Identify top N pathogens by occurrence
#'
#' Ranks pathogens by the number of records (rows) in the dataset, optionally
#' filtered to a specific syndrome, specimen, facility, or outcome. Returns
#' the top N pathogens overall or broken down per facility.
#'
#' Use this to decide which pathogens to focus on before calling
#' \code{calculate_P_Lk_prime_BSI()}, \code{calculate_cfr_lk()}, or
#' \code{calculate_YLD()}.
#'
#' @param data          Data frame of facility-level records.
#' @param pathogen_col  Character. Pathogen column.
#' @param n             Integer. Number of top pathogens to return. Default 5.
#' @param syndrome_col  Character or NULL. Filter to this syndrome column.
#' @param syndrome_name Character or NULL. Syndrome value to filter to.
#' @param specimen_col  Character or NULL. Filter to this specimen column.
#' @param specimen_name Character or NULL. Specimen value to filter to.
#' @param outcome_col   Character or NULL. Filter to this outcome column.
#' @param outcome_name  Character or NULL. Outcome value to filter to
#'   (e.g., \code{"Died"} or \code{"Discharged"}).
#' @param facility_col  Character or NULL. Facility identifier column.
#'   When provided without \code{facility_name}, returns top N per facility.
#' @param facility_name Character or NULL. If provided, restricts to that
#'   facility only before ranking.
#'
#' @return Data frame with columns: \code{pathogen_col}, \code{n_records},
#'   \code{rank} (1 = most common), and \code{facility_col} if supplied.
#' @export
get_top_pathogens <- function(data,
                              pathogen_col,
                              n = 5L,
                              syndrome_col = NULL,
                              syndrome_name = NULL,
                              specimen_col = NULL,
                              specimen_name = NULL,
                              outcome_col = NULL,
                              outcome_name = NULL,
                              facility_col = NULL,
                              facility_name = NULL) {
  # ── Input validation ──────────────────────────────────────────────────────
  if (!pathogen_col %in% names(data)) {
    stop(sprintf("pathogen_col '%s' not found in data.", pathogen_col))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  df <- data

  # ── Optional filters ──────────────────────────────────────────────────────
  if (!is.null(syndrome_col) && !is.null(syndrome_name)) {
    df <- df %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
  }
  if (!is.null(specimen_col) && !is.null(specimen_name)) {
    df <- df %>% dplyr::filter(.data[[specimen_col]] == specimen_name)
  }
  if (!is.null(outcome_col) && !is.null(outcome_name)) {
    df <- df %>% dplyr::filter(.data[[outcome_col]] == outcome_name)
  }
  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
  }

  if (nrow(df) == 0) {
    warning("No records remain after applying filters.")
    return(data.frame())
  }

  # ── Count and rank per facility (if requested) or overall ─────────────────
  group_vars <- if (!is.null(facility_col)) {
    c(facility_col, pathogen_col)
  } else {
    pathogen_col
  }

  ranked <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(n_records = dplyr::n(), .groups = "drop")

  if (!is.null(facility_col)) {
    ranked <- ranked %>%
      dplyr::group_by(.data[[facility_col]]) %>%
      dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(n_records))) %>%
      dplyr::filter(rank <= n) %>%
      dplyr::arrange(.data[[facility_col]], rank) %>%
      dplyr::ungroup()
  } else {
    ranked <- ranked %>%
      dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(n_records))) %>%
      dplyr::filter(rank <= n) %>%
      dplyr::arrange(rank)
  }

  message(sprintf(
    "Top %d pathogen(s)%s returned.",
    n,
    if (!is.null(facility_col) && is.null(facility_name)) {
      sprintf(" per facility (%d facilities)", dplyr::n_distinct(df[[facility_col]]))
    } else {
      ""
    }
  ))

  return(ranked)
}
# ── YLD per pathogen ──────────────────────────────────────────────────────────

#' Calculate YLD per pathogen
#'
#' Computes Years Lived with Disability (YLD) attributable to each pathogen K
#' within syndrome L:
#'
#' \deqn{YLD_K = \text{Incidence}_L \times P'_{LK} \times \text{DW\_sepsis}}
#'
#' The YLD weight per incident case is drawn from a GBD-derived reference table
#' stratified by Indian state (\code{yld_ref}) when \code{DW_sepsis} is not
#' provided. If \code{DW_sepsis} is provided directly, that scalar is used for
#' all rows and \code{yld_ref} is ignored.
#'
#' \strong{Facility-level mode} (when \code{facility_col} is supplied and
#' \code{facility_name} is NULL): \code{P_Lk_prime_tbl} must contain a
#' \code{facility_col} column (i.e., the \code{facility_level} element from
#' \code{calculate_P_Lk_prime()}). \code{incidence_data} must be a data frame
#' with \code{facility_col} and an incidence count column. A
#' \code{facility_state_map} data frame linking each facility to its state is
#' required only when \code{DW_sepsis} is not supplied.
#'
#' \strong{Pooled / no-facility mode}: \code{incidence_data} is a single
#' numeric scalar and \code{P_Lk_prime_tbl} has no facility column.
#' \code{state_name} selects the YLD weight; defaults to "India" if NULL.
#'
#' @param incidence_data  Numeric scalar (pooled) or data frame with
#'   \code{facility_col} + \code{incidence_col} columns (facility-level).
#' @param P_Lk_prime_tbl  Data frame from \code{calculate_P_Lk_prime()}.
#'   Use the \code{P_Lk_prime} element for pooled mode or the
#'   \code{facility_level} element for facility-level mode.
#' @param yld_ref         Data frame with columns \code{location_name} and
#'   \code{DW_sepsis}. Loaded from
#'   \file{inst/extdata/Proxy_YLD_per_case.xlsx}. Optional if
#'   \code{DW_sepsis} is provided.
#' @param DW_sepsis Numeric scalar or NULL. If provided, this value is used
#'   directly for all rows and \code{yld_ref} is ignored.
#' @param avg_los_years Numeric scalar or NULL. Overall average length of stay
#'   in years across all patients (resistant and susceptible combined), used to
#'   convert DW to a duration-weighted value:
#'   \code{effective_DW = DW_sepsis * avg_los_years}.
#'   If NULL, DW_sepsis is used as-is (caller is responsible for duration weighting).
#' @param state_name      Character or NULL. State for YLD weight in pooled
#'   mode. NULL uses the India row. Ignored in facility-level mode when
#'   \code{DW_sepsis} is provided.
#' @param facility_col    Character or NULL. Facility identifier column.
#' @param facility_name   Character or NULL. Restrict to one facility only.
#' @param facility_state_map Data frame with \code{facility_col} and
#'   \code{state_col} columns, mapping each facility to its state.
#'   Required in facility-level mode only when \code{DW_sepsis} is not provided.
#' @param state_col       Character. Column in \code{facility_state_map}
#'   containing state names that match \code{yld_ref$location_name}.
#'   Default \code{"state"}.
#' @param pathogen_col    Character. Pathogen column in \code{P_Lk_prime_tbl}.
#'   Default \code{"pathogen"}.
#' @param pathogen_name   Character vector or NULL. If provided, restricts
#'   output to those pathogen(s) only.
#' @param plk_col         Character. P'LK column in \code{P_Lk_prime_tbl}.
#'   Default \code{"P_Lk_prime"}.
#' @param incidence_col   Character. Incidence column when
#'   \code{incidence_data} is a data frame. Default \code{"n_cases"}.
#'
#' @return Data frame with columns: \code{pathogen_col}, \code{P_Lk_prime},
#'   \code{incidence_L}, \code{DW_sepsis}, \code{YLD},
#'   and \code{facility_col} / \code{state_col} when in facility-level mode.
#' @export
calculate_YLD <- function(incidence_data,
                          P_Lk_prime_tbl,
                          yld_ref = NULL,
                          DW_sepsis = NULL,
                          avg_los_years = NULL,
                          state_name = NULL,
                          facility_col = NULL,
                          facility_name = NULL,
                          facility_state_map = NULL,
                          state_col = "state",
                          pathogen_col = "pathogen",
                          pathogen_name = NULL,
                          plk_col = "P_Lk_prime",
                          incidence_col = "n_cases") {
  # ── Input validation ──────────────────────────────────────────────────────
  use_scalar_proxy <- !is.null(DW_sepsis)

  if (use_scalar_proxy) {
    if (!is.numeric(DW_sepsis) || length(DW_sepsis) != 1 ||
      is.na(DW_sepsis)) {
      stop("DW_sepsis must be a single non-missing numeric value.")
    }
  } else {
    if (is.null(yld_ref) ||
      !all(c("location_name", "DW_sepsis") %in% names(yld_ref))) {
      stop(
        "Provide either DW_sepsis as a numeric scalar, or yld_ref ",
        "with columns: 'location_name', 'DW_sepsis'."
      )
    }
  }

  # Validate avg_los_years
  if (!is.null(avg_los_years)) {
    if (!is.numeric(avg_los_years) || length(avg_los_years) != 1 || is.na(avg_los_years) || avg_los_years <= 0) {
      stop("avg_los_years must be a single positive numeric value (mean LOS in years).")
    }
  }

  if (!pathogen_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf("pathogen_col '%s' not found in P_Lk_prime_tbl.", pathogen_col))
  }
  if (!plk_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf("plk_col '%s' not found in P_Lk_prime_tbl.", plk_col))
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  pool_by_facility <- !is.null(facility_col) && is.null(facility_name) &&
    is.data.frame(incidence_data)

  # ── Pathogen filter ───────────────────────────────────────────────────────
  if (!is.null(pathogen_name)) {
    P_Lk_prime_tbl <- P_Lk_prime_tbl %>%
      dplyr::filter(.data[[pathogen_col]] %in% pathogen_name)
    if (nrow(P_Lk_prime_tbl) == 0) {
      stop(sprintf(
        "No rows in P_Lk_prime_tbl for pathogen(s): %s",
        paste(pathogen_name, collapse = ", ")
      ))
    }
  }

  # ── Single facility restriction ───────────────────────────────────────────
  if (!is.null(facility_name) && !is.null(facility_col)) {
    if (facility_col %in% names(P_Lk_prime_tbl)) {
      P_Lk_prime_tbl <- P_Lk_prime_tbl %>%
        dplyr::filter(.data[[facility_col]] == facility_name)
    }
    if (is.data.frame(incidence_data) &&
      facility_col %in% names(incidence_data)) {
      incidence_data <- incidence_data %>%
        dplyr::filter(.data[[facility_col]] == facility_name)
    }
  }

  # ═════════════════════════════════════════════════════════════════════════
  # POOLED / NO-FACILITY MODE
  # ═════════════════════════════════════════════════════════════════════════
  if (!pool_by_facility) {
    # Resolve scalar incidence
    if (is.data.frame(incidence_data)) {
      if (!incidence_col %in% names(incidence_data)) {
        stop(sprintf(
          "incidence_col '%s' not found in incidence_data.",
          incidence_col
        ))
      }
      incidence_L <- sum(incidence_data[[incidence_col]], na.rm = TRUE)
    } else {
      incidence_L <- as.numeric(incidence_data)
    }

    # Resolve YLD weight
    if (use_scalar_proxy) {
      loc <- "user_input"
      yld_weight <- as.numeric(DW_sepsis)
    } else {
      loc <- if (is.null(state_name)) "India" else state_name
      yld_weight <- yld_ref %>%
        dplyr::filter(location_name == loc) %>%
        dplyr::pull(DW_sepsis)

      if (length(yld_weight) == 0) {
        stop(sprintf(
          "Location '%s' not found in yld_ref. Available: %s",
          loc,
          paste(head(yld_ref$location_name, 10), collapse = ", ")
        ))
      }
    }

    effective_dw <- if (!is.null(avg_los_years)) yld_weight * avg_los_years else yld_weight

    result <- P_Lk_prime_tbl %>%
      dplyr::select(dplyr::all_of(c(pathogen_col, plk_col))) %>%
      dplyr::mutate(
        incidence_L = incidence_L,
        DW_sepsis = yld_weight,
        avg_los_years = if (!is.null(avg_los_years)) avg_los_years else NA_real_,
        effective_DW = effective_dw,
        YLD = incidence_L * .data[[plk_col]] * effective_dw
      )

    message(sprintf(
      "YLD computed (pooled, location='%s'): %d pathogen(s), effective_DW=%.6f, total YLD = %.2f.",
      loc, nrow(result), effective_dw, sum(result$YLD, na.rm = TRUE)
    ))

    return(result)
  }

  # ═════════════════════════════════════════════════════════════════════════
  # FACILITY-LEVEL MODE
  # ═════════════════════════════════════════════════════════════════════════

  if (!use_scalar_proxy) {
    if (is.null(facility_state_map)) {
      stop("facility_state_map is required in facility-level mode.")
    }
    if (!all(c(facility_col, state_col) %in% names(facility_state_map))) {
      stop(sprintf(
        "facility_state_map must have columns '%s' and '%s'.",
        facility_col, state_col
      ))
    }
  }

  if (!facility_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf(
      "facility_col '%s' not found in P_Lk_prime_tbl. ",
      facility_col,
      "Pass the 'facility_level' element from calculate_P_Lk_prime()."
    ))
  }
  if (!incidence_col %in% names(incidence_data)) {
    stop(sprintf(
      "incidence_col '%s' not found in incidence_data.",
      incidence_col
    ))
  }

  # Join P'LK with incidence
  result <- P_Lk_prime_tbl %>%
    dplyr::select(dplyr::all_of(c(facility_col, pathogen_col, plk_col))) %>%
    dplyr::left_join(
      incidence_data %>%
        dplyr::select(dplyr::all_of(c(facility_col, incidence_col))),
      by = facility_col
    )

  if (use_scalar_proxy) {
    result <- result %>%
      dplyr::mutate(DW_sepsis = as.numeric(.env$DW_sepsis))
  } else {
    result <- result %>%
      dplyr::left_join(
        facility_state_map %>%
          dplyr::select(dplyr::all_of(c(facility_col, state_col))) %>%
          dplyr::distinct(),
        by = facility_col
      ) %>%
      dplyr::left_join(
        yld_ref %>%
          dplyr::select(location_name, DW_sepsis) %>%
          dplyr::rename(!!state_col := location_name),
        by = state_col
      )

    # Warn if any facilities couldn't be matched to a state YLD weight
    missing_yld <- result %>%
      dplyr::filter(is.na(DW_sepsis)) %>%
      dplyr::pull(.data[[facility_col]]) %>%
      unique()

    if (length(missing_yld) > 0) {
      warning(sprintf(
        "No YLD weight found for facility/state: %s. Using India-wide fallback.",
        paste(missing_yld, collapse = ", ")
      ))
      india_yld <- yld_ref %>%
        dplyr::filter(location_name == "India") %>%
        dplyr::pull(DW_sepsis)

      result <- result %>%
        dplyr::mutate(
          DW_sepsis = dplyr::if_else(
            is.na(DW_sepsis), india_yld, DW_sepsis
          )
        )
    }
  }

  result <- result %>%
    dplyr::mutate(
      avg_los_years = if (!is.null(avg_los_years)) avg_los_years else NA_real_,
      effective_DW = if (!is.null(avg_los_years)) DW_sepsis * avg_los_years else DW_sepsis,
      YLD = .data[[incidence_col]] * .data[[plk_col]] * effective_DW
    )

  message(sprintf(
    "YLD computed (facility-level): %d facility/facilities, %d pathogen(s), total YLD = %.2f.",
    dplyr::n_distinct(result[[facility_col]]),
    dplyr::n_distinct(result[[pathogen_col]]),
    sum(result$YLD, na.rm = TRUE)
  ))

  return(result)
}


# ── CR_L : CFR adjustment factor ───────────────────────────────────────────────

#' Calculate the CFR adjustment factor (CR_L)
#'
#' Computes CR_L, the factor that adjusts a hospital-derived CFR to account for
#' infection cases managed outside the inpatient setting. The adjustment type
#' for each syndrome is looked up from the \code{adjustment_ref} table
#' (loaded from \file{inst/extdata/adjustment_for_CFR}).
#'
#' Three adjustment types are supported:
#' \describe{
#'   \item{None}{CR_L = 1. The hospital CFR applies directly (e.g., BSI,
#'     Meningitis, hospital-acquired infections).}
#'   \item{Inpatient ratio}{
#'     CR_L = (patients with \eqn{\ge} 1 inpatient visit) / (all patients).
#'     Used when community cases are captured partly in outpatient data
#'     (e.g., community-acquired LRI, UTI).}
#'   \item{Outpatient to inpatient ratio}{
#'     CR_L = (patients with \eqn{\ge} 1 outpatient AND \eqn{\ge} 1 inpatient
#'     visit) / (patients with \eqn{\ge} 1 outpatient visit).
#'     Used for syndromes where OP-to-IP transition captures disease severity
#'     (e.g., STI, Skin, Eye, Oral, Bone/joint infections).}
#' }
#'
#' @param data             Data frame of facility-level records.
#' @param syndrome_col     Character. Column containing infectious syndrome labels.
#' @param syndrome_name    Character. Syndrome to compute CR_L for.
#' @param patient_col      Character. Unique patient identifier column.
#' @param visit_type_col   Character. Column indicating visit type per record
#'   (inpatient / outpatient).
#' @param inpatient_value  Character. Value in \code{visit_type_col} that denotes
#'   an inpatient visit. Default \code{"Inpatient"}.
#' @param outpatient_value Character. Value in \code{visit_type_col} that denotes
#'   an outpatient visit. Default \code{"Outpatient"}.
#' @param adjustment_ref   Data frame with columns \code{infectious_syndrome} and
#'   \code{adjustment_factor_on_CFR}. Load from
#'   \file{inst/extdata/adjustment_for_CFR}.
#' @param facility_col     Character or NULL. Facility identifier column. When
#'   provided (and \code{facility_name} is NULL), CR_L is returned per facility.
#' @param facility_name    Character or NULL. If provided, restricts to that
#'   facility only.
#'
#' @return Data frame with columns \code{syndrome} (= \code{syndrome_name}),
#'   \code{adjustment_type}, \code{CR_L}, and (when \code{facility_col} is
#'   supplied) \code{facility_col}. Additional columns
#'   (\code{n_inpatient}, \code{n_total}, etc.) give the raw counts used.
#' @export

calculate_CR_L <- function(data,
                           syndrome_col,
                           syndrome_name,
                           patient_col,
                           visit_type_col,
                           inpatient_value = "Inpatient",
                           outpatient_value = "Outpatient",
                           adjustment_ref,
                           facility_col = NULL,
                           facility_name = NULL) {
  # ── Input validation ──────────────────────────────────────────────────────
  required_cols <- c(syndrome_col, patient_col, visit_type_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  if (!all(c("infectious_syndrome", "adjustment_factor_on_CFR") %in%
    names(adjustment_ref))) {
    stop("adjustment_ref must have columns: infectious_syndrome, adjustment_factor_on_CFR.")
  }
  if (!is.null(facility_name) && is.null(facility_col)) {
    stop("facility_col must be provided when facility_name is specified.")
  }

  # ── Look up adjustment type ───────────────────────────────────────────────
  adj_type <- adjustment_ref %>%
    dplyr::filter(.data[["infectious_syndrome"]] == syndrome_name) %>%
    dplyr::pull(.data[["adjustment_factor_on_CFR"]])

  if (length(adj_type) == 0) {
    warning(sprintf(
      "Syndrome '%s' not found in adjustment_ref. Returning CR_L = 1.",
      syndrome_name
    ))
    return(data.frame(
      syndrome        = syndrome_name,
      adjustment_type = "Unknown",
      CR_L            = 1
    ))
  }
  adj_type <- trimws(adj_type[1])

  if (adj_type == "None") {
    message(sprintf("CR_L for '%s': adjustment type = None → CR_L = 1.", syndrome_name))
    return(data.frame(
      syndrome        = syndrome_name,
      adjustment_type = "None",
      CR_L            = 1
    ))
  }

  # ── Filter to syndrome ────────────────────────────────────────────────────
  df <- data %>%
    dplyr::filter(.data[[syndrome_col]] == syndrome_name)

  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
  }

  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  # ── Inpatient ratio ───────────────────────────────────────────────────────
  if (adj_type == "Inpatient ratio") {
    result <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(fac_grp, patient_col)))) %>%
      dplyr::summarise(
        has_ip = any(.data[[visit_type_col]] == inpatient_value),
        .groups = "drop"
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
      dplyr::summarise(
        n_inpatient = sum(has_ip),
        n_total     = dplyr::n(),
        CR_L        = dplyr::if_else(n_total > 0, n_inpatient / n_total, NA_real_),
        .groups     = "drop"
      ) %>%
      dplyr::mutate(
        syndrome        = syndrome_name,
        adjustment_type = "Inpatient ratio"
      )

    # ── Outpatient to inpatient ratio ─────────────────────────────────────────
  } else if (adj_type == "Outpatient to inpatient ratio") {
    result <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(fac_grp, patient_col)))) %>%
      dplyr::summarise(
        has_op = any(.data[[visit_type_col]] == outpatient_value),
        has_ip = any(.data[[visit_type_col]] == inpatient_value),
        .groups = "drop"
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(fac_grp))) %>%
      dplyr::summarise(
        n_op_patients = sum(has_op),
        n_op_to_ip = sum(has_op & has_ip),
        CR_L = dplyr::if_else(
          n_op_patients > 0,
          n_op_to_ip / n_op_patients,
          NA_real_
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        syndrome        = syndrome_name,
        adjustment_type = "Outpatient to inpatient ratio"
      )
  } else {
    warning(sprintf("Unknown adjustment type '%s'. Returning CR_L = 1.", adj_type))
    return(data.frame(
      syndrome        = syndrome_name,
      adjustment_type = adj_type,
      CR_L            = 1
    ))
  }

  message(sprintf(
    "CR_L for '%s' (%s):%s mean = %.3f.",
    syndrome_name,
    adj_type,
    if (!is.null(facility_col)) sprintf(" %d facility/facilities,", nrow(result)) else "",
    mean(result$CR_L, na.rm = TRUE)
  ))
  return(result)
}


# ── Incidence (formula-based) ──────────────────────────────────────────────────

#' Calculate syndrome incidence from deaths, CFR, and CR_L (formula-based)
#'
#' Estimates the number of incident cases of syndrome L using:
#'
#' \deqn{I_L = \frac{D_L}{\text{CFR}_L \times \text{CR}_L}}
#'
#' where the syndrome-level CFR is the pathogen-weighted average:
#'
#' \deqn{\text{CFR}_L = \sum_k P'_{Lk} \times \text{CFR}_{Lk}}
#'
#' Use this when you have population- or facility-level death counts and want
#' to back-calculate incidence. For a direct patient count from facility data,
#' use \code{count_incident_cases()} instead.
#'
#' @param deaths_L       Numeric scalar (pooled mode) or data frame with a
#'   \code{facility_col} column and a \code{deaths_col} column
#'   (facility-level mode).
#' @param cfr_lk_tbl     Data frame with at minimum columns \code{pathogen_col}
#'   and \code{cfr_col}. Typically the \code{cfr_table} element from
#'   \code{calculate_cfr_lk()}.
#' @param P_Lk_prime_tbl Data frame with \code{pathogen_col} and \code{plk_col}.
#'   Use the \code{P_Lk_prime} (pooled) element from
#'   \code{calculate_P_Lk_prime_BSI()} or \code{calculate_P_Lk_prime()}.
#' @param CR_L           Numeric scalar. CFR adjustment factor from
#'   \code{calculate_CR_L()}. Default \code{1} (no adjustment).
#' @param pathogen_col   Character. Pathogen column in both tables.
#'   Default \code{"pathogen"}.
#' @param cfr_col        Character. CFR column in \code{cfr_lk_tbl}.
#'   Default \code{"CFR_LK"}.
#' @param plk_col        Character. P'LK column in \code{P_Lk_prime_tbl}.
#'   Default \code{"P_Lk_prime"}.
#' @param facility_col   Character or NULL. Facility identifier. When provided,
#'   \code{cfr_lk_tbl} and \code{P_Lk_prime_tbl} must each contain
#'   \code{facility_col}, and \code{deaths_L} must be a data frame with
#'   \code{facility_col} + \code{deaths_col}.
#' @param deaths_col     Character. Column in \code{deaths_L} data frame
#'   containing death counts. Default \code{"deaths"}. Ignored when
#'   \code{deaths_L} is a scalar.
#'
#' @return Data frame with columns \code{deaths}, \code{CFR_L}, \code{CR_L},
#'   \code{I_L} (incident cases), and \code{facility_col} when applicable.
#' @export

calculate_incidence_L <- function(deaths_L,
                                  cfr_lk_tbl,
                                  P_Lk_prime_tbl,
                                  CR_L = 1,
                                  pathogen_col = "pathogen",
                                  cfr_col = "CFR_LK",
                                  plk_col = "P_Lk_prime",
                                  facility_col = NULL,
                                  deaths_col = "deaths") {
  # ── Input validation ──────────────────────────────────────────────────────
  for (tbl_name in c("cfr_lk_tbl", "P_Lk_prime_tbl")) {
    tbl <- get(tbl_name)
    if (!pathogen_col %in% names(tbl)) {
      stop(sprintf("'%s' not found in %s.", pathogen_col, tbl_name))
    }
  }
  if (!cfr_col %in% names(cfr_lk_tbl)) {
    stop(sprintf("cfr_col '%s' not found in cfr_lk_tbl.", cfr_col))
  }
  if (!plk_col %in% names(P_Lk_prime_tbl)) {
    stop(sprintf("plk_col '%s' not found in P_Lk_prime_tbl.", plk_col))
  }
  if (!is.null(facility_col)) {
    if (!facility_col %in% names(cfr_lk_tbl)) {
      stop(sprintf("facility_col '%s' not found in cfr_lk_tbl.", facility_col))
    }
    if (!facility_col %in% names(P_Lk_prime_tbl)) {
      stop(sprintf("facility_col '%s' not found in P_Lk_prime_tbl.", facility_col))
    }
    if (!is.data.frame(deaths_L)) {
      stop("deaths_L must be a data frame with facility_col when facility_col is provided.")
    }
    if (!facility_col %in% names(deaths_L)) {
      stop(sprintf("facility_col '%s' not found in deaths_L.", facility_col))
    }
    if (!deaths_col %in% names(deaths_L)) {
      stop(sprintf("deaths_col '%s' not found in deaths_L.", deaths_col))
    }
  }

  # ── Join P'LK and CFR_LK ─────────────────────────────────────────────────
  join_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  merged <- dplyr::inner_join(P_Lk_prime_tbl, cfr_lk_tbl, by = join_cols)

  n_unmatched <- nrow(P_Lk_prime_tbl) - nrow(merged)
  if (n_unmatched > 0) {
    warning(sprintf(
      "%d pathogen(s) in P_Lk_prime_tbl have no matching CFR and are excluded from CFR_L.",
      n_unmatched
    ))
  }

  # ── CFR_L = sum_K( P'_LK × CFR_LK ) per facility ────────────────────────
  grp <- if (!is.null(facility_col)) facility_col else character(0)

  cfr_L_tbl <- merged %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      CFR_L = sum(.data[[plk_col]] * .data[[cfr_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # ── Attach deaths and compute I_L ─────────────────────────────────────────
  if (!is.null(facility_col)) {
    result <- dplyr::left_join(cfr_L_tbl, deaths_L, by = facility_col) %>%
      dplyr::rename(deaths = dplyr::all_of(deaths_col)) %>%
      dplyr::mutate(
        CR_L = CR_L,
        I_L = dplyr::if_else(
          CFR_L * CR_L > 0,
          deaths / (CFR_L * CR_L),
          NA_real_
        )
      )
  } else {
    if (is.data.frame(deaths_L)) {
      stop("deaths_L must be a scalar when facility_col is NULL.")
    }
    result <- cfr_L_tbl %>%
      dplyr::mutate(
        deaths = deaths_L,
        CR_L = CR_L,
        I_L = dplyr::if_else(
          CFR_L * CR_L > 0,
          deaths / (CFR_L * CR_L),
          NA_real_
        )
      )
  }

  message(sprintf(
    "Incidence I_L computed: CFR_L = %.4f, CR_L = %.4f → I_L = %.1f%s.",
    mean(result$CFR_L, na.rm = TRUE),
    mean(result$CR_L, na.rm = TRUE),
    mean(result$I_L, na.rm = TRUE),
    if (!is.null(facility_col)) sprintf(" (mean across %d facility/facilities)", nrow(result)) else ""
  ))
  return(result)
}


# ══════════════════════════════════════════════════════════════════════════════
# LOS-BASED RR AND PAF FOR YLD ATTRIBUTABLE TO AMR
# ══════════════════════════════════════════════════════════════════════════════
#
# Implements two procedures for estimating RR_LOS(k, c):
#
#   Procedure 1 — fit_los_rr_nima()
#     Distribution fitting (Weibull / Lognormal / Gamma) on drug-level R vs S
#     LOS vectors per centre. Produces one overall RR per pathogen. Validation.
#
#   Procedure 2 — fit_los_rr_poisson()
#     Quasi-Poisson regression on class-level binary wide matrix, with HAI as
#     covariate. Produces per-class RR(k, c) with 95% CI. Primary PAF input.
#     Two model options:
#       "pooled_fe"  (default) — one model across all centres with centre FE
#       "per_centre"           — per-centre models, RRs pooled afterwards
#
#   LOS computation:
#     HAI: LOS = date_discharge - date_of_first_positive_culture
#     CAI: LOS = date_discharge - date_of_admission
#     HAI/CAI derived from type_of_infection; NULL / "Not known" rows are
#     classified by the gap (culture - admission): <= threshold -> CAI, else HAI.
#
#   Downstream:
#     assign_rr_to_profiles() -- max rule: RR_kd = max RR_kc for c in C_R(d)
#     compute_paf_los()       -- PAF_LOS(k,d) and overall PAF_k
#
# NOTE (stated limitation): when syndrome_name is supplied, RR_LOS is
#   syndrome-specific. It is then applied to all profiles of pathogen k,
#   assuming syndrome-invariant LOS prolongation across infection sources.
#   Set syndrome_name = NULL for a universal RR pooled over all syndromes.
#
# References:
#   Antimicrobial Resistance Collaborators. Lancet. 2022.

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Compute analytical mean LOS from a fitdistrplus fit object
#' @keywords internal
.compute_mean_from_fit <- function(fit, dist) {
  tryCatch(
    {
      if (dist == "weibull") {
        k <- fit$estimate["shape"]
        l <- fit$estimate["scale"]
        l * gamma(1 + 1 / k)
      } else if (dist == "lnorm") {
        mu <- fit$estimate["meanlog"]
        sg <- fit$estimate["sdlog"]
        exp(mu + sg^2 / 2)
      } else if (dist == "gamma") {
        a <- fit$estimate["shape"]
        r <- fit$estimate["rate"]
        a / r
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
}

#' @keywords internal
.safe_fitdist <- function(x, dist) safe_fit(x, dist)


# ── Step 1a ───────────────────────────────────────────────────────────────────

#' Derive Infection Type (HAI / CAI) per Patient
#'
#' Classifies each row as HAI or CAI. Uses \code{infection_type_col} when it
#' contains a valid value. For rows where that column is \code{NA},
#' \code{"Not known"}, or \code{"NULL"}, derives the classification from the
#' gap between \code{date_culture_col} and \code{date_admission_col}:
#' \itemize{
#'   \item gap <= \code{hai_threshold_hours} -> \strong{CAI}
#'   \item gap >  \code{hai_threshold_hours} -> \strong{HAI}
#' }
#'
#' @param data Data frame.
#' @param infection_type_col Character. Raw infection type column.
#'   Default \code{"type_of_infection"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Date of first positive culture.
#'   Default \code{"date_of_first_positive_culture"}.
#' @param hai_threshold_hours Numeric. Gap threshold in hours. Default \code{48}.
#'
#' @return \code{data} with column \code{infection_type_derived}
#'   (\code{"HAI"} / \code{"CAI"} / \code{"Unknown"}).
#' @export
derive_infection_type <- function(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
) {
  required <- c(infection_type_col, date_admission_col, date_culture_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  # Identify patients with ambiguous infection type AND at least one missing date.
  # These patients cannot have HAI/CAI inferred and are assigned "Not Known".
  ambiguous_mask <- {
    inf_raw_check <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[infection_type_col]])
    ))
    inf_raw_check %in% c("NOT KNOWN", "NOT_KNOWN", "UNKNOWN", "NULL", "NA", "") |
      is.na(data[[infection_type_col]])
  }
  missing_admit <- is.na(data[[date_admission_col]])
  missing_culture <- is.na(data[[date_culture_col]])
  cannot_infer <- ambiguous_mask & (missing_admit | missing_culture)

  if (any(cannot_infer)) {
    n_cannot <- sum(cannot_infer)
    message(sprintf(
      paste0(
        "Cannot infer infection type for %d patient(s): ",
        "infection type is unknown/null AND at least one date is missing. ",
        "Assigning 'Not Known'."
      ),
      n_cannot
    ))
    if (patient_id_col %in% names(data)) {
      flagged <- data[cannot_infer, patient_id_col, drop = TRUE]
      message("  Patient IDs with missing date(s):")
      message(paste0("    ", paste(flagged, collapse = ", ")))
    } else {
      message(
        "  Row indices with missing date(s): ",
        paste(which(cannot_infer), collapse = ", ")
      )
    }
    n_miss_admit <- sum(ambiguous_mask & missing_admit)
    n_miss_culture <- sum(ambiguous_mask & missing_culture)
    message(sprintf(
      "  Breakdown: missing %s = %d | missing %s = %d",
      date_admission_col, n_miss_admit,
      date_culture_col,   n_miss_culture
    ))
  }

  data <- data %>%
    dplyr::mutate(
      .inf_raw = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[infection_type_col]])
      )),
      .gap_h = as.numeric(difftime(
        as.Date(.data[[date_culture_col]]),
        as.Date(.data[[date_admission_col]]),
        units = "hours"
      )),
      .cannot_infer = .inf_raw %in% c(
        "NOT KNOWN", "NOT_KNOWN", "UNKNOWN",
        "NULL", "NA", ""
      ) |
        is.na(.data[[infection_type_col]]),
      infection_type_derived = dplyr::case_when(
        .inf_raw %in% c(
          "HAI", "HOSPITAL ACQUIRED",
          "HOSPITAL-ACQUIRED", "HOSPITAL_ACQUIRED"
        ) |
          grepl("HOSPITAL.ACQUIRED|HEALTH.CARE.ASSOCIATED|HEALTHCARE.ASSOCIATED|\\bHAI\\b",
            .inf_raw,
            perl = TRUE
          ) ~ "HAI",
        .inf_raw %in% c(
          "CAI", "COMMUNITY ACQUIRED",
          "COMMUNITY-ACQUIRED", "COMMUNITY_ACQUIRED"
        ) |
          grepl("COMMUNITY.ACQUIRED|\\bCAI\\b",
            .inf_raw,
            perl = TRUE
          ) ~ "CAI",
        # Ambiguous type (Not Known/NA) but BOTH dates present: infer from gap
        .cannot_infer &
          !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours, "CAI", "HAI"),
        # Ambiguous type and at least one date missing: cannot infer
        .cannot_infer ~ "Not Known",
        # Unrecognised label but BOTH dates present: infer from gap
        !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours, "CAI", "HAI"),
        # Unrecognised label and dates missing: cannot resolve
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::select(-".inf_raw", -".gap_h", -".cannot_infer")

  n_hai <- sum(data$infection_type_derived == "HAI", na.rm = TRUE)
  n_cai <- sum(data$infection_type_derived == "CAI", na.rm = TRUE)
  n_not_known <- sum(data$infection_type_derived == "Not Known", na.rm = TRUE)
  n_unknown <- sum(data$infection_type_derived == "Unknown", na.rm = TRUE)
  message(sprintf(
    "Infection type derived: HAI = %d | CAI = %d | Not Known = %d | Unknown = %d.",
    n_hai, n_cai, n_not_known, n_unknown
  ))
  return(data)
}


# ── Step 1b ───────────────────────────────────────────────────────────────────

#' Compute Patient-Level Post-Infection LOS
#'
#' Computes LOS with infection-type-specific clock start:
#' \itemize{
#'   \item \strong{CAI}: LOS = date_discharge - date_admission
#'   \item \strong{HAI}: LOS = date_discharge - date_culture
#'   \item \strong{Unknown / Not known}: LOS = date_discharge - date_culture
#'     (if culture date present), else date_discharge - date_admission. If
#'     neither reference date is available the patient is excluded.
#' }
#' Only discharged patients are retained. A checkpoint reports how many are
#' missing a discharge date; those patients fall back to \code{los_col} when
#' provided, otherwise they are excluded. Rows with LOS <= 0 or
#' LOS > \code{max_los} are dropped. Returns one row per patient.
#'
#' @param data Data frame (after \code{derive_infection_type()} has been run).
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param organism_col Character. Column containing organism/pathogen names.
#'   Used for episode grouping: a new episode is only created when the culture
#'   date gap exceeds 14 days AND the organism differs from the episode's first
#'   organism. Default \code{"organism_name"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_discharge_col Character. Default \code{"final_outcome_date"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_value Character. Default \code{"Discharged"}.
#' @param infection_type_derived_col Character. Column from
#'   \code{derive_infection_type()}. Default \code{"infection_type_derived"}.
#' @param syndrome_col Character. Syndrome column name. Only used when
#'   \code{syndrome_name} is not \code{NULL}. Default \code{"syndrome"}.
#' @param syndrome_name Character or \code{NULL}. If provided, only patients
#'   with this syndrome are retained before LOS computation. Default \code{NULL}.
#' @param los_col Character or \code{NULL}. Optional name of a pre-computed
#'   LOS column (e.g. recorded directly in the data). Used only for patients
#'   whose discharge date is missing. Values that are \code{NA}, zero, or
#'   negative are treated as invalid and those patients are excluded.
#'   Default \code{NULL}.
#' @param max_los Numeric. Upper cap on LOS in days; patients exceeding this
#'   are excluded. Default \code{200}.
#'
#' @return Data frame: one row per patient with \code{patient_id_col},
#'   \code{facility_col}, \code{infection_type_derived_col}, \code{LOS_days}.
#' @export
compute_patient_los <- function(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  syndrome_name = NULL,
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  infection_type_derived_col = "infection_type_derived",
  los_col = NULL,
  max_los = 200
) {
  # ── Column validation ──────────────────────────────────────────────────────
  required <- c(
    patient_id_col, facility_col, organism_col,
    date_admission_col, date_discharge_col, date_culture_col,
    final_outcome_col, infection_type_derived_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  if (!is.null(los_col) && !los_col %in% names(data)) {
    stop(sprintf("los_col '%s' not found in data.", los_col))
  }

  # ── Facility filtering ─────────────────────────────────────────────────────
  if (!is.null(facility_name)) {
    avail <- unique(data[[facility_col]])
    if (!facility_name %in% avail) {
      stop(sprintf(
        "facility_name '%s' not found in column '%s'. Available: %s",
        facility_name, facility_col,
        paste(sort(as.character(avail)), collapse = ", ")
      ))
    }
    data <- data[data[[facility_col]] == facility_name, ]
    message(sprintf(
      "Filtering to facility '%s': %d row(s) retained.",
      facility_name, nrow(data)
    ))
  } else {
    message(sprintf(
      "Computing LOS across all facilities; '%s' column included in output.",
      facility_col
    ))
  }

  # ── Syndrome filtering ────────────────────────────────────────────────────
  if (!is.null(syndrome_name)) {
    avail_syn <- unique(data[[syndrome_col]])
    if (!syndrome_name %in% avail_syn) {
      stop(sprintf(
        "syndrome_name '%s' not found in column '%s'. Available: %s",
        syndrome_name, syndrome_col,
        paste(sort(as.character(avail_syn)), collapse = ", ")
      ))
    }
    data <- data[data[[syndrome_col]] == syndrome_name, ]
    message(sprintf(
      "Filtering to syndrome '%s': %d row(s) retained.",
      syndrome_name, nrow(data)
    ))
  }

  # ── Checkpoint: discharge date coverage ───────────────────────────────────
  discharged_data <- data[data[[final_outcome_col]] == final_outcome_value, ]
  n_total <- nrow(discharged_data)
  n_no_disc <- sum(is.na(as.Date(discharged_data[[date_discharge_col]])))

  if (n_no_disc > 0L) {
    message(sprintf(
      "Checkpoint: %d of %d discharged patient(s) are missing a discharge date.",
      n_no_disc, n_total
    ))
    if (!is.null(los_col)) {
      message(sprintf(
        "  -> Will fall back to '%s' column for those patients (values <= 0 or NA will be excluded).",
        los_col
      ))
    } else {
      message(
        "  -> No los_col provided; patients with missing discharge date will be excluded."
      )
    }
  }

  # ── Pre-compute episode-minimum culture dates ──────────────────────────────
  # A new episode starts only when BOTH conditions are met:
  #   1. Gap from the current episode's start date > 14 days
  #   2. Organism differs from the current episode's first organism
  # If the same organism returns after > 14 days it is treated as the same
  # infection event; the original (first) culture date is used for LOS.
  episode_min_cult <- data %>%
    dplyr::filter(
      .data[[final_outcome_col]] == final_outcome_value,
      !is.na(.data[[date_culture_col]])
    ) %>%
    dplyr::distinct(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(date_culture_col),
      !!rlang::sym(organism_col)
    ) %>%
    dplyr::mutate(.raw_cult = as.Date(.data[[date_culture_col]])) %>%
    dplyr::arrange(!!rlang::sym(patient_id_col), .raw_cult) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::mutate(
      .ep_id = {
        n <- dplyr::n()
        ep <- integer(n)
        org_vec <- .data[[organism_col]]
        ep_start <- .raw_cult[1L]
        ep_org <- org_vec[1L]
        current_ep <- 0L
        for (i in seq_len(n)) {
          if (i > 1L &&
            as.numeric(.raw_cult[i] - ep_start) > 14L &&
            !identical(org_vec[i], ep_org)) {
            current_ep <- current_ep + 1L
            ep_start <- .raw_cult[i]
            ep_org <- org_vec[i]
          }
          ep[i] <- current_ep
        }
        ep
      }
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col), .ep_id) %>%
    dplyr::mutate(.cult_ep_min = min(.raw_cult)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(date_culture_col),
      !!rlang::sym(organism_col),
      .cult_ep_min
    )

  # ── LOS computation ───────────────────────────────────────────────────────
  unknown_types_lc <- c("unknown", "not known")

  los_data <- data %>%
    dplyr::filter(.data[[final_outcome_col]] == final_outcome_value) %>%
    dplyr::left_join(episode_min_cult,
      by = c(patient_id_col, date_culture_col, organism_col)
    ) %>%
    dplyr::mutate(
      .adm = as.Date(.data[[date_admission_col]]),
      .disc = as.Date(.data[[date_discharge_col]]),
      .cult = .cult_ep_min,
      .inf_lc = tolower(trimws(as.character(.data[[infection_type_derived_col]]))),
      LOS_days = dplyr::case_when(
        # ── Patients WITH a discharge date ────────────────────────────
        # HAI: clock starts at culture date
        !is.na(.disc) & .inf_lc == "hai" ~
          as.numeric(.disc - .cult),
        # CAI: clock starts at admission date
        !is.na(.disc) & .inf_lc == "cai" ~
          as.numeric(.disc - .adm),
        # Unknown / Not known: prefer culture date, fall back to admission date
        !is.na(.disc) & .inf_lc %in% unknown_types_lc & !is.na(.cult) ~
          as.numeric(.disc - .cult),
        !is.na(.disc) & .inf_lc %in% unknown_types_lc & is.na(.cult) & !is.na(.adm) ~
          as.numeric(.disc - .adm),
        # Unknown / Not known with no reference date: exclude
        !is.na(.disc) & .inf_lc %in% unknown_types_lc ~
          NA_real_,
        # All other infection types with discharge date: fallback to admission
        !is.na(.disc) ~
          as.numeric(.disc - .adm),
        # ── Patients WITHOUT a discharge date: resolved below ─────────
        TRUE ~ NA_real_
      )
    )

  # For patients missing a discharge date, fall back to los_col if provided.
  # Values <= 0 or NA in los_col are treated as invalid and excluded.
  if (!is.null(los_col)) {
    los_data <- los_data %>%
      dplyr::mutate(
        LOS_days = dplyr::if_else(
          is.na(.disc),
          dplyr::if_else(
            !is.na(.data[[los_col]]) & as.numeric(.data[[los_col]]) > 0,
            as.numeric(.data[[los_col]]),
            NA_real_
          ),
          LOS_days
        )
      )
  }

  # ── Final filtering and deduplication ─────────────────────────────────────
  los_data <- los_data %>%
    dplyr::select(-".adm", -".disc", -".cult", -".cult_ep_min", -".inf_lc") %>%
    dplyr::filter(!is.na(LOS_days), LOS_days > 0L, LOS_days <= max_los) %>%
    dplyr::distinct(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(facility_col),
      !!rlang::sym(infection_type_derived_col),
      LOS_days
    )

  message(sprintf(
    "Patient LOS computed: %d patients retained (LOS > 0 and <= %d days).",
    nrow(los_data), max_los
  ))
  return(los_data)
}


# ── Step 1c (internal) ────────────────────────────────────────────────────────

#' Collapse drug-level data to antibiotic-class binary wide matrix
#'
#' Standardises antibiotic values (Intermediate -> R/S per GBD rules),
#' collapses to class level (class = 1 if ANY drug in class = R), and pivots
#' wide: one row per patient, one column per antibiotic class (0/1).
#' Column names are sanitised with make.names(); the original-to-safe name
#' mapping is stored as the "class_name_map" attribute.
#'
#' @keywords internal
.build_class_resistance_wide <- function(
  data,
  patient_id_col = "PatientInformation_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  untested_fill = 0L
) {
  abx_std <- data %>%
    dplyr::filter(
      !is.na(.data[[antibiotic_value_col]]),
      !is.na(.data[[antibiotic_class_col]])
    ) %>%
    dplyr::mutate(
      .abx_name = stringr::str_to_lower(stringr::str_trim(
        as.character(.data[[antibiotic_name_col]])
      )),
      .abx_val = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[antibiotic_value_col]])
      )),
      .abx_val = dplyr::case_when(
        .abx_name == "colistin" & .abx_val == "I" ~ "S",
        .abx_val == "I" ~ "R",
        .abx_val %in% c("R", "RESISTANT") ~ "R",
        .abx_val %in% c("S", "SUSCEPTIBLE", "SENSITIVE") ~ "S",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.abx_val))

  class_long <- abx_std %>%
    dplyr::group_by(
      .data[[patient_id_col]],
      .data[[antibiotic_class_col]]
    ) %>%
    dplyr::summarise(
      class_R = as.integer(any(.abx_val == "R")),
      .groups = "drop"
    )

  orig_classes <- sort(unique(class_long[[antibiotic_class_col]]))
  safe_classes <- make.names(orig_classes)
  class_name_map <- setNames(orig_classes, safe_classes)

  class_wide <- class_long %>%
    dplyr::mutate(
      .class_safe = make.names(.data[[antibiotic_class_col]])
    ) %>%
    tidyr::pivot_wider(
      id_cols     = !!rlang::sym(patient_id_col),
      names_from  = ".class_safe",
      values_from = "class_R",
      values_fill = untested_fill
    )

  attr(class_wide, "class_name_map") <- class_name_map
  return(class_wide)
}


# ── Procedure 1 ───────────────────────────────────────────────────────────────

#' Estimate Per-Profile LOS Relative Risk via Distribution Fitting (Nima Procedure)
#'
#' For each pathogen k and each resistance profile delta from
#' \code{compute_resistance_profiles()}, collects the LOS vector of patients
#' whose class-level resistance pattern matches delta, fits Weibull, Lognormal,
#' and Gamma distributions (best by AIC), and computes:
#'   RR_LOS_k_delta = E[LOS | delta] / E[LOS | delta_0]
#' where delta_0 is the all-susceptible profile ("SSS...S").
#'
#' LOS uses the HAI/CAI-specific clock (via \code{compute_patient_los}):
#'   HAI: date_discharge - date_of_first_positive_culture
#'   CAI: date_discharge - date_of_admission
#' HAI/CAI is derived from \code{infection_type_col}; rows where it is
#' NA / "Not known" are classified by the culture-admission gap.
#'
#' Resistance is classified at antibiotic class level (class = R if any drug
#' in that class is R), matching the profile classes from
#' \code{compute_resistance_profiles()}. Patients not tested for all profile
#' classes are excluded from profile assignment.
#'
#' When \code{facility_name} is supplied the data are pre-filtered to that
#' facility; the \code{resistance_profiles} argument should correspondingly be
#' the facility-specific output of \code{compute_resistance_profiles()}.
#' When \code{facility_name = NULL} all facilities are pooled.
#'
#' @param data Data frame. Merged AMR dataset at isolate x antibiotic level.
#' @param resistance_profiles Named list returned by
#'   \code{compute_resistance_profiles()}. One entry per pathogen.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param facility_name Character or \code{NULL}. When supplied, data are
#'   filtered to this facility before computation. Default \code{NULL}.
#' @param organism_col Character. Default \code{"organism_name"}.
#' @param syndrome_col Character. Default \code{"syndrome"}.
#' @param infection_type_col Character. Raw infection-type column used for
#'   HAI/CAI derivation. Default \code{"type_of_infection"}.
#' @param antibiotic_class_col Character. Default \code{"antibiotic_class"}.
#' @param antibiotic_name_col Character. Default \code{"antibiotic_name"}.
#' @param antibiotic_value_col Character. Default \code{"antibiotic_value"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_discharge_col Character. Default \code{"final_outcome_date"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_value Character. Default \code{"Discharged"}.
#' @param syndrome_name Character or NULL. Filter to syndrome. NULL = all.
#' @param organism_name Character or NULL. Pathogen(s) to process. NULL = all.
#' @param distributions Character vector. Distributions to try.
#'   Default \code{c("weibull","lnorm","gamma")}. Ignored when
#'   \code{los_summary = "median"}.
#' @param los_summary Character. How to summarise each profile's LOS vector
#'   before computing the RR. \code{"mean"} (default) fits the best
#'   parametric distribution (Weibull / Lognormal / Gamma by AIC) and
#'   returns the analytical mean. \code{"median"} uses the empirical
#'   median directly — no distribution fitting, more robust to outliers
#'   but ignores the shape of the LOS distribution.
#' @param hai_threshold_hours Numeric. Gap threshold (hours) for HAI/CAI
#'   derivation. Default \code{48}.
#' @param max_los Numeric. Cap on LOS (days). Default \code{200}.
#' @param min_n Integer. Minimum patients required per profile to fit a
#'   distribution. Profiles below this threshold return \code{NA} for RR.
#'   Default \code{10}.
#'
#' @return Named list (one entry per pathogen k):
#' \describe{
#'   \item{\code{RR_k_delta}}{Data frame with columns: \code{profile},
#'     \code{probability}, \code{n_patients}, \code{best_dist},
#'     \code{mean_LOS}, \code{RR_LOS} (relative to all-S profile).}
#'   \item{\code{classes}}{Character vector of antibiotic class names
#'     (from resistance_profiles).}
#'   \item{\code{reference_profile}}{The all-susceptible profile label
#'     (e.g. \code{"SSSS"}).}
#'   \item{\code{mean_LOS_S}}{Fitted mean LOS for the all-S reference profile.}
#'   \item{\code{best_dist_S}}{Best-fit distribution for the all-S profile.}
#'   \item{\code{syndrome_scope}}{Syndrome filter applied, or \code{"all"}.}
#'   \item{\code{facility_scope}}{Facility filter applied, or \code{"all"}.}
#' }
#' @export
fit_los_rr_nima <- function(
  data,
  resistance_profiles,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  facility_name = NULL,
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  syndrome_name = NULL,
  organism_name = NULL,
  distributions = c("weibull", "lnorm", "gamma"),
  los_summary = c("mean", "median"),
  hai_threshold_hours = 48,
  max_los = 365,
  min_n = 10L
) {
  # ── Input validation ───────────────────────────────────────────────────────
  los_summary <- match.arg(los_summary)
  if (los_summary == "mean" && !requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required for los_summary='mean': install.packages('fitdistrplus')")
  }
  if (!is.list(resistance_profiles) || length(resistance_profiles) == 0L) {
    stop("resistance_profiles must be a non-empty named list from compute_resistance_profiles().")
  }

  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col,
    date_culture_col, final_outcome_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  # ── Facility filtering ─────────────────────────────────────────────────────
  if (!is.null(facility_name)) {
    avail <- unique(data[[facility_col]])
    if (!facility_name %in% avail) {
      stop(sprintf(
        "facility_name '%s' not found in column '%s'. Available: %s",
        facility_name, facility_col,
        paste(sort(as.character(avail)), collapse = ", ")
      ))
    }
    data <- data[data[[facility_col]] == facility_name, ]
    message(sprintf(
      "Filtering to facility '%s': %d row(s) retained.",
      facility_name, nrow(data)
    ))
  } else {
    message(sprintf(
      "Computing relative LOS (Nima) across all facilities; '%s' column present.",
      facility_col
    ))
  }

  # ── Derive infection type (HAI / CAI) ──────────────────────────────────────
  data <- derive_infection_type(
    data,
    infection_type_col  = infection_type_col,
    date_admission_col  = date_admission_col,
    date_culture_col    = date_culture_col,
    hai_threshold_hours = hai_threshold_hours
  )

  # ── Filter to discharged patients (+ optional syndrome) ───────────────────
  df <- dplyr::filter(data, .data[[final_outcome_col]] == final_outcome_value)
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .data[[syndrome_col]] == syndrome_name)
  }

  if (nrow(df) == 0L) {
    warning(sprintf(
      "No discharged patients remain after filters (outcome='%s'%s). Returning empty list.",
      final_outcome_value,
      if (!is.null(syndrome_name)) sprintf(", syndrome='%s'", syndrome_name) else ""
    ))
    return(list())
  }

  # ── Checkpoint: discharge date coverage & HAI/CAI split ───────────────────
  n_total <- nrow(df)
  n_no_disc <- sum(is.na(as.Date(df[[date_discharge_col]])))
  if (n_no_disc > 0L) {
    message(sprintf(
      "Checkpoint: %d of %d discharged patient(s) are missing a discharge date — excluded.",
      n_no_disc, n_total
    ))
  }
  inf_types_lc <- tolower(trimws(as.character(df[["infection_type_derived"]])))
  message(sprintf(
    "Checkpoint infection type: HAI=%d | CAI=%d | Other=%d (of %d discharged).",
    sum(inf_types_lc == "hai", na.rm = TRUE),
    sum(inf_types_lc == "cai", na.rm = TRUE),
    sum(!inf_types_lc %in% c("hai", "cai"), na.rm = TRUE),
    n_total
  ))

  # ── Resolve pathogens ──────────────────────────────────────────────────────
  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  no_profile <- setdiff(pathogens, names(resistance_profiles))
  if (length(no_profile) > 0L) {
    message(sprintf(
      "Note: %d pathogen(s) skipped — not in resistance_profiles: %s",
      length(no_profile), paste(no_profile, collapse = ", ")
    ))
  }
  pathogens <- intersect(pathogens, names(resistance_profiles))
  if (length(pathogens) == 0L) {
    stop("No pathogens remain after intersecting with resistance_profiles names.")
  }

  # ── Inner helper: LOS point estimate (mean via distribution fit OR median) ──
  # Returns list(estimate, best_dist).
  # los_summary == "mean"  : fit best parametric distribution by AIC,
  #                          return the analytical mean of that distribution.
  # los_summary == "median": return the empirical median directly;
  #                          best_dist is NA (no fitting performed).
  .compute_los_estimate <- function(los_vec, dists, method) {
    if (method == "median") {
      med <- stats::median(los_vec, na.rm = TRUE)
      return(list(
        estimate = if (is.finite(med)) med else NA_real_,
        best_dist = NA_character_
      ))
    }
    # method == "mean": parametric distribution fitting
    fits <- Filter(
      Negate(is.null),
      setNames(lapply(dists, .safe_fitdist, x = los_vec), dists)
    )
    if (length(fits) == 0L) {
      if (length(unique(los_vec)) == 1L) {
        message(sprintf(
          "Distribution fitting failed: all %d LOS values are identical (%.4g). Consider los_summary='median'.",
          length(los_vec), los_vec[1L]
        ))
      } else {
        message(sprintf(
          "Distribution fitting failed for all %d candidate distribution(s) on %d observations.",
          length(dists), length(los_vec)
        ))
      }
      return(list(estimate = NA_real_, best_dist = NA_character_))
    }
    aics <- sapply(fits, `[[`, "aic")
    best <- names(which.min(aics))
    list(
      estimate = .compute_mean_from_fit(fits[[best]], best),
      best_dist = best
    )
  }

  out <- list()

  for (path in pathogens) {
    path_df <- df %>%
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
          stringr::str_to_lower(stringr::str_trim(path))
      )
    if (nrow(path_df) == 0L) {
      message(sprintf("'%s': no data after filters, skipping.", path))
      next
    }

    # ── HAI/CAI-specific LOS (delegates to compute_patient_los) ───────────
    los_pat <- compute_patient_los(
      path_df,
      patient_id_col             = patient_id_col,
      facility_col               = facility_col,
      facility_name              = NULL, # already filtered above
      date_admission_col         = date_admission_col,
      date_discharge_col         = date_discharge_col,
      date_culture_col           = date_culture_col,
      final_outcome_col          = final_outcome_col,
      final_outcome_value        = final_outcome_value,
      infection_type_derived_col = "infection_type_derived",
      max_los                    = max_los
    )

    if (nrow(los_pat) == 0L) {
      message(sprintf("'%s': no patients with valid LOS, skipping.", path))
      next
    }

    # ── Class-level resistance wide matrix (patient × class, binary) ───────
    resist_wide <- .build_class_resistance_wide(
      path_df,
      patient_id_col       = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col  = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col
    )

    # ── Map profile class names → safe column names in resist_wide ─────────
    # compute_resistance_profiles() uses original names; resist_wide uses
    # make.names()-sanitised names in the same alphabetical order.
    prof_classes_orig <- resistance_profiles[[path]]$classes
    prof_classes_safe <- make.names(prof_classes_orig)

    missing_cls <- setdiff(prof_classes_safe, names(resist_wide))
    if (length(missing_cls) > 0L) {
      message(sprintf(
        "'%s': %d class(es) from profiles absent in patient resistance data (%s) — skipping.",
        path, length(missing_cls), paste(missing_cls, collapse = ", ")
      ))
      next
    }

    # ── Join LOS + class-resistance, assign profile label per patient ──────
    # Profile label = "R"/"S" per class in prof_classes_safe order.
    # Patients with any untested class (NA in resist_wide) are excluded.
    model_data <- dplyr::inner_join(
      los_pat, resist_wide,
      by = patient_id_col
    ) %>%
      dplyr::mutate(
        patient_profile = apply(
          dplyr::select(., dplyr::all_of(prof_classes_safe)),
          1L,
          function(row) {
            if (anyNA(row)) {
              return(NA_character_)
            }
            paste(ifelse(row == 1L, "R", "S"), collapse = "")
          }
        )
      ) %>%
      dplyr::filter(!is.na(patient_profile))

    if (nrow(model_data) == 0L) {
      message(sprintf(
        "'%s': no patients with complete class-resistance data after joining — skipping.",
        path
      ))
      next
    }

    # ── Reference: all-susceptible profile ────────────────────────────────
    all_s_label <- paste(rep("S", length(prof_classes_safe)), collapse = "")
    los_S_ref <- model_data %>%
      dplyr::filter(patient_profile == all_s_label) %>%
      dplyr::pull(LOS_days)

    if (length(los_S_ref) == 0L) {
      warning(sprintf(
        "'%s': no all-susceptible patients found (profile label '%s') — relative LOS framework requires a susceptible reference; skipping.",
        path, all_s_label
      ))
      next
    }
    if (length(los_S_ref) < min_n) {
      warning(sprintf(
        "'%s': all-S profile has only %d patient(s) (min_n=%d) — insufficient reference sample; skipping.",
        path, length(los_S_ref), min_n
      ))
      next
    }

    est_S_ref <- .compute_los_estimate(los_S_ref, distributions, los_summary)
    if (is.na(est_S_ref$estimate) || est_S_ref$estimate <= 0) {
      reason <- if (is.na(est_S_ref$estimate)) {
        "estimation returned NA (distribution fitting may have failed)"
      } else {
        sprintf("%s LOS estimate is %.4g (must be > 0)", los_summary, est_S_ref$estimate)
      }
      warning(sprintf(
        "'%s': all-S reference LOS could not be used as denominator — %s; skipping.",
        path, reason
      ))
      next
    }
    los_S_estimate <- est_S_ref$estimate

    # ── Per-profile RR_LOS_k_delta computation ────────────────────────────
    prof_df <- resistance_profiles[[path]]$profiles
    profile_rr_rows <- vector("list", nrow(prof_df))

    for (i in seq_len(nrow(prof_df))) {
      delta_label <- prof_df$profile[i]
      delta_prob <- prof_df$probability[i]

      los_delta <- model_data %>%
        dplyr::filter(patient_profile == delta_label) %>%
        dplyr::pull(LOS_days)

      n_delta <- length(los_delta)

      if (n_delta == 0L) {
        message(sprintf(
          "'%s' profile '%s': 0 patients matched — possible label mismatch or no observations for this profile.",
          path, delta_label
        ))
        profile_rr_rows[[i]] <- data.frame(
          profile = delta_label,
          probability = round(delta_prob, 6L),
          n_patients = 0L,
          best_dist = NA_character_,
          LOS_estimate = NA_real_,
          RR_LOS = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }
      if (n_delta < min_n) {
        message(sprintf(
          "'%s' profile '%s': only %d patient(s) (min_n=%d) — relative LOS set to NA.",
          path, delta_label, n_delta, min_n
        ))
        profile_rr_rows[[i]] <- data.frame(
          profile = delta_label,
          probability = round(delta_prob, 6L),
          n_patients = n_delta,
          best_dist = NA_character_,
          LOS_estimate = NA_real_,
          RR_LOS = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      est_delta <- .compute_los_estimate(los_delta, distributions, los_summary)
      rr_delta <- if (!is.na(est_delta$estimate) && est_delta$estimate > 0) {
        est_delta$estimate / los_S_estimate
      } else {
        NA_real_
      }

      if (!is.na(rr_delta)) {
        if (rr_delta > 10) {
          message(sprintf(
            "'%s' profile '%s': extreme relative LOS = %.2f (LOS_delta=%.2f / LOS_S=%.2f) — check data.",
            path, delta_label, rr_delta,
            est_delta$estimate, los_S_estimate
          ))
        } else if (rr_delta < 1) {
          message(sprintf(
            "'%s' profile '%s': relative LOS = %.4f < 1 (resistant patients have shorter LOS than susceptible) — verify.",
            path, delta_label, rr_delta
          ))
        }
      }

      profile_rr_rows[[i]] <- data.frame(
        profile = delta_label,
        probability = round(delta_prob, 6L),
        n_patients = n_delta,
        best_dist = est_delta$best_dist, # NA when los_summary="median"
        LOS_estimate = round(est_delta$estimate, 3L),
        RR_LOS = round(rr_delta, 4L),
        stringsAsFactors = FALSE
      )
    }

    rr_k_delta <- dplyr::bind_rows(profile_rr_rows)

    class_cols <- setdiff(names(prof_df), c("profile", "probability"))

    rr_k_delta <- rr_k_delta %>%
      dplyr::left_join(
        prof_df[, c("profile", class_cols)],
        by = "profile"
      )

    n_rr_valid <- sum(!is.na(rr_k_delta$RR_LOS))

    if (n_rr_valid == 0L) {
      warning(sprintf(
        "'%s': all %d profile(s) have NA relative LOS — no usable relative LOS values produced for this pathogen.",
        path, nrow(rr_k_delta)
      ))
    }

    if (all(rr_k_delta$n_patients == 0L)) {
      warning(sprintf(
        "'%s': every profile matched 0 patients — profile labels in resistance_profiles may not align with patient_profile labels derived from class columns (check make.names() sanitisation).",
        path
      ))
    }

    message(sprintf(
      "'%s': %d/%d profile(s) with relative LOS (method=%s). Reference '%s': %s_LOS=%.2f%s, n=%d.",
      path, n_rr_valid, nrow(rr_k_delta), los_summary,
      all_s_label, los_summary, round(los_S_estimate, 2L),
      if (los_summary == "mean") {
        sprintf(" (%s)", est_S_ref$best_dist)
      } else {
        ""
      },
      length(los_S_ref)
    ))

    out[[path]] <- list(
      RR_k_delta        = rr_k_delta,
      classes           = prof_classes_orig,
      reference_profile = all_s_label,
      LOS_S_estimate    = round(los_S_estimate, 3L),
      best_dist_S       = est_S_ref$best_dist, # NA when los_summary="median"
      los_summary       = los_summary,
      syndrome_scope    = if (!is.null(syndrome_name)) syndrome_name else "all",
      facility_scope    = if (!is.null(facility_name)) facility_name else "all"
    )
  }

  return(out)
}


# ── Procedure 2 ───────────────────────────────────────────────────────────────

#' Estimate Per-Class LOS Relative Risk via Quasi-Poisson Regression
#'
#' For each pathogen k and antibiotic class c, fits:
#'   log(E[LOS_i]) = b0 + b_c * class_ci + b_HAI * HAI_i + sum_f gf * I(centre=f)
#' and returns RR_kc_LOS = exp(b_c) with 95% Wald CI.
#'
#' LOS is computed with HAI/CAI-specific clock (see compute_patient_los).
#' Resistance is classified at antibiotic class level (class = 1 if any drug
#' in class = R). HAI enters as a binary covariate to absorb residual baseline
#' severity differences beyond the LOS clock adjustment. No PreDays covariate
#' is used — the LOS clock change from admission to culture date for HAI
#' patients already removes the pre-infection period from the outcome.
#'
#' Stated limitation: when syndrome_name is supplied, RR_kc_LOS is
#' syndrome-specific. Downstream PAF applies this RR to all profiles of
#' pathogen k, assuming syndrome-invariant LOS prolongation across infection
#' sources. Set syndrome_name = NULL for a universal RR.
#'
#' @param data Data frame. Merged AMR dataset at isolate x antibiotic level.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param isolate_id_col Character. Default \code{"isolate_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param organism_col Character. Default \code{"organism_name"}.
#' @param syndrome_col Character. Default \code{"syndrome"}.
#' @param infection_type_col Character. Default \code{"type_of_infection"}.
#' @param antibiotic_class_col Character. Default \code{"antibiotic_class"}.
#' @param antibiotic_name_col Character. Default \code{"antibiotic_name"}.
#' @param antibiotic_value_col Character. Default \code{"antibiotic_value"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_discharge_col Character. Default \code{"final_outcome_date"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_value Character. Default \code{"Discharged"}.
#' @param syndrome_name Character or NULL. Filter to syndrome (e.g. "BSI").
#'   NULL = all syndromes (universal RR).
#' @param organism_name Character or NULL. Pathogen(s) to process. NULL = all.
#' @param hai_threshold_hours Numeric. HAI/CAI gap threshold. Default \code{48}.
#' @param centre_model Character. \code{"pooled_fe"} (default) or
#'   \code{"per_centre"}.
#' @param pool_method Character. Used when centre_model = "per_centre".
#'   \code{"weighted_mean"} (default) or \code{"median"}.
#' @param max_los Numeric. Default \code{200}.
#' @param min_n Integer. Minimum patients per model. Default \code{10}.
#'
#' @return Data frame: pathogen, antibiotic_class, RR_LOS, CI_lower, CI_upper,
#'   n_patients, centre_model, syndrome_scope. When centre_model = "per_centre",
#'   attribute "per_centre_rr" contains per-centre estimates.
#' @export
fit_los_rr_poisson <- function(
  data,
  patient_id_col = "PatientInformation_id",
  isolate_id_col = "isolate_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  syndrome_name = NULL,
  organism_name = NULL,
  hai_threshold_hours = 48,
  centre_model = c("pooled_fe", "per_centre"),
  pool_method = c("weighted_mean", "median"),
  max_los = 365,
  min_n = 10L
) {
  centre_model <- match.arg(centre_model)
  pool_method <- match.arg(pool_method)

  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col,
    date_culture_col, final_outcome_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }

  # Step A: derive infection type
  data <- derive_infection_type(
    data,
    infection_type_col  = infection_type_col,
    date_admission_col  = date_admission_col,
    date_culture_col    = date_culture_col,
    hai_threshold_hours = hai_threshold_hours
  )

  df <- dplyr::filter(data, .data[[final_outcome_col]] == final_outcome_value)
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .data[[syndrome_col]] == syndrome_name)
  }

  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  all_rr <- list()
  per_ctr_rr <- list()

  # ── inner helpers ─────────────────────────────────────────────────────────
  .fit_one_class <- function(sub_data, cls_col, fac_col, use_fac_fe) {
    vals <- sub_data[[cls_col]]
    if (length(unique(stats::na.omit(vals))) < 2L) {
      return(NULL)
    }
    if (nrow(sub_data) < min_n) {
      return(NULL)
    }
    fmla <- if (use_fac_fe &&
      length(unique(sub_data[[fac_col]])) > 1L) {
      stats::as.formula(
        sprintf("LOS_days ~ `%s` + HAI + %s", cls_col, fac_col)
      )
    } else {
      stats::as.formula(sprintf("LOS_days ~ `%s` + HAI", cls_col))
    }
    tryCatch(
      suppressWarnings(
        stats::glm(fmla,
          data = sub_data,
          family = stats::quasipoisson(link = "log")
        )
      ),
      error = function(e) NULL
    )
  }

  .extract_rr <- function(model, cls_col) {
    if (is.null(model)) {
      return(NULL)
    }
    coefs <- stats::coef(model)
    vcmat <- stats::vcov(model)
    pat <- paste0("^`?", gsub(".", "\\.", cls_col, fixed = TRUE), "`?$")
    cname <- grep(pat, names(coefs), value = TRUE)
    if (length(cname) == 0L) {
      return(NULL)
    }
    beta <- coefs[[cname[1L]]]
    se <- sqrt(vcmat[cname[1L], cname[1L]])
    data.frame(
      RR = exp(beta),
      CI_lower = exp(beta - 1.96 * se),
      CI_upper = exp(beta + 1.96 * se),
      n = nrow(model$model)
    )
  }

  # ── main loop ─────────────────────────────────────────────────────────────
  for (path in pathogens) {
    path_df <- df %>%
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
          stringr::str_to_lower(stringr::str_trim(path))
      )
    if (nrow(path_df) == 0L) {
      message(sprintf("'%s': no data after filters, skipping.", path))
      next
    }

    # Step B: patient LOS
    los_data <- compute_patient_los(
      path_df,
      patient_id_col             = patient_id_col,
      facility_col               = facility_col,
      date_admission_col         = date_admission_col,
      date_discharge_col         = date_discharge_col,
      date_culture_col           = date_culture_col,
      final_outcome_col          = final_outcome_col,
      final_outcome_value        = final_outcome_value,
      infection_type_derived_col = "infection_type_derived",
      max_los                    = max_los
    )

    # Step C: class resistance wide matrix
    resist_wide <- .build_class_resistance_wide(
      path_df,
      patient_id_col       = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col  = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col
    )
    class_name_map <- attr(resist_wide, "class_name_map")
    class_safe <- setdiff(names(resist_wide), patient_id_col)

    model_data <- dplyr::inner_join(
      los_data, resist_wide,
      by = patient_id_col
    ) %>%
      dplyr::mutate(
        HAI = as.integer(infection_type_derived == "HAI")
      )

    if (nrow(model_data) < min_n) {
      message(sprintf(
        "'%s': only %d patients after merging, skipping.",
        path, nrow(model_data)
      ))
      next
    }

    # Step D: Poisson per class
    if (centre_model == "pooled_fe") {
      for (cls in class_safe) {
        orig_class <- class_name_map[[cls]]
        model <- .fit_one_class(model_data, cls, facility_col, TRUE)
        rr_df <- .extract_rr(model, cls)
        if (is.null(rr_df)) {
          message(sprintf(
            "'%s' | class '%s': model failed or no variation, skipping.",
            path, orig_class
          ))
          next
        }
        all_rr[[length(all_rr) + 1L]] <- data.frame(
          pathogen = path,
          antibiotic_class = orig_class,
          RR_LOS = round(rr_df$RR, 4L),
          CI_lower = round(rr_df$CI_lower, 4L),
          CI_upper = round(rr_df$CI_upper, 4L),
          n_patients = rr_df$n,
          centre_model = "pooled_fe",
          syndrome_scope = if (!is.null(syndrome_name)) {
            syndrome_name
          } else {
            "all"
          },
          stringsAsFactors = FALSE
        )
      }
    } else { # per_centre

      centres <- sort(unique(model_data[[facility_col]]))
      path_ctr_list <- list()

      for (cls in class_safe) {
        orig_class <- class_name_map[[cls]]
        ctr_rr_rows <- list()

        for (ctr in centres) {
          ctr_df <- dplyr::filter(
            model_data,
            .data[[facility_col]] == ctr
          )
          model <- .fit_one_class(ctr_df, cls, facility_col, FALSE)
          rr_df <- .extract_rr(model, cls)
          if (is.null(rr_df)) next
          ctr_rr_rows[[length(ctr_rr_rows) + 1L]] <- data.frame(
            pathogen         = path,
            antibiotic_class = orig_class,
            centre           = ctr,
            RR_centre        = rr_df$RR,
            CI_lower_centre  = rr_df$CI_lower,
            CI_upper_centre  = rr_df$CI_upper,
            n_centre         = rr_df$n,
            stringsAsFactors = FALSE
          )
        }

        if (length(ctr_rr_rows) == 0L) {
          message(sprintf(
            "'%s' | class '%s': no centre model succeeded, skipping.",
            path, orig_class
          ))
          next
        }

        ctr_rr_df <- dplyr::bind_rows(ctr_rr_rows)
        path_ctr_list[[orig_class]] <- ctr_rr_df

        # Pool on log scale
        log_rr <- log(ctr_rr_df$RR_centre)
        log_se <- (log(ctr_rr_df$CI_upper_centre) -
          log(ctr_rr_df$CI_lower_centre)) / (2 * 1.96)
        wts <- ctr_rr_df$n_centre

        if (pool_method == "weighted_mean") {
          ok <- !is.na(log_rr) & !is.na(log_se)
          pool_log_rr <- sum(wts[ok] * log_rr[ok]) / sum(wts[ok])
          pool_log_se <- sqrt(sum(wts[ok]^2 * log_se[ok]^2)) /
            sum(wts[ok])
        } else {
          pool_log_rr <- stats::median(log_rr, na.rm = TRUE)
          pool_log_se <- stats::median(log_se, na.rm = TRUE)
        }

        all_rr[[length(all_rr) + 1L]] <- data.frame(
          pathogen = path,
          antibiotic_class = orig_class,
          RR_LOS = round(exp(pool_log_rr), 4L),
          CI_lower = round(exp(pool_log_rr -
            1.96 * pool_log_se), 4L),
          CI_upper = round(exp(pool_log_rr +
            1.96 * pool_log_se), 4L),
          n_patients = sum(ctr_rr_df$n_centre),
          centre_model = "per_centre",
          syndrome_scope = if (!is.null(syndrome_name)) {
            syndrome_name
          } else {
            "all"
          },
          stringsAsFactors = FALSE
        )
      }
      per_ctr_rr[[path]] <- path_ctr_list
    }

    n_cls_fitted <- sum(sapply(all_rr, function(x) x$pathogen == path))
    message(sprintf(
      "'%s': %d class relative LOS(s) estimated (model=%s, syndrome=%s).",
      path, n_cls_fitted, centre_model,
      if (!is.null(syndrome_name)) syndrome_name else "all"
    ))
  }

  if (length(all_rr) == 0L) {
    warning("No relative LOS values could be estimated. Check data, filters, and min_n.")
    return(data.frame())
  }

  result <- dplyr::bind_rows(all_rr)
  if (length(per_ctr_rr) > 0L) {
    attr(result, "per_centre_rr") <- per_ctr_rr
  }
  return(result)
}


# ── Step 2b ───────────────────────────────────────────────────────────────────

#' Estimate Per-Class LOS Relative Risk via Parametric Distribution Fitting
#'
#' For each pathogen x antibiotic class, splits patients into resistant (R) and
#' susceptible (S) groups, fits a parametric distribution (default: gamma) to
#' each group's LOS, and returns:
#'   RR_LOS = mean_LOS(R) / mean_LOS(S)
#'
#' This is an alternative to \code{fit_los_rr_poisson()} that avoids regression
#' and models the full LOS distribution per group. The output flat data frame
#' is fully compatible with \code{assign_rr_to_profiles()} and
#' \code{filter_profiles_to_rr_classes()}.
#'
#' @param data Data frame of isolate/patient rows (long format).
#' @param patient_id_col Character. Patient identifier column.
#' @param facility_col Character. Facility/centre column.
#' @param organism_col Character. Organism name column.
#' @param syndrome_col Character. Syndrome column (used only when
#'   \code{syndrome_name} is not \code{NULL}).
#' @param infection_type_col Character. Column used to derive HAI/CAI.
#' @param antibiotic_class_col Character. Antibiotic class column.
#' @param antibiotic_name_col Character. Antibiotic name column.
#' @param antibiotic_value_col Character. Antibiotic value column (S/I/R).
#' @param date_admission_col Character. Date of admission column.
#' @param date_discharge_col Character. Date of discharge column.
#' @param date_culture_col Character. Date of culture column.
#' @param final_outcome_col Character. Final outcome column.
#' @param final_outcome_value Character. Value indicating discharge.
#' @param syndrome_name Character or \code{NULL}. Restrict to this syndrome.
#' @param organism_name Character vector or \code{NULL}. Restrict to these
#'   pathogen(s); otherwise all pathogens in the filtered data.
#' @param hai_threshold_hours Numeric. Hours after admission before a culture
#'   is classified as HAI. Default \code{48}.
#' @param distributions Character vector. Candidate distributions to fit.
#'   Default \code{"gamma"}.
#' @param max_los Numeric. Maximum plausible LOS in days. Default \code{200}.
#' @param min_n Integer. Minimum patients required in both R and S groups.
#'   Default \code{10L}.
#'
#' @return A data frame with one row per pathogen x class:
#'   \describe{
#'     \item{pathogen}{Pathogen name.}
#'     \item{antibiotic_class}{Antibiotic class name.}
#'     \item{RR_LOS}{Fitted mean ratio: mean_LOS(R) / mean_LOS(S).}
#'     \item{mean_LOS_R}{Fitted mean LOS for resistant patients.}
#'     \item{mean_LOS_S}{Fitted mean LOS for susceptible patients.}
#'     \item{best_dist_R}{Distribution fitted to the R group.}
#'     \item{best_dist_S}{Distribution fitted to the S group.}
#'     \item{n_R}{Number of resistant patients.}
#'     \item{n_S}{Number of susceptible patients.}
#'     \item{syndrome_scope}{Syndrome filter applied, or "all".}
#'   }
#'
#' @seealso \code{\link{fit_los_rr_poisson}}, \code{\link{assign_rr_to_profiles}}
#' @export
fit_los_rr_distribution <- function(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  date_admission_col = "date_of_admission",
  date_discharge_col = "final_outcome_date",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_value = "Discharged",
  syndrome_name = NULL,
  organism_name = NULL,
  facility_name = NULL,
  hai_threshold_hours = 48,
  distributions = "gamma",
  max_los = 365,
  min_n = 10L
) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required: install.packages('fitdistrplus')")
  }

  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_discharge_col,
    date_culture_col, final_outcome_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # ── Derive infection type (HAI / CAI) ──────────────────────────────────────
  data <- derive_infection_type(
    data,
    infection_type_col  = infection_type_col,
    date_admission_col  = date_admission_col,
    date_culture_col    = date_culture_col,
    hai_threshold_hours = hai_threshold_hours
  )

  # ── Filter to discharged patients + optional syndrome + optional facility ──
  df <- dplyr::filter(data, .data[[final_outcome_col]] == final_outcome_value)
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .data[[syndrome_col]] == syndrome_name)
  }
  if (!is.null(facility_name)) {
    df <- dplyr::filter(df, .data[[facility_col]] == facility_name)
  }

  if (nrow(df) == 0L) {
    warning(sprintf(
      "No discharged patients remain after filters (outcome='%s'%s%s). Returning empty data frame.",
      final_outcome_value,
      if (!is.null(syndrome_name)) sprintf(", syndrome='%s'", syndrome_name) else "",
      if (!is.null(facility_name)) sprintf(", facility='%s'", facility_name) else ""
    ))
    return(data.frame())
  }

  # ── Resolve pathogens ──────────────────────────────────────────────────────
  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  # ── Inner helper: fit distribution(s) and return analytical mean ───────────
  .fit_dist_mean <- function(los_vec) {
    fits <- Filter(
      Negate(is.null),
      setNames(
        lapply(distributions, .safe_fitdist, x = los_vec),
        distributions
      )
    )
    if (length(fits) == 0L) {
      return(list(mean = NA_real_, dist = NA_character_))
    }
    aics <- sapply(fits, `[[`, "aic")
    best <- names(which.min(aics))
    list(mean = .compute_mean_from_fit(fits[[best]], best), dist = best)
  }

  all_rr <- list()

  for (path in pathogens) {
    path_df <- df %>%
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(.data[[organism_col]])) ==
          stringr::str_to_lower(stringr::str_trim(path))
      )
    if (nrow(path_df) == 0L) {
      message(sprintf("'%s': no data after filters, skipping.", path))
      next
    }

    # ── Patient LOS ────────────────────────────────────────────────────────
    los_pat <- compute_patient_los(
      path_df,
      patient_id_col             = patient_id_col,
      facility_col               = facility_col,
      facility_name              = NULL,
      date_admission_col         = date_admission_col,
      date_discharge_col         = date_discharge_col,
      date_culture_col           = date_culture_col,
      final_outcome_col          = final_outcome_col,
      final_outcome_value        = final_outcome_value,
      infection_type_derived_col = "infection_type_derived",
      max_los                    = max_los
    )
    if (nrow(los_pat) == 0L) {
      message(sprintf("'%s': no patients with valid LOS, skipping.", path))
      next
    }

    # ── Class resistance wide matrix (patient x class, binary 0/1) ──────────
    # untested_fill = NA_integer_: patients not tested for a class get NA,
    # so they are excluded from both n_R and n_S for that class (not treated
    # as susceptible). Only patients actually tested contribute to the fit.
    resist_wide <- .build_class_resistance_wide(
      path_df,
      patient_id_col       = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col  = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col,
      untested_fill        = NA_integer_
    )
    class_name_map <- attr(resist_wide, "class_name_map")
    class_safe <- setdiff(names(resist_wide), patient_id_col)

    model_data <- dplyr::inner_join(los_pat, resist_wide, by = patient_id_col)

    if (nrow(model_data) < min_n) {
      message(sprintf(
        "'%s': only %d patients after merging, skipping.",
        path, nrow(model_data)
      ))
      next
    }

    n_fitted <- 0L

    # ── Per-class R vs S distribution fit ──────────────────────────────────
    for (cls in class_safe) {
      orig_class <- class_name_map[[cls]]

      los_r <- model_data$LOS_days[
        !is.na(model_data[[cls]]) & model_data[[cls]] == 1L
      ]
      los_s <- model_data$LOS_days[
        !is.na(model_data[[cls]]) & model_data[[cls]] == 0L
      ]

      n_r <- length(los_r)
      n_s <- length(los_s)

      if (n_r < min_n || n_s < min_n) {
        message(sprintf(
          "'%s' | class '%s': n_R=%d / n_S=%d — one group below min_n=%d, skipping.",
          path, orig_class, n_r, n_s, min_n
        ))
        next
      }

      est_r <- .fit_dist_mean(los_r)
      est_s <- .fit_dist_mean(los_s)

      if (is.na(est_r$mean) || is.na(est_s$mean) || est_s$mean <= 0) {
        message(sprintf(
          "'%s' | class '%s': distribution fitting failed (mean_R=%s, mean_S=%s), skipping.",
          path, orig_class,
          if (is.na(est_r$mean)) "NA" else round(est_r$mean, 3L),
          if (is.na(est_s$mean)) "NA" else round(est_s$mean, 3L)
        ))
        next
      }

      rr_los <- est_r$mean / est_s$mean

      if (rr_los > 10) {
        message(sprintf(
          "'%s' | class '%s': extreme RR_LOS = %.2f (mean_R=%.2f / mean_S=%.2f) — check data.",
          path, orig_class, rr_los, est_r$mean, est_s$mean
        ))
      } else if (rr_los < 1) {
        message(sprintf(
          "'%s' | class '%s': RR_LOS = %.4f < 1 (resistant patients have shorter LOS than susceptible) — verify.",
          path, orig_class, rr_los
        ))
      }

      all_rr[[length(all_rr) + 1L]] <- data.frame(
        pathogen         = path,
        antibiotic_class = orig_class,
        facility_name    = if (!is.null(facility_name)) facility_name else NA_character_,
        RR_LOS           = round(rr_los, 4L),
        mean_LOS_R       = round(est_r$mean, 3L),
        mean_LOS_S       = round(est_s$mean, 3L),
        mean_LOS_R_years = round(est_r$mean / 365, 6L),
        mean_LOS_S_years = round(est_s$mean / 365, 6L),
        best_dist_R      = est_r$dist,
        best_dist_S      = est_s$dist,
        n_R              = n_r,
        n_S              = n_s,
        syndrome_scope   = if (!is.null(syndrome_name)) syndrome_name else "all",
        stringsAsFactors = FALSE
      )
      n_fitted <- n_fitted + 1L
    }

    message(sprintf(
      "'%s': %d class relative LOS(s) estimated (gamma distribution, syndrome=%s).",
      path, n_fitted,
      if (!is.null(syndrome_name)) syndrome_name else "all"
    ))
  }

  if (length(all_rr) == 0L) {
    warning("No relative LOS values could be estimated. Check data, filters, and min_n.")
    return(data.frame(
      pathogen = character(), antibiotic_class = character(),
      facility_name = character(), RR_LOS = numeric(),
      mean_LOS_R = numeric(), mean_LOS_S = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  return(dplyr::bind_rows(all_rr))
}


# ── Step 3 ────────────────────────────────────────────────────────────────────

#' Assign Per-Class LOS RR to Resistance Profiles (Max Rule)
#'
#' For each resistance profile delta (from compute_resistance_profiles()),
#' determines the profile-level RR_kd_LOS using the GBD max rule:
#'   RR_kd_LOS = max over c in C_R(d) of RR_kc_LOS   [if C_R(d) non-empty]
#'             = 1                                      [if d = all-susceptible]
#' where C_R(d) = {c : d_c = 1}.
#' The CI reported for each profile is that of its dominant (max-RR) class.
#'
#' @param profiles_output Named list from compute_resistance_profiles().
#' @param rr_table Data frame from fit_los_rr_poisson() or fit_los_rr_nima().
#'   Must have columns pathogen_col, class_col, rr_col, and optionally
#'   CI_lower / CI_upper.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param rr_col Character. Default \code{"RR_LOS"}.
#' @param fallback_rr Numeric. RR for resistant classes with no match.
#'   Default \code{1} (no attributable effect).
#'
#' @return Named list (one entry per pathogen): original profiles data frame
#'   augmented with RR_LOS_profile, dominant_class, and (if available)
#'   CI_lower_profile / CI_upper_profile.
#' @export
assign_rr_to_profiles <- function(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1
) {
  if (!is.list(profiles_output) ||
    !all(sapply(
      profiles_output,
      function(x) all(c("profiles", "classes") %in% names(x))
    ))) {
    stop("profiles_output must be the list returned by compute_resistance_profiles().")
  }
  missing_rr <- setdiff(c(pathogen_col, class_col, rr_col), names(rr_table))
  if (length(missing_rr) > 0L) {
    stop(sprintf(
      "Column(s) not found in rr_table: %s",
      paste(missing_rr, collapse = ", ")
    ))
  }

  has_ci <- all(c("CI_lower", "CI_upper") %in% names(rr_table))
  out <- list()

  for (path in names(profiles_output)) {
    prof_list <- profiles_output[[path]]
    profiles <- prof_list$profiles
    classes <- prof_list$classes # alphabetical; matches profile cols

    rr_k <- rr_table[rr_table[[pathogen_col]] == path, ]
    rr_lookup <- setNames(rr_k[[rr_col]], rr_k[[class_col]])
    ci_lo_lookup <- if (has_ci) setNames(rr_k$CI_lower, rr_k[[class_col]]) else NULL
    ci_hi_lookup <- if (has_ci) setNames(rr_k$CI_upper, rr_k[[class_col]]) else NULL

    n_prof <- nrow(profiles)
    rr_profile <- numeric(n_prof)
    dom_class <- character(n_prof)
    ci_lo_prof <- if (has_ci) numeric(n_prof) else NULL
    ci_hi_prof <- if (has_ci) numeric(n_prof) else NULL

    for (i in seq_len(n_prof)) {
      resist_cls <- classes[as.integer(profiles[i, classes]) == 1L]

      if (length(resist_cls) == 0L) {
        rr_profile[i] <- 1.0
        dom_class[i] <- "all_susceptible"
        if (has_ci) {
          ci_lo_prof[i] <- 1.0
          ci_hi_prof[i] <- 1.0
        }
        next
      }

      rrs_c <- sapply(resist_cls, function(cc) {
        ifelse(cc %in% names(rr_lookup), rr_lookup[[cc]], fallback_rr)
      })
      max_idx <- which.max(rrs_c)
      rr_profile[i] <- rrs_c[max_idx]
      dom_class[i] <- resist_cls[max_idx]

      if (has_ci) {
        dc <- dom_class[i]
        ci_lo_prof[i] <- ifelse(dc %in% names(ci_lo_lookup),
          ci_lo_lookup[[dc]], fallback_rr
        )
        ci_hi_prof[i] <- ifelse(dc %in% names(ci_hi_lookup),
          ci_hi_lookup[[dc]], fallback_rr
        )
      }
    }

    profiles$RR_LOS_profile <- round(rr_profile, 4L)
    profiles$dominant_class <- dom_class
    if (has_ci) {
      profiles$CI_lower_profile <- round(ci_lo_prof, 4L)
      profiles$CI_upper_profile <- round(ci_hi_prof, 4L)
    }

    out[[path]] <- profiles
    message(sprintf(
      "'%s': relative LOS assigned to %d profiles. Max relative LOS = %.4f (dominant class: %s).",
      path, n_prof, max(rr_profile), dom_class[which.max(rr_profile)]
    ))
  }

  return(out)
}


# ── Step 3b ───────────────────────────────────────────────────────────────────

#' Filter Profiles to Classes with Actual RR Estimates
#'
#' After \code{assign_rr_to_profiles()}, some profiles may be resistant only to
#' antibiotic classes that have no RR estimate in the 7-centres LOS data (e.g.
#' a class present in the susceptibility dataset but not tested in the LOS
#' cohort). Those profiles receive \code{fallback_rr = 1} and contribute zero
#' to the PAF numerator while their probability mass still enters the
#' denominator, silently distorting the estimate.
#'
#' This function drops any profile where \strong{every} resistant class
#' (\eqn{\delta_c = 1}) is absent from the per-pathogen RR table, then
#' re-normalises the surviving profile probabilities to sum to 1.
#'
#' @param profiles_with_rr Named list from \code{assign_rr_to_profiles()}.
#'   Each element is a data frame with class indicator columns and
#'   \code{RR_LOS_profile}.
#' @param rr_table Data frame from \code{fit_los_rr_poisson()}.
#'   Must have columns \code{pathogen_col} and \code{class_col}.
#' @param pathogen_col Character. Default \code{"pathogen"}.
#' @param class_col Character. Default \code{"antibiotic_class"}.
#' @param probability_col Character. Default \code{"probability"}.
#' @param fallback_rr Numeric. The fallback value used in
#'   \code{assign_rr_to_profiles()}. Default \code{1}.
#'
#' @return Named list with the same structure as \code{profiles_with_rr} but
#'   with unmatched profiles removed and probabilities re-normalised.
#' @export
filter_profiles_to_rr_classes <- function(
  profiles_with_rr,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  probability_col = "probability",
  fallback_rr = 1
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }
  if (!all(c(pathogen_col, class_col) %in% names(rr_table))) {
    stop(sprintf(
      "Columns '%s' and '%s' must be present in rr_table.",
      pathogen_col, class_col
    ))
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    # Classes that have a real RR estimate for this pathogen
    rr_classes <- rr_table[rr_table[[pathogen_col]] == path, class_col]
    rr_classes <- unique(rr_classes[!is.na(rr_classes)])

    # Identify class indicator columns present in this profile data frame
    non_class_cols <- c(
      "profile", probability_col, "RR_LOS_profile",
      "dominant_class", "CI_lower_profile",
      "CI_upper_profile", "numerator", "PAF_LOS",
      "denominator"
    )
    class_cols <- setdiff(names(df), non_class_cols)

    # Matchable classes: those present as indicator columns AND in rr_table
    matchable <- intersect(class_cols, rr_classes)

    if (length(matchable) == 0L) {
      warning(sprintf(
        "'%s': no profile classes overlap with relative LOS table — all profiles dropped.",
        path
      ))
      next
    }

    # Keep profiles where:
    #   (a) the profile is all-susceptible (reference, always kept), OR
    #   (b) at least one matchable class is resistant (== 1)
    # This ensures the all-susceptible reference profile is never dropped.
    all_class_cols <- setdiff(names(df), non_class_cols)
    is_all_s <- if (length(all_class_cols) > 0L) {
      rowSums(df[, all_class_cols, drop = FALSE] == 1L, na.rm = TRUE) == 0L
    } else {
      rep(TRUE, nrow(df))
    }

    if (length(matchable) == 1L) {
      has_matchable_r <- df[[matchable]] == 1L
    } else {
      has_matchable_r <- rowSums(df[, matchable, drop = FALSE] == 1L) > 0L
    }

    keep <- is_all_s | has_matchable_r

    n_before <- nrow(df)
    df_kept <- df[keep, , drop = FALSE]
    n_after <- nrow(df_kept)
    n_dropped <- n_before - n_after

    if (n_after == 0L) {
      warning(sprintf("'%s': all profiles dropped after matching — skipping.", path))
      next
    }

    # Re-normalise probabilities
    prob_sum <- sum(df_kept[[probability_col]], na.rm = TRUE)
    if (prob_sum <= 0) {
      warning(sprintf("'%s': probability sum is zero after filtering — skipping.", path))
      next
    }
    df_kept[[probability_col]] <- df_kept[[probability_col]] / prob_sum

    out[[path]] <- df_kept

    message(sprintf(
      "'%s': %d -> %d profiles kept (%d dropped, no matchable relative LOS class). Prob re-normalised.",
      path, n_before, n_after, n_dropped
    ))
  }

  return(out)
}


# ── Step 4 ────────────────────────────────────────────────────────────────────

#' Compute LOS Population Attributable Fraction per Resistance Profile
#'
#' Computes PAF_kd_LOS for each pathogen k and profile d using the GBD
#' multi-exposure Levin formula:
#'
#'   PAF_kd_LOS = R'_kd * (RR_kd - 1) / [1 + sum_d R'_kd * (RR_kd - 1)]
#'
#' where R'_kd is the profile probability from compute_resistance_profiles()
#' and RR_kd is from assign_rr_to_profiles(). The denominator equals
#' E_d[RR_kd] (expected RR over all profiles) and is shared across profiles
#' for pathogen k. The all-susceptible profile contributes 0 (RR = 1).
#'
#' Overall PAF for pathogen k:
#'   PAF_k = sum_d PAF_kd_LOS
#'         = [sum_d R'_kd(RR_kd-1)] / [1 + sum_d R'_kd(RR_kd-1)]
#'
#' @param profiles_with_rr Named list from assign_rr_to_profiles().
#' @param probability_col Character. Profile probability column.
#'   Default \code{"probability"}.
#' @param rr_profile_col Character. Profile-level RR column.
#'   Default \code{"RR_LOS_profile"}.
#' @param profile_col Character. Profile label column. Default \code{"profile"}.
#'
#' @return Named list (one entry per pathogen):
#'   per_profile: profiles data frame augmented with numerator, PAF_LOS,
#'     denominator.
#'   PAF_k: overall PAF for pathogen k.
#'   denominator: 1 + sum_d R'_kd(RR_kd-1).
#' @export
compute_paf_los <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col, profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_kd  (sums to 1)
    rr <- df[[rr_profile_col]] # RR_LOS_kd
    numerator_vec <- p * (rr - 1.0) # R'_kd * (RR_kd - 1)
    denom <- 1.0 + sum(numerator_vec, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) {
      warning(sprintf(
        "'%s': PAF denominator = %.6g (must be > 0) — all relative LOS may be < 1 or NA; skipping.",
        path, denom
      ))
      next
    }
    paf_vec <- numerator_vec / denom

    df$numerator <- round(numerator_vec, 6L)
    df$PAF_LOS <- round(paf_vec, 6L)
    df$denominator <- round(denom, 6L)

    paf_k <- sum(paf_vec)

    out[[path]] <- list(
      per_profile = df,
      PAF_k       = round(paf_k, 6L),
      denominator = round(denom, 6L)
    )
    message(sprintf(
      "'%s': PAF_k = %.4f | E[relative LOS] (denominator) = %.4f | %d profiles.",
      path, paf_k, denom, nrow(df)
    ))
  }

  return(out)
}


# ── Step 5 ────────────────────────────────────────────────────────────────────

#' Compute Associated-Burden Fractions per Resistance Profile
#'
#' Computes the fraction of total expected YLD burden that *occurs in*
#' infections with each resistance profile delta (YLDs associated with
#' resistance).  This is a **burden partition**, not a counterfactual.
#'
#' For a single drug-class d the formula simplifies to:
#'
#'   Fraction_assoc_Kd = R'_Kd * RR_Kd / [(1 - R'_Kd) + R'_Kd * RR_Kd]
#'
#' With resistance profiles the denominator becomes the expected RR across
#' ALL profiles (including the all-susceptible profile with RR = 1):
#'
#'   E_RR_k = sum_delta  R'_K_delta * RR_K_delta
#'          = 1 + sum_delta  R'_K_delta * (RR_K_delta - 1)   [equivalent]
#'
#' Per-profile associated fraction:
#'   fraction_K_delta = R'_K_delta * RR_K_delta / E_RR_k
#'
#' Overall associated fraction (all resistant profiles combined):
#'   Fraction_k = sum_{delta != 0}  fraction_K_delta
#'
#' where delta != 0 denotes profiles with at least one resistant class.
#'
#' Note: E_RR_k is numerically identical to the `denominator` produced by
#' compute_paf_los() — both equal 1 + sum_d R'_kd*(RR_kd - 1).
#'
#' @param profiles_with_rr Named list from assign_rr_to_profiles() or
#'   filter_profiles_to_rr_classes().  Each entry is a profile data frame.
#' @param probability_col Character.  Profile probability column.
#'   Default \code{"probability"}.
#' @param rr_profile_col Character.  Profile-level RR column.
#'   Default \code{"RR_LOS_profile"}.
#'
#' @return Named list (one entry per pathogen) containing:
#'   \itemize{
#'     \item \code{per_profile}: profile data frame augmented with
#'       \code{numerator_assoc} (= p * rr) and \code{fraction_assoc}
#'       (= p * rr / E_RR_k).
#'     \item \code{Fraction_k}: overall associated fraction for the pathogen
#'       (sum of \code{fraction_assoc} over all resistant profiles).
#'     \item \code{E_RR_k}: expected RR = sum_delta R'_K_delta * RR_K_delta.
#'   }
#' @export
compute_fraction_associated <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  # Columns that are metadata / computed — NOT class indicator columns
  non_class_cols <- c(
    "profile", probability_col, rr_profile_col,
    "dominant_class", "CI_lower_profile", "CI_upper_profile",
    "numerator", "PAF_LOS", "denominator",
    "numerator_assoc", "fraction_assoc"
  )

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_K_delta  (sums to 1 after normalisation)
    rr <- df[[rr_profile_col]] # RR_LOS_K_delta (1.0 for all-susceptible profile)

    # E[RR_k] = sum_delta  R'_K_delta * RR_K_delta
    # Equivalent: 1 + sum_delta R'_K_delta*(RR-1) = PAF denominator
    E_RR_k <- sum(p * rr)

    if (E_RR_k <= 0) {
      warning(sprintf("'%s': E[relative LOS] <= 0 — skipping.", path))
      next
    }

    # Per-profile numerator and fraction
    numerator_assoc <- p * rr
    fraction_assoc_vec <- numerator_assoc / E_RR_k

    df$numerator_assoc <- round(numerator_assoc, 6L)
    df$fraction_assoc <- round(fraction_assoc_vec, 6L)

    # Resistant profiles: at least one binary class indicator column == 1.
    # The all-susceptible profile has every class column = 0 and RR = 1.
    class_cols <- setdiff(names(df), non_class_cols)
    if (length(class_cols) > 0L) {
      is_resistant <- rowSums(
        df[, class_cols, drop = FALSE] == 1L,
        na.rm = TRUE
      ) > 0L
    } else {
      # Fallback when class columns absent: use RR > 1 as proxy
      is_resistant <- rr > 1.0
    }

    Fraction_k <- sum(fraction_assoc_vec[is_resistant])

    out[[path]] <- list(
      per_profile = df,
      Fraction_k  = round(Fraction_k, 6L),
      E_RR_k      = round(E_RR_k, 6L)
    )

    message(sprintf(
      "'%s': Fraction_k (associated) = %.4f | E[relative LOS] = %.4f | %d profiles (%d resistant).",
      path, Fraction_k, E_RR_k, nrow(df), sum(is_resistant)
    ))
  }

  return(out)
}


# ── Step 5b ───────────────────────────────────────────────────────────────────

#' Compute Fatal Prevalence of Resistance (R_kd)
#'
#' Computes the fatal prevalence of resistance R_kd for each pathogen k and
#' resistance profile delta, using per-profile mortality odds ratios
#' (OR_death) from \code{fit_mortality_rr_logistic()} in place of LOS
#' relative risks.  The formula is identical to
#' \code{compute_fraction_associated()} — only the RR source changes.
#'
#' \deqn{E[\text{OR}_k] = \sum_\delta R'_{K\delta} \cdot \text{OR}_{K\delta}}
#'
#' \deqn{R_{kd} = \frac{R'_{K\delta} \cdot \text{OR}_{K\delta}}{E[\text{OR}_k]}}
#'
#' \deqn{R_{k} = \sum_{\delta \neq \delta_0} R_{kd}}
#'
#' where \eqn{R'_{K\delta}} is the profile probability from
#' \code{compute_resistance_profiles()} and \eqn{\text{OR}_{K\delta}} is
#' assigned by \code{assign_rr_to_profiles()} with \code{rr_col = "OR_death"}
#' (the all-susceptible profile carries OR = 1).
#'
#' \strong{Usage pipeline:}
#' \preformatted{
#'   # 1. Fit mortality OR per class
#'   mort_or <- fit_mortality_rr_logistic(data, ...)
#'
#'   # 2. Assign OR to profiles via max rule
#'   profiles_with_or <- assign_rr_to_profiles(
#'       profiles_output,
#'       rr_table = mort_or,
#'       rr_col   = "OR_death"
#'   )
#'
#'   # 3. Compute R_kd (fatal prevalence of resistance)
#'   rkd <- compute_R_kd_fatal(profiles_with_or)
#' }
#'
#' @param profiles_with_rr Named list from \code{assign_rr_to_profiles()}
#'   called with \code{rr_col = "OR_death"}.
#' @param probability_col Character.  Profile probability column.
#'   Default \code{"probability"}.
#' @param rr_profile_col Character.  Profile-level OR column as produced by
#'   \code{assign_rr_to_profiles()}.  Default \code{"RR_LOS_profile"}.
#'
#' @return Named list (one entry per pathogen) containing:
#'   \itemize{
#'     \item \code{per_profile}: profile data frame augmented with
#'       \code{numerator_R_kd} (= \eqn{R'_{K\delta} \cdot \text{OR}_{K\delta}})
#'       and \code{R_kd} (= numerator / E[OR_k]) per profile.
#'     \item \code{R_k}: overall fatal prevalence of resistance for pathogen k
#'       (sum of \code{R_kd} over resistant profiles only).
#'     \item \code{E_OR_k}: expected OR =
#'       \eqn{\sum_\delta R'_{K\delta} \cdot \text{OR}_{K\delta}}.
#'   }
#' @export
compute_R_kd_fatal <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  # Columns that are metadata / computed — NOT class indicator columns
  non_class_cols <- c(
    "profile", probability_col, rr_profile_col,
    "dominant_class", "CI_lower_profile", "CI_upper_profile",
    "numerator", "PAF_LOS", "denominator",
    "numerator_assoc", "fraction_assoc",
    "numerator_R_kd", "R_kd"
  )

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_K_delta  (sums to 1)
    or <- df[[rr_profile_col]] # OR_death_K_delta (1.0 for all-susceptible)

    # E[OR_k] = sum_delta  R'_K_delta * OR_K_delta
    E_OR_k <- sum(p * or, na.rm = TRUE)

    if (!is.finite(E_OR_k) || E_OR_k <= 0) {
      warning(sprintf(
        "'%s': E[OR_death] = %.6g (must be > 0) — all mortality ORs may be <= 1 or NA; skipping.",
        path, E_OR_k
      ))
      next
    }

    numerator_R_kd <- p * or
    R_kd_vec <- numerator_R_kd / E_OR_k

    df$numerator_R_kd <- round(numerator_R_kd, 6L)
    df$R_kd <- round(R_kd_vec, 6L)

    # Resistant profiles: at least one binary class indicator column == 1.
    # The all-susceptible profile has every class column = 0 and OR = 1.
    class_cols <- setdiff(names(df), non_class_cols)
    if (length(class_cols) > 0L) {
      is_resistant <- rowSums(
        df[, class_cols, drop = FALSE] == 1L,
        na.rm = TRUE
      ) > 0L
    } else {
      is_resistant <- or > 1.0
    }

    R_k <- sum(R_kd_vec[is_resistant], na.rm = TRUE)

    out[[path]] <- list(
      per_profile = df,
      R_k         = round(R_k, 6L),
      E_OR_k      = round(E_OR_k, 6L)
    )

    message(sprintf(
      "'%s': R_k (fatal resistance prevalence) = %.4f | E[OR_death] = %.4f | %d profiles (%d resistant).",
      path, R_k, E_OR_k, nrow(df), sum(is_resistant)
    ))
  }

  return(out)
}

# ── Step 6 ────────────────────────────────────────────────────────────────────

#' Compute YLDs Associated with Resistance
#'
#' Multiplies \code{YLD_k} (from \code{calculate_YLD()}) by the associated
#' burden score/fraction (from \code{compute_fraction_associated()}).
#'
#'   YLD_associated_k = YLD_k * Fraction_k
#'
#' @param yld_k_tbl Data frame from \code{calculate_YLD()} containing at
#'   least a pathogen column and a YLD column.
#' @param fraction_assoc_list Named list from \code{compute_fraction_associated()}.
#' @param pathogen_col Character. Pathogen column in \code{yld_k_tbl}.
#'   Default \code{"pathogen"}.
#' @param yld_col Character. YLD column in \code{yld_k_tbl}.
#'   Default \code{"YLD"}.
#' @param probability_col Character. Profile probability column in
#'   \code{per_profile}. Default \code{"probability"}.
#' @param rr_profile_col Character. Profile LOS/RR column in \code{per_profile}.
#'   Default \code{"RR_LOS_profile"}.
#'
#' @return \code{yld_k_tbl} augmented with columns \code{Fraction_k},
#'   \code{R_k_delta}, \code{LOS_k_delta}, and \code{YLD_associated}.
#' @export
compute_yld_associated <- function(
  yld_k_tbl,
  fraction_assoc_list,
  pathogen_col = "pathogen",
  yld_col = "YLD",
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
) {
  if (!is.data.frame(yld_k_tbl)) {
    stop("yld_k_tbl must be a data frame (output of calculate_YLD()).")
  }
  if (!is.list(fraction_assoc_list)) {
    stop("fraction_assoc_list must be the list returned by compute_fraction_associated().")
  }
  if (!pathogen_col %in% names(yld_k_tbl)) {
    stop(sprintf("pathogen_col '%s' not found in yld_k_tbl.", pathogen_col))
  }
  if (!yld_col %in% names(yld_k_tbl)) {
    stop(sprintf("yld_col '%s' not found in yld_k_tbl.", yld_col))
  }

  # Build flat lookup: one row per pathogen
  frac_rows <- lapply(names(fraction_assoc_list), function(k) {
    res <- fraction_assoc_list[[k]]
    if (is.null(res)) {
      return(NULL)
    }

    # Defaults if per_profile unavailable
    R_k_delta_txt <- NA_character_
    LOS_k_delta_txt <- NA_character_

    if (!is.null(res$per_profile) && is.data.frame(res$per_profile)) {
      df <- res$per_profile

      # Classify resistant profiles the same way as Step 5
      non_class_cols <- c(
        "profile", probability_col, rr_profile_col,
        "dominant_class", "CI_lower_profile", "CI_upper_profile",
        "numerator", "PAF_LOS", "denominator",
        "numerator_assoc", "fraction_assoc"
      )

      class_cols <- setdiff(names(df), non_class_cols)
      rr <- if (rr_profile_col %in% names(df)) df[[rr_profile_col]] else rep(NA_real_, nrow(df))

      if (length(class_cols) > 0L) {
        is_resistant <- rowSums(df[, class_cols, drop = FALSE] == 1L, na.rm = TRUE) > 0L
      } else {
        is_resistant <- rr > 1.0
      }

      if (probability_col %in% names(df) && rr_profile_col %in% names(df)) {
        p_res <- df[[probability_col]][is_resistant]
        rr_res <- df[[rr_profile_col]][is_resistant]

        if (length(p_res) > 0L) {
          R_k_delta_txt <- paste(round(p_res, 6L), collapse = ",")
        }
        if (length(rr_res) > 0L) {
          LOS_k_delta_txt <- paste(round(rr_res, 6L), collapse = ",")
        }
      }
    }

    data.frame(
      .pathogen = k,
      Fraction_k = res$Fraction_k,
      R_k_delta = R_k_delta_txt,
      LOS_k_delta = LOS_k_delta_txt,
      stringsAsFactors = FALSE
    )
  })

  frac_df <- do.call(rbind, frac_rows)

  if (is.null(frac_df) || nrow(frac_df) == 0L) {
    stop("fraction_assoc_list is empty — no fractions to join.")
  }

  names(frac_df)[names(frac_df) == ".pathogen"] <- pathogen_col

  out <- merge(yld_k_tbl, frac_df, by = pathogen_col, all.x = TRUE)

  out$YLD_associated <- out[[yld_col]] * out$Fraction_k

  n_unmatched <- sum(is.na(out$Fraction_k))
  if (n_unmatched > 0L) {
    warning(sprintf(
      "%d row(s) in yld_k_tbl had no matching Fraction_k; YLD_associated set to NA.",
      n_unmatched
    ))
  }

  message(sprintf(
    "YLD_associated computed: %d pathogen(s), total YLD_associated = %.4f.",
    sum(!is.na(out$YLD_associated)),
    sum(out$YLD_associated, na.rm = TRUE)
  ))

  return(out)
}


# ── Step 7 ────────────────────────────────────────────────────────────────────

#' Compute YLDs Attributable to Resistance
#'
#' Multiplies \code{YLD_k} (from \code{calculate_YLD()}) by the LOS-based
#' PAF (from \code{compute_paf_los()}).
#'
#'   YLD_attributable_k = YLD_k * PAF_k
#'
#' Answers: "How much disability burden exists *only because* infections were
#' resistant instead of susceptible?"  This is a counterfactual — it measures
#' the excess burden driven purely by resistance.
#'
#' Note: YLD_attributable_k < YLD_associated_k always, because
#'   PAF_k = Fraction_k * (1 - 1/E_RR_k)  <  Fraction_k.
#'
#' @param yld_k_tbl Data frame from \code{calculate_YLD()} containing at
#'   least a pathogen column and a YLD column.
#' @param paf_los_list Named list from \code{compute_paf_los()}.
#' @param pathogen_col Character.  Pathogen column in \code{yld_k_tbl}.
#'   Default \code{"pathogen"}.
#' @param yld_col Character.  YLD column in \code{yld_k_tbl}.
#'   Default \code{"YLD"}.
#'
#' @return \code{yld_k_tbl} augmented with columns \code{PAF_k},
#'   \code{denominator}, and \code{YLD_attributable}.
#' @export
compute_yld_attributable <- function(
  yld_k_tbl,
  paf_los_list,
  pathogen_col = "pathogen",
  yld_col = "YLD"
) {
  if (!is.data.frame(yld_k_tbl)) {
    stop("yld_k_tbl must be a data frame (output of calculate_YLD()).")
  }
  if (!is.list(paf_los_list)) {
    stop("paf_los_list must be the list returned by compute_paf_los().")
  }
  if (!pathogen_col %in% names(yld_k_tbl)) {
    stop(sprintf("pathogen_col '%s' not found in yld_k_tbl.", pathogen_col))
  }
  if (!yld_col %in% names(yld_k_tbl)) {
    stop(sprintf("yld_col '%s' not found in yld_k_tbl.", yld_col))
  }

  # Build flat lookup: one row per pathogen
  paf_rows <- lapply(names(paf_los_list), function(k) {
    res <- paf_los_list[[k]]
    if (is.null(res)) {
      return(NULL)
    }
    data.frame(
      .pathogen = k,
      PAF_k = res$PAF_k,
      denominator = res$denominator,
      stringsAsFactors = FALSE
    )
  })
  paf_df <- do.call(rbind, paf_rows)

  if (is.null(paf_df) || nrow(paf_df) == 0L) {
    stop("paf_los_list is empty — no PAFs to join.")
  }

  names(paf_df)[names(paf_df) == ".pathogen"] <- pathogen_col

  out <- merge(yld_k_tbl, paf_df, by = pathogen_col, all.x = TRUE)

  out$YLD_attributable <- out[[yld_col]] * out$PAF_k

  n_unmatched <- sum(is.na(out$PAF_k))
  if (n_unmatched > 0L) {
    warning(sprintf(
      "%d row(s) in yld_k_tbl had no matching PAF_k; YLD_attributable set to NA.",
      n_unmatched
    ))
  }

  message(sprintf(
    "YLD_attributable computed: %d pathogen(s), total YLD_attributable = %.4f.",
    sum(!is.na(out$YLD_attributable)),
    sum(out$YLD_attributable, na.rm = TRUE)
  ))

  return(out)
}


# ══════════════════════════════════════════════════════════════════════════════
# MORTALITY-BASED RELATIVE RISK FOR AMR BURDEN
# ══════════════════════════════════════════════════════════════════════════════
#
# Mixed-effects logistic regression per antibiotic class per pathogen:
#
#   logit(P(death_i)) = β₀ + β₁·Resistance_ci + β₂·Age_i + β₃·Sex_i
#                     + β₄·HAI_i + β₅·ICU_i + β₆·Comorbidity_i
#                     + u_facility   [u_facility ~ N(0, σ²)]
#
#   OR_death_kc = exp(β₁) : odds ratio for death, resistant vs susceptible,
#                controlling for age, sex, infection acquisition route,
#                disease severity (ICU), comorbidity, and facility clustering.
#
# Pipeline:
#   M1. derive_infection_type_for_mortality()  -- HAI/CAI with death-specific
#                                                 data-quality reporting
#   M2. .derive_icu_binary()                  -- "ever in ICU" binary per patient
#   M3. .encode_comorbidity_mortality()       -- standardise comorbidity column
#   M4. .check_hai_icu_collinearity()         -- collinearity guard (phi coeff)
#   M5. .build_class_resistance_wide()        -- binary class matrix (reused)
#   fit_mortality_rr_logistic()               -- main glmer loop
#
# References:
#   Antimicrobial Resistance Collaborators. Lancet. 2022.


# ── Step M1 ───────────────────────────────────────────────────────────────────

#' Derive Infection Type (HAI / CAI) for Mortality RR Model
#'
#' Variant of \code{derive_infection_type()} tailored for the mortality
#' relative-risk model. Key differences from the LOS version:
#' \enumerate{
#'   \item Processes \strong{all} patients (not just discharged), because the
#'         mortality model requires both dead and surviving patients.
#'   \item Recognises explicit HAI/CAI labels first (and common synonyms such
#'         as "hospital-acquired infection" / "community acquired infection").
#'         Date-gap derivation is applied only when the label is absent or
#'         ambiguous.
#'   \item Performs a dedicated death-date data-quality check: patients whose
#'         \code{final_outcome} equals \code{death_value} but whose outcome
#'         date is missing — yet admission or culture date is present — are
#'         listed. HAI/CAI derivation is unaffected (it uses admission vs
#'         culture date), but the missing death date is surfaced for review.
#' }
#'
#' @param data Data frame.
#' @param infection_type_col Character. Raw infection type column.
#'   Default \code{"type_of_infection"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_date_col Character. Date of final outcome (death date
#'   for deceased patients). Default \code{"final_outcome_date"}.
#' @param death_value Character. Value in \code{final_outcome_col} that
#'   indicates death. Default \code{"Death"}.
#' @param hai_threshold_hours Numeric. Gap threshold in hours between
#'   admission and culture date. Default \code{48}.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with column \code{infection_type_derived}
#'   (\code{"HAI"} / \code{"CAI"} / \code{"Not Known"}).
#'   Attribute \code{"missing_death_date_patients"} is a data frame of rows
#'   where death is confirmed but the outcome date is absent.
#' @export
derive_infection_type_for_mortality <- function(
  data,
  infection_type_col = "type_of_infection",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  death_value = "Death",
  hai_threshold_hours = 48,
  patient_id_col = "PatientInformation_id"
) {
  # ── Column validation ──────────────────────────────────────────────────────
  required <- c(
    infection_type_col, date_admission_col,
    date_culture_col, final_outcome_col
  )
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  has_outcome_date <- final_outcome_date_col %in% names(data)

  # ── Death-date data-quality check ─────────────────────────────────────────
  # Flag deceased patients who are missing their outcome date but do have
  # admission or culture date recorded.  HAI/CAI derivation proceeds normally
  # (it only needs admission vs culture date); this is purely a data-quality
  # report so the analyst can trace incomplete records.
  if (has_outcome_date) {
    is_dead <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[final_outcome_col]])
    )) == stringr::str_to_upper(death_value)

    missing_outcome_dt <- is.na(suppressWarnings(
      as.Date(as.character(data[[final_outcome_date_col]]))
    ))
    has_admit_or_culture <- !is.na(data[[date_admission_col]]) |
      !is.na(data[[date_culture_col]])

    dead_no_date <- is_dead & missing_outcome_dt & has_admit_or_culture

    if (any(dead_no_date)) {
      n_flag <- sum(dead_no_date)
      message(sprintf(
        paste0(
          "[Mortality] Data-quality flag: %d patient row(s) have ",
          "final_outcome='%s' but no outcome date recorded ",
          "(admission/culture dates are present). ",
          "HAI/CAI derivation is unaffected. ",
          "Review these records for date completeness."
        ),
        n_flag, death_value
      ))
      if (patient_id_col %in% names(data)) {
        flagged_ids <- unique(
          data[dead_no_date, patient_id_col, drop = TRUE]
        )
        message(sprintf(
          "  Flagged patient IDs (first 20 of %d): %s",
          length(flagged_ids),
          paste(head(flagged_ids, 20L), collapse = ", ")
        ))
      }
      flagged_df <- data[dead_no_date, , drop = FALSE]
    } else {
      flagged_df <- data[0L, , drop = FALSE]
    }
  } else {
    message(sprintf(
      "[Mortality] Column '%s' not found; skipping death-date quality check.",
      final_outcome_date_col
    ))
    flagged_df <- data[0L, , drop = FALSE]
  }

  # ── Identify ambiguous labels ──────────────────────────────────────────────
  ambiguous_mask <- {
    inf_up <- stringr::str_to_upper(stringr::str_trim(
      as.character(data[[infection_type_col]])
    ))
    inf_up %in% c("NOT KNOWN", "NOT_KNOWN", "UNKNOWN", "NULL", "NA", "") |
      is.na(data[[infection_type_col]])
  }
  missing_admit <- is.na(data[[date_admission_col]])
  missing_culture <- is.na(data[[date_culture_col]])
  cannot_infer <- ambiguous_mask & (missing_admit | missing_culture)

  if (any(cannot_infer)) {
    n_cannot <- sum(cannot_infer)
    message(sprintf(
      paste0(
        "[Mortality] Cannot infer HAI/CAI for %d patient row(s): ",
        "label is ambiguous AND >=1 date is missing. ",
        "Assigned 'Not Known' — these rows are excluded from the ",
        "mortality regression (missing HAI covariate)."
      ),
      n_cannot
    ))
    if (patient_id_col %in% names(data)) {
      ids_flag <- unique(data[cannot_infer, patient_id_col, drop = TRUE])
      message(sprintf(
        "  Patient IDs (first 20 of %d): %s",
        length(ids_flag),
        paste(head(ids_flag, 20L), collapse = ", ")
      ))
    }
    message(sprintf(
      "  Breakdown — missing %s: %d | missing %s: %d",
      date_admission_col, sum(ambiguous_mask & missing_admit),
      date_culture_col,   sum(ambiguous_mask & missing_culture)
    ))
  }

  # ── Derive infection type ──────────────────────────────────────────────────
  # Explicit labels take priority; date-gap only used when label is ambiguous.
  data <- data %>%
    dplyr::mutate(
      .inf_raw = stringr::str_to_upper(stringr::str_trim(
        as.character(.data[[infection_type_col]])
      )),
      .gap_h = as.numeric(difftime(
        as.Date(suppressWarnings(
          as.character(.data[[date_culture_col]])
        )),
        as.Date(suppressWarnings(
          as.character(.data[[date_admission_col]])
        )),
        units = "hours"
      )),
      .cannot_infer = (.inf_raw %in% c(
        "NOT KNOWN", "NOT_KNOWN",
        "UNKNOWN", "NULL", "NA", ""
      ) |
        is.na(.data[[infection_type_col]])),
      infection_type_derived = dplyr::case_when(
        # ── Explicit HAI labels ────────────────────────────────────────
        .inf_raw %in% c(
          "HAI",
          "HOSPITAL ACQUIRED",
          "HOSPITAL-ACQUIRED",
          "HOSPITAL_ACQUIRED",
          "NOSOCOMIAL",
          "HOSPITAL ACQUIRED INFECTION",
          "HOSPITAL-ACQUIRED INFECTION"
        ) ~ "HAI",
        # ── Explicit CAI labels ────────────────────────────────────────
        .inf_raw %in% c(
          "CAI",
          "COMMUNITY ACQUIRED",
          "COMMUNITY-ACQUIRED",
          "COMMUNITY_ACQUIRED",
          "COMMUNITY ACQUIRED INFECTION",
          "COMMUNITY-ACQUIRED INFECTION"
        ) ~ "CAI",
        # ── Date-gap derivation (both dates present) ──────────────────
        .cannot_infer &
          !is.na(.data[[date_admission_col]]) &
          !is.na(.data[[date_culture_col]]) ~
          dplyr::if_else(.gap_h <= hai_threshold_hours,
            "CAI", "HAI"
          ),
        # ── Cannot infer ──────────────────────────────────────────────
        .cannot_infer ~ "Not Known",
        # ── Unrecognised non-empty label ──────────────────────────────
        TRUE ~ "Not Known"
      )
    ) %>%
    dplyr::select(-".inf_raw", -".gap_h", -".cannot_infer")

  n_hai <- sum(data$infection_type_derived == "HAI", na.rm = TRUE)
  n_cai <- sum(data$infection_type_derived == "CAI", na.rm = TRUE)
  n_not_known <- sum(data$infection_type_derived == "Not Known", na.rm = TRUE)
  message(sprintf(
    "[Mortality] Infection type: HAI = %d | CAI = %d | Not Known = %d (excluded from model).",
    n_hai, n_cai, n_not_known
  ))

  attr(data, "missing_death_date_patients") <- flagged_df
  return(data)
}


# ── Step M2 ───────────────────────────────────────────────────────────────────

#' Derive ICU Binary Flag per Patient
#'
#' Collapses unit-type data (one row per drug test per patient) to a
#' patient-level binary ICU indicator using the \strong{"ever in ICU"} rule:
#' if a patient has \emph{any} row with an ICU-type location during that
#' admission, \code{ICU = 1}; otherwise \code{ICU = 0}.
#'
#' Missing unit-type values are coded as \code{NA}. When the proportion of
#' patients with \emph{all} unit-type rows missing exceeds
#' \code{missing_threshold}, the function returns a three-level ordered factor
#' (\code{"ICU"} / \code{"Ward"} / \code{"Unknown"}) and emits a message
#' recommending that the \code{"Unknown"} level be included in the model.
#'
#' @param data Data frame at drug-test level (before patient collapse).
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param unit_type_col Character. Column containing ward/ICU labels.
#'   Default \code{"unit_type"}.
#' @param icu_values Character vector. Values (matched case-insensitively) that
#'   indicate ICU admission.
#'   Default \code{c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU")}.
#' @param missing_threshold Numeric in [0, 1]. Proportion of patients with
#'   entirely missing unit-type above which a 3-level factor is returned.
#'   Default \code{0.10}.
#'
#' @return Patient-level data frame: \code{patient_id_col} and \code{ICU}
#'   (integer 0/1, or factor \code{"ICU"/"Ward"/"Unknown"} when missing is
#'   high).
#' @keywords internal
.derive_icu_binary <- function(
  data,
  patient_id_col = "PatientInformation_id",
  unit_type_col = "unit_type",
  icu_values = c(
    "ICU", "Intensive Care", "Critical Care",
    "PICU", "NICU"
  ),
  missing_threshold = 0.10
) {
  if (!unit_type_col %in% names(data)) {
    message(sprintf(
      "[ICU] Column '%s' not found — ICU covariate will be omitted.",
      unit_type_col
    ))
    return(
      dplyr::distinct(data, !!rlang::sym(patient_id_col)) %>%
        dplyr::mutate(ICU = NA_integer_)
    )
  }

  icu_pattern <- paste(stringr::str_to_lower(icu_values), collapse = "|")

  per_patient <- data %>%
    dplyr::mutate(
      .unit_norm = stringr::str_to_lower(stringr::str_trim(
        as.character(.data[[unit_type_col]])
      )),
      .is_missing = is.na(.data[[unit_type_col]]) |
        .unit_norm %in% c("na", "", "null"),
      .is_icu = !.is_missing &
        grepl(icu_pattern, .unit_norm, fixed = FALSE)
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::summarise(
      .any_icu     = any(.is_icu, na.rm = TRUE),
      .all_missing = all(.is_missing),
      .groups      = "drop"
    )

  n_patients <- nrow(per_patient)
  prop_missing <- sum(per_patient$.all_missing) / n_patients

  if (prop_missing > missing_threshold) {
    message(sprintf(
      paste0(
        "[ICU] %.1f%% of patients have entirely missing unit_type. ",
        "Returning 3-level factor (ICU / Ward / Unknown). ",
        "Consider including 'Unknown' as a separate level in the model."
      ),
      prop_missing * 100
    ))
    per_patient <- per_patient %>%
      dplyr::mutate(
        ICU = dplyr::case_when(
          .all_missing ~ "Unknown",
          .any_icu ~ "ICU",
          TRUE ~ "Ward"
        ),
        ICU = factor(ICU, levels = c("Ward", "ICU", "Unknown"))
      )
    n_icu <- sum(per_patient$ICU == "ICU", na.rm = TRUE)
    n_ward <- sum(per_patient$ICU == "Ward", na.rm = TRUE)
    n_unk <- sum(per_patient$ICU == "Unknown", na.rm = TRUE)
  } else {
    per_patient <- per_patient %>%
      dplyr::mutate(
        ICU = dplyr::case_when(
          .all_missing ~ NA_integer_,
          .any_icu ~ 1L,
          TRUE ~ 0L
        )
      )
    n_icu <- sum(per_patient$ICU == 1L, na.rm = TRUE)
    n_ward <- sum(per_patient$ICU == 0L, na.rm = TRUE)
    n_unk <- sum(is.na(per_patient$ICU))
  }

  message(sprintf(
    "[ICU] Per-patient ICU flag: ICU = %d | Ward = %d | Unknown/NA = %d.",
    n_icu, n_ward, n_unk
  ))

  per_patient %>%
    dplyr::select(!!rlang::sym(patient_id_col), ICU)
}


# ── Step M3 ───────────────────────────────────────────────────────────────────

#' Encode Comorbidity Column for Mortality Model
#'
#' Standardises a free-text or numeric comorbidity column to a consistent
#' coding for use as a covariate in \code{fit_mortality_rr_logistic()}.
#'
#' Three encoding strategies are applied in order:
#' \enumerate{
#'   \item \strong{Numeric} (Charlson / Elixhauser index already present):
#'         used as-is.
#'   \item \strong{Binary text} (\code{"present"} / \code{"none"}, etc.):
#'         recoded to integer 0 / 1.
#'   \item \strong{Ordinal text} (\code{"none"} / \code{"mild"} /
#'         \code{"moderate"} / \code{"severe"}): recoded to an ordered factor.
#' }
#' Missing / unknown values are set to \code{NA} in all cases.
#'
#' @param data Patient-level data frame.
#' @param comorbidity_col Character. Default \code{"comorbidities"}.
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#'
#' @return \code{data} with column \code{comorbidity_encoded} added.
#'   Attribute \code{"comorbidity_encoding"} records the strategy used:
#'   \code{"numeric"}, \code{"binary"}, \code{"ordinal"}, or \code{"absent"}.
#' @keywords internal
.encode_comorbidity_mortality <- function(
  data,
  comorbidity_col = "comorbidities",
  patient_id_col = "PatientInformation_id"
) {
  if (!comorbidity_col %in% names(data)) {
    message(sprintf(
      "[Comorbidity] Column '%s' not found — covariate will be omitted.",
      comorbidity_col
    ))
    data$comorbidity_encoded <- NA_real_
    attr(data, "comorbidity_encoding") <- "absent"
    return(data)
  }

  raw <- data[[comorbidity_col]]

  # Strategy 1: already a numeric index (Charlson / Elixhauser)
  if (is.numeric(raw)) {
    data$comorbidity_encoded <- raw
    attr(data, "comorbidity_encoding") <- "numeric"
    message(sprintf(
      "[Comorbidity] Numeric index detected (range %.1f–%.1f). Used as-is.",
      min(raw, na.rm = TRUE), max(raw, na.rm = TRUE)
    ))
    return(data)
  }

  # Normalise to uppercase for pattern matching
  raw_up <- stringr::str_to_upper(stringr::str_trim(as.character(raw)))
  # Sentinel missing values
  na_vals <- c(
    "NA", "NULL", "", "UNKNOWN", "NOT KNOWN", "NOT_KNOWN",
    "MISSING", "N/A", "NONE RECORDED"
  )
  raw_up[raw_up %in% na_vals | is.na(raw)] <- NA_character_

  # Strategy 2: binary present / none
  binary_present <- c(
    "PRESENT", "YES", "1", "TRUE", "POSITIVE",
    "COMORBID", "COMORBIDITY PRESENT"
  )
  binary_none <- c(
    "NONE", "NO", "0", "FALSE", "ABSENT",
    "NO COMORBIDITY", "COMORBIDITY ABSENT"
  )
  non_na_vals <- raw_up[!is.na(raw_up)]

  if (length(non_na_vals) > 0L &&
    all(non_na_vals %in% c(binary_present, binary_none))) {
    data$comorbidity_encoded <- dplyr::case_when(
      raw_up %in% binary_present ~ 1L,
      raw_up %in% binary_none ~ 0L,
      TRUE ~ NA_integer_
    )
    n_present <- sum(data$comorbidity_encoded == 1L, na.rm = TRUE)
    n_none <- sum(data$comorbidity_encoded == 0L, na.rm = TRUE)
    n_miss <- sum(is.na(data$comorbidity_encoded))
    attr(data, "comorbidity_encoding") <- "binary"
    message(sprintf(
      "[Comorbidity] Binary encoding: present = %d | none = %d | NA = %d.",
      n_present, n_none, n_miss
    ))
    return(data)
  }

  # Strategy 3: ordinal none / mild / moderate / severe
  ord_none <- c("NONE", "NO COMORBIDITY", "ABSENT", "0")
  ord_mild <- c("MILD", "LOW", "MINOR", "MINIMAL", "1")
  ord_moderate <- c("MODERATE", "MEDIUM", "MODERATE COMORBIDITY", "2")
  ord_severe <- c(
    "SEVERE", "HIGH", "MAJOR", "HEAVY",
    "SIGNIFICANT", "MULTIPLE", "3"
  )

  data$comorbidity_encoded <- dplyr::case_when(
    raw_up %in% ord_none ~ "none",
    raw_up %in% ord_mild ~ "mild",
    raw_up %in% ord_moderate ~ "moderate",
    raw_up %in% ord_severe ~ "severe",
    TRUE ~ NA_character_
  )

  n_unrecog <- sum(
    !is.na(raw_up) &
      !(raw_up %in% c(ord_none, ord_mild, ord_moderate, ord_severe))
  )
  if (n_unrecog > 0L) {
    message(sprintf(
      "[Comorbidity] %d unrecognised value(s) set to NA. Examples: %s",
      n_unrecog,
      paste(head(unique(
        raw_up[!is.na(raw_up) &
          !(raw_up %in% c(
            ord_none, ord_mild,
            ord_moderate, ord_severe
          ))]
      ), 5L), collapse = ", ")
    ))
  }

  data$comorbidity_encoded <- factor(
    data$comorbidity_encoded,
    levels  = c("none", "mild", "moderate", "severe"),
    ordered = TRUE
  )
  attr(data, "comorbidity_encoding") <- "ordinal"
  message(sprintf(
    "[Comorbidity] Ordinal encoding: %s.",
    paste(
      paste0(
        levels(data$comorbidity_encoded), " = ",
        as.integer(table(data$comorbidity_encoded))
      ),
      collapse = " | "
    )
  ))
  return(data)
}


# ── Step M4 ───────────────────────────────────────────────────────────────────

#' Check Collinearity Between HAI and ICU Covariates
#'
#' Computes the phi (Pearson) correlation coefficient for the 2×2 contingency
#' table of HAI × ICU. When ICU is hospital-acquired, the two covariates can
#' be highly correlated, making regression coefficients unstable.
#'
#' Perfect separation (any 2×2 marginal = 0) is detected and reported
#' separately, as it guarantees model non-convergence.
#'
#' @param df Patient-level data frame with integer 0/1 columns.
#' @param hai_col Character. Default \code{"HAI"}.
#' @param icu_col Character. Default \code{"ICU"}.
#' @param phi_threshold Numeric. Absolute phi above which a warning is issued.
#'   Default \code{0.7}.
#'
#' @return Named list: \code{phi}, \code{tbl} (2×2 table),
#'   \code{warning_issued} (logical).
#' @keywords internal
.check_hai_icu_collinearity <- function(df,
                                        hai_col = "HAI",
                                        icu_col = "ICU",
                                        phi_threshold = 0.7) {
  if (!hai_col %in% names(df) || !icu_col %in% names(df)) {
    return(list(phi = NA_real_, tbl = NULL, warning_issued = FALSE))
  }

  h <- as.integer(df[[hai_col]])
  ic <- as.integer(df[[icu_col]])
  ok <- !is.na(h) & !is.na(ic)
  h <- h[ok]
  ic <- ic[ok]

  if (length(h) < 10L) {
    message("[Collinearity] Too few complete observations to assess HAI/ICU correlation.")
    return(list(phi = NA_real_, tbl = NULL, warning_issued = FALSE))
  }

  tbl <- table(HAI = h, ICU = ic)
  message("[Collinearity] HAI \u00d7 ICU cross-table:")
  print(tbl)

  # Perfect separation: any marginal is zero
  if (any(rowSums(tbl) == 0L) || any(colSums(tbl) == 0L)) {
    warning(
      "[Collinearity] Perfect separation detected between HAI and ICU: ",
      "at least one cell combination has zero patients. ",
      "Model will likely fail to converge. ",
      "Consider dropping one of these two covariates."
    )
    return(list(phi = NA_real_, tbl = tbl, warning_issued = TRUE))
  }

  phi <- tryCatch(
    stats::cor(h, ic, method = "pearson"),
    error = function(e) NA_real_
  )

  warn_issued <- FALSE
  if (!is.na(phi) && abs(phi) >= phi_threshold) {
    warning(sprintf(
      paste0(
        "[Collinearity] HAI and ICU are highly correlated ",
        "(phi = %.3f \u2265 threshold %.2f). ",
        "Regression coefficients may be unstable. ",
        "Consider removing one covariate or checking your data."
      ),
      phi, phi_threshold
    ))
    warn_issued <- TRUE
  } else if (!is.na(phi)) {
    message(sprintf(
      "[Collinearity] HAI/ICU phi = %.3f (below threshold %.2f \u2014 OK).",
      phi, phi_threshold
    ))
  }

  list(phi = phi, tbl = tbl, warning_issued = warn_issued)
}


# ── Main function ─────────────────────────────────────────────────────────────

#' Estimate Per-Class Mortality Odds Ratio via Mixed-Effects Logistic Regression
#'
#' For each pathogen k and antibiotic class c, fits:
#' \preformatted{
#'   logit(P(death_i)) = b0 + b1*Resistance_ci + b2*Age_i + b3*Sex_i
#'                     + b4*HAI_i + b5*ICU_i + b6*Comorbidity_i
#'                     + u_facility   [u_facility ~ N(0, sigma^2)]
#' }
#' and returns \code{OR_death_kc = exp(b1)} with 95\% Wald CI.
#'
#' \strong{Resistance} is classified at antibiotic class level (class = 1 if
#' any drug in the class is R), matching
#' \code{.build_class_resistance_wide()}.
#'
#' \strong{HAI/CAI} is derived by
#' \code{derive_infection_type_for_mortality()}, which prioritises explicit
#' labels over date-gap inference. Patients labelled "Not Known" are excluded
#' from each model (missing HAI covariate).
#'
#' \strong{ICU flag} is derived using the "ever in ICU" rule by
#' \code{.derive_icu_binary()}: a patient is flagged ICU = 1 if any row for
#' that admission records an ICU unit type.
#'
#' \strong{Comorbidity} is standardised by
#' \code{.encode_comorbidity_mortality()} (numeric Charlson, binary
#' present/none, or ordinal none/mild/moderate/severe).
#'
#' \strong{Pre-fitting checks} per pathogen × class:
#' \itemize{
#'   \item \code{min_n} patients with complete required covariates.
#'   \item \code{min_deaths} deaths (ensures outcome variation).
#'   \item Variation in resistance (both R and S patients present).
#'   \item HAI/ICU collinearity (phi coefficient; \code{phi_threshold}).
#' }
#'
#' @param data Data frame. Merged AMR dataset at isolate \eqn{\times}
#'   antibiotic level (one row per patient \eqn{\times} drug test).
#' @param patient_id_col Character. Default \code{"PatientInformation_id"}.
#' @param facility_col Character. Default \code{"center_name"}.
#' @param organism_col Character. Default \code{"organism_name"}.
#' @param syndrome_col Character. Default \code{"syndrome"}.
#' @param infection_type_col Character. Raw infection type column (HAI/CAI
#'   labels before derivation). Default \code{"type_of_infection"}.
#' @param antibiotic_class_col Character. Default \code{"antibiotic_class"}.
#' @param antibiotic_name_col Character. Default \code{"antibiotic_name"}.
#' @param antibiotic_value_col Character. Default \code{"antibiotic_value"}.
#' @param unit_type_col Character. Ward/ICU column. Default \code{"unit_type"}.
#' @param date_admission_col Character. Default \code{"date_of_admission"}.
#' @param date_culture_col Character. Default
#'   \code{"date_of_first_positive_culture"}.
#' @param final_outcome_col Character. Default \code{"final_outcome"}.
#' @param final_outcome_date_col Character. Death date column (used only for
#'   data-quality reporting). Default \code{"final_outcome_date"}.
#' @param age_col Character. Default \code{"Age"}.
#' @param sex_col Character. Default \code{"Gender"}.
#' @param comorbidity_col Character or \code{NULL}. Name of the comorbidity
#'   column to include as a covariate. Set to \code{NULL} (default) to omit
#'   comorbidity from the model entirely. Supply the column name (e.g.
#'   \code{"comorbidities"}) to include it — the column is then standardised
#'   automatically by \code{.encode_comorbidity_mortality()}.
#' @param death_value Character. Value in \code{final_outcome_col} that
#'   identifies a death event. Default \code{"Death"}.
#' @param syndrome_name Character or \code{NULL}. Filter to one syndrome
#'   before fitting. \code{NULL} = all syndromes. Default \code{NULL}.
#' @param organism_name Character vector or \code{NULL}. Pathogen(s) to
#'   process. \code{NULL} = all found in data. Default \code{NULL}.
#' @param hai_threshold_hours Numeric. Admission-to-culture gap threshold (h)
#'   used when HAI/CAI is inferred from dates. Default \code{48}.
#' @param icu_values Character vector. Unit-type values treated as ICU
#'   admission (case-insensitive). Default
#'   \code{c("ICU", "Intensive Care", "Critical Care", "PICU", "NICU")}.
#' @param phi_threshold Numeric. HAI/ICU phi above which a collinearity
#'   warning is issued. Default \code{0.7}.
#' @param min_n Integer. Minimum patients per model. Default \code{10L}.
#' @param min_deaths Integer. Minimum deaths per model. Default \code{5L}.
#'
#' @return Data frame with one row per (pathogen, antibiotic_class):
#'   \code{pathogen}, \code{antibiotic_class}, \code{OR_death},
#'   \code{CI_lower}, \code{CI_upper}, \code{n_patients}, \code{n_deaths},
#'   \code{convergence_warning}, \code{syndrome_scope},
#'   \code{comorbidity_encoding}.
#' @export
fit_mortality_rr_logistic <- function(
  data,
  patient_id_col = "PatientInformation_id",
  facility_col = "center_name",
  organism_col = "organism_name",
  syndrome_col = "syndrome",
  infection_type_col = "type_of_infection",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  unit_type_col = "unit_type",
  date_admission_col = "date_of_admission",
  date_culture_col = "date_of_first_positive_culture",
  final_outcome_col = "final_outcome",
  final_outcome_date_col = "final_outcome_date",
  age_col = "Age",
  sex_col = "Gender",
  comorbidity_col = NULL,
  death_value = "Death",
  syndrome_name = NULL,
  organism_name = NULL,
  hai_threshold_hours = 48,
  icu_values = c(
    "ICU", "Intensive Care", "Critical Care",
    "PICU", "NICU"
  ),
  phi_threshold = 0.7,
  min_n = 10L,
  min_deaths = 5L
) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop(
      "Package 'lme4' is required for fit_mortality_rr_logistic(). ",
      "Install with: install.packages('lme4')"
    )
  }

  # ── Column validation ──────────────────────────────────────────────────────
  required <- c(
    patient_id_col, facility_col, organism_col,
    infection_type_col, antibiotic_class_col,
    antibiotic_name_col, antibiotic_value_col,
    date_admission_col, date_culture_col,
    final_outcome_col, age_col, sex_col
  )
  if (!is.null(syndrome_name)) required <- c(required, syndrome_col)
  # syndrome used as covariate when no syndrome_name filter is applied
  use_syndrome_covariate <- is.null(syndrome_name) && syndrome_col %in% names(data)
  missing_req <- setdiff(required, names(data))
  if (length(missing_req) > 0L) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing_req, collapse = ", ")
    ))
  }

  # ── Step A: derive infection type (mortality-specific) ─────────────────────
  data <- derive_infection_type_for_mortality(
    data,
    infection_type_col     = infection_type_col,
    date_admission_col     = date_admission_col,
    date_culture_col       = date_culture_col,
    final_outcome_col      = final_outcome_col,
    final_outcome_date_col = final_outcome_date_col,
    death_value            = death_value,
    hai_threshold_hours    = hai_threshold_hours,
    patient_id_col         = patient_id_col
  )

  # ── Step B: optional syndrome filter ──────────────────────────────────────
  df <- data
  if (!is.null(syndrome_name)) {
    df <- dplyr::filter(df, .data[[syndrome_col]] == syndrome_name)
  }

  pathogens <- if (!is.null(organism_name)) {
    organism_name
  } else {
    sort(unique(df[[organism_col]]))
  }

  # ── Step C: ICU binary (patient level, computed once globally) ─────────────
  icu_tbl <- .derive_icu_binary(
    df,
    patient_id_col    = patient_id_col,
    unit_type_col     = unit_type_col,
    icu_values        = icu_values
  )
  use_icu <- unit_type_col %in% names(df) && !all(is.na(icu_tbl$ICU))

  # ── Step D: patient-level covariate table (one row per patient) ────────────
  # Variables that are the same for every drug row of a patient are taken
  # from the first row.  The binary death outcome and HAI flag are derived
  # here and then re-used for every class model.
  patient_covars <- df %>%
    dplyr::mutate(
      death = as.integer(
        stringr::str_to_upper(stringr::str_trim(
          as.character(.data[[final_outcome_col]])
        )) == stringr::str_to_upper(death_value)
      ),
      HAI = as.integer(infection_type_derived == "HAI")
    ) %>%
    dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      !!rlang::sym(patient_id_col),
      !!rlang::sym(facility_col),
      !!rlang::sym(organism_col),
      !!rlang::sym(age_col),
      !!rlang::sym(sex_col),
      death,
      HAI,
      infection_type_derived
    ) %>%
    dplyr::mutate(
      !!rlang::sym(age_col) := suppressWarnings(
        as.numeric(.data[[age_col]])
      ),
      !!rlang::sym(sex_col) := factor(.data[[sex_col]])
    )

  # ── Syndrome covariate (patient level) ────────────────────────────────────
  if (use_syndrome_covariate) {
    syndrome_raw <- df %>%
      dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        !!rlang::sym(patient_id_col),
        !!rlang::sym(syndrome_col)
      ) %>%
      dplyr::mutate(
        !!rlang::sym(syndrome_col) := factor(.data[[syndrome_col]])
      )
    patient_covars <- dplyr::left_join(
      patient_covars, syndrome_raw,
      by = patient_id_col
    )
    message(sprintf(
      "[Syndrome] '%s' added as covariate (%d levels: %s).",
      syndrome_col,
      nlevels(patient_covars[[syndrome_col]]),
      paste(levels(patient_covars[[syndrome_col]]), collapse = ", ")
    ))
  }

  # ── Step E: comorbidity encoding (patient level) ───────────────────────────
  use_comorbidity <- FALSE
  comorbidity_encoding <- "absent"

  if (!is.null(comorbidity_col) && comorbidity_col %in% names(df)) {
    comorbid_raw <- df %>%
      dplyr::group_by(!!rlang::sym(patient_id_col)) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        !!rlang::sym(patient_id_col),
        !!rlang::sym(comorbidity_col)
      )

    comorbid_enc <- .encode_comorbidity_mortality(
      comorbid_raw,
      comorbidity_col = comorbidity_col,
      patient_id_col  = patient_id_col
    )
    comorbidity_encoding <- attr(comorbid_enc, "comorbidity_encoding")

    patient_covars <- dplyr::left_join(
      patient_covars,
      dplyr::select(
        comorbid_enc,
        !!rlang::sym(patient_id_col),
        comorbidity_encoded
      ),
      by = patient_id_col
    )
    use_comorbidity <- comorbidity_encoding != "absent" &&
      !all(is.na(patient_covars$comorbidity_encoded))
  } else {
    patient_covars$comorbidity_encoded <- NA_real_
    if (!is.null(comorbidity_col)) {
      message(sprintf(
        "[Comorbidity] Column '%s' not found — covariate omitted.",
        comorbidity_col
      ))
    }
  }

  # Join ICU flag
  patient_covars <- dplyr::left_join(patient_covars, icu_tbl,
    by = patient_id_col
  )

  # ── Main loop: one model per pathogen × class ──────────────────────────────
  all_or <- list()

  for (path in pathogens) {
    # Patients for this pathogen
    path_ids <- df %>%
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(
          .data[[organism_col]]
        )) == stringr::str_to_lower(stringr::str_trim(path))
      ) %>%
      dplyr::distinct(!!rlang::sym(patient_id_col)) %>%
      dplyr::pull(!!rlang::sym(patient_id_col))

    if (length(path_ids) == 0L) {
      message(sprintf("'%s': no patients found, skipping.", path))
      next
    }

    # Class resistance wide matrix for this pathogen
    path_df_raw <- df %>%
      dplyr::filter(
        stringr::str_to_lower(stringr::str_trim(
          .data[[organism_col]]
        )) == stringr::str_to_lower(stringr::str_trim(path))
      )

    resist_wide <- .build_class_resistance_wide(
      path_df_raw,
      patient_id_col       = patient_id_col,
      antibiotic_class_col = antibiotic_class_col,
      antibiotic_name_col  = antibiotic_name_col,
      antibiotic_value_col = antibiotic_value_col
    )
    class_name_map <- attr(resist_wide, "class_name_map")
    class_safe <- setdiff(names(resist_wide), patient_id_col)

    # Join covariates + resistance; exclude "Not Known" HAI rows and
    # patients with NA in any required covariate
    model_data <- patient_covars %>%
      dplyr::filter(.data[[patient_id_col]] %in% path_ids) %>%
      dplyr::inner_join(resist_wide, by = patient_id_col) %>%
      dplyr::filter(
        infection_type_derived != "Not Known",
        !is.na(.data[[age_col]]),
        !is.na(.data[[sex_col]]),
        !is.na(HAI)
      )
    if (use_syndrome_covariate) {
      model_data <- dplyr::filter(
        model_data, !is.na(.data[[syndrome_col]])
      )
    }

    if (nrow(model_data) < min_n) {
      message(sprintf(
        "'%s': only %d patients with complete covariates (min_n=%d), skipping.",
        path, nrow(model_data), min_n
      ))
      next
    }

    n_deaths_path <- sum(model_data$death, na.rm = TRUE)
    if (n_deaths_path < min_deaths) {
      message(sprintf(
        "'%s': only %d death(s) in cohort (min_deaths=%d), skipping.",
        path, n_deaths_path, min_deaths
      ))
      next
    }

    # HAI/ICU collinearity check — once per pathogen before the class loop
    if (use_icu) {
      invisible(.check_hai_icu_collinearity(
        model_data,
        hai_col = "HAI", icu_col = "ICU",
        phi_threshold = phi_threshold
      ))
    }

    n_facilities <- length(unique(model_data[[facility_col]]))

    # ── Per-class model fit ────────────────────────────────────────────────
    for (cls in class_safe) {
      orig_class <- class_name_map[[cls]]

      sub <- model_data %>%
        dplyr::filter(!is.na(.data[[cls]]))

      # Check resistance variation
      if (length(unique(stats::na.omit(sub[[cls]]))) < 2L) {
        message(sprintf(
          "'%s' | class '%s': no variation in resistance (all R or all S), skipping.",
          path, orig_class
        ))
        next
      }

      if (nrow(sub) < min_n) {
        message(sprintf(
          "'%s' | class '%s': %d patients after resistance filter (min_n=%d), skipping.",
          path, orig_class, nrow(sub), min_n
        ))
        next
      }

      n_deaths_cls <- sum(sub$death, na.rm = TRUE)
      n_surv_cls <- sum(sub$death == 0L, na.rm = TRUE)
      if (n_deaths_cls < min_deaths || n_surv_cls < min_deaths) {
        message(sprintf(
          "'%s' | class '%s': deaths = %d, survivors = %d (min_deaths=%d), skipping.",
          path, orig_class, n_deaths_cls, n_surv_cls, min_deaths
        ))
        next
      }

      # Build formula: fixed effects added only when covariate is usable
      fixed_terms <- c(
        sprintf("`%s`", cls),
        sprintf("`%s`", age_col),
        sprintf("`%s`", sex_col),
        "HAI"
      )
      if (use_icu &&
        length(unique(stats::na.omit(as.character(sub$ICU)))) > 1L) {
        fixed_terms <- c(fixed_terms, "ICU")
      }

      if (use_comorbidity &&
        !all(is.na(sub$comorbidity_encoded)) &&
        length(unique(stats::na.omit(
          as.character(sub$comorbidity_encoded)
        ))) > 1L) {
        fixed_terms <- c(fixed_terms, "comorbidity_encoded")
      }

      if (use_syndrome_covariate &&
        length(unique(stats::na.omit(
          as.character(sub[[syndrome_col]])
        ))) > 1L) {
        fixed_terms <- c(fixed_terms, sprintf("`%s`", syndrome_col))
      }

      # ── EPV check: ≥10 deaths per fixed predictor ──────────────────
      n_fixed_preds <- length(fixed_terms)
      epv <- n_deaths_cls / n_fixed_preds
      if (epv < 10) {
        message(sprintf(
          "'%s' | class '%s': EPV = %.1f (deaths=%d / predictors=%d) < 10 — estimates may be unreliable.",
          path, orig_class, epv, n_deaths_cls, n_fixed_preds
        ))
      }

      # ── Complete separation check ───────────────────────────────────
      ct_sep <- table(sub[[cls]], sub$death)
      if (any(ct_sep == 0L)) {
        message(sprintf(
          "'%s' | class '%s': zero cell(s) in resistance × death table (possible complete separation), skipping.",
          path, orig_class
        ))
        next
      }

      # Random intercept only when >1 facility present
      re_term <- if (n_facilities > 1L) {
        sprintf("(1 | `%s`)", facility_col)
      } else {
        NULL
      }

      fmla <- stats::as.formula(sprintf(
        "death ~ %s%s",
        paste(fixed_terms, collapse = " + "),
        if (!is.null(re_term)) paste0(" + ", re_term) else ""
      ))

      # Fit model with tryCatch; capture convergence warnings
      conv_warn <- FALSE
      model <- tryCatch(
        {
          fit <- suppressWarnings(
            lme4::glmer(
              fmla,
              data = sub,
              family = stats::binomial,
              control = lme4::glmerControl(
                optimizer = "bobyqa",
                optCtrl   = list(maxfun = 2e5)
              )
            )
          )
          # Singular fit (random-effect variance collapsed to 0)
          if (lme4::isSingular(fit)) {
            message(sprintf(
              "'%s' | class '%s': singular fit (random-effect variance ≈ 0).",
              path, orig_class
            ))
            conv_warn <- TRUE
          }
          # Optimizer convergence warnings
          conv_msgs <- fit@optinfo$conv$lme4$messages
          if (!is.null(conv_msgs) && length(conv_msgs) > 0L) {
            message(sprintf(
              "'%s' | class '%s': convergence warning — %s",
              path, orig_class, paste(conv_msgs, collapse = "; ")
            ))
            conv_warn <- TRUE
          }
          fit
        },
        error = function(e) {
          message(sprintf(
            "'%s' | class '%s': model failed — %s",
            path, orig_class, conditionMessage(e)
          ))
          NULL
        }
      )

      if (is.null(model)) next

      # Extract OR and 95% Wald CI for the resistance term
      coefs <- lme4::fixef(model)
      vcmat <- stats::vcov(model)
      pat <- paste0(
        "^`?",
        gsub(".", "\\.", cls, fixed = TRUE),
        "`?$"
      )
      cname <- grep(pat, names(coefs), value = TRUE)

      if (length(cname) == 0L) {
        message(sprintf(
          "'%s' | class '%s': resistance term not found in model output, skipping.",
          path, orig_class
        ))
        next
      }

      beta <- coefs[[cname[1L]]]
      se <- sqrt(vcmat[cname[1L], cname[1L]])

      # ── VIF check for multicollinearity ────────────────────────────
      if (requireNamespace("car", quietly = TRUE)) {
        vif_vals <- tryCatch(car::vif(model), error = function(e) NULL)
        if (!is.null(vif_vals)) {
          vif_num <- if (is.matrix(vif_vals)) {
            vif_vals[, "GVIF"]
          } else {
            as.numeric(vif_vals)
          }
          high_vif <- names(vif_num)[vif_num > 10]
          if (length(high_vif) > 0) {
            message(sprintf(
              "'%s' | class '%s': high VIF (>10) for: %s — check collinearity.",
              path, orig_class, paste(high_vif, collapse = ", ")
            ))
          }
        }
      }

      all_or[[length(all_or) + 1L]] <- data.frame(
        pathogen = path,
        antibiotic_class = orig_class,
        OR_death = formatC(
          exp(beta),
          format = "f",
          digits = 4L
        ),
        CI_lower = formatC(
          exp(beta - 1.96 * se),
          format = "f",
          digits = 4L
        ),
        CI_upper = formatC(
          exp(beta + 1.96 * se),
          format = "f",
          digits = 4L
        ),
        n_patients = nrow(sub),
        n_deaths = n_deaths_cls,
        convergence_warning = conv_warn,
        syndrome_scope = if (!is.null(syndrome_name)) syndrome_name else if (use_syndrome_covariate) "covariate" else "all",
        comorbidity_encoding = comorbidity_encoding,
        stringsAsFactors = FALSE
      )
    } # end class loop

    n_cls_fitted <- sum(sapply(all_or, function(x) x$pathogen == path))
    message(sprintf(
      "'%s': %d class OR(s) estimated (syndrome=%s, comorbidity=%s).",
      path, n_cls_fitted,
      if (!is.null(syndrome_name)) syndrome_name else "all",
      comorbidity_encoding
    ))
  } # end pathogen loop

  if (length(all_or) == 0L) {
    warning(
      "No mortality OR values could be estimated. ",
      "Check data, filters, min_n, and min_deaths."
    )
    return(data.frame())
  }

  dplyr::bind_rows(all_or)
}


# ── YLL Associated ────────────────────────────────────────────────────────────

#' Load and Parse India Life Expectancy Lookup Table
#'
#' Reads the statewise life expectancy Excel file and returns a tidy data frame
#' with columns \code{sex}, \code{age_bin}, and \code{life_expectancy} for the
#' India column only.  Three sex sections are parsed: \code{"Combined"},
#' \code{"Male"}, and \code{"Female"}.
#'
#' @param le_path Character.  Full path to the life expectancy xlsx file.
#'
#' @return Data frame with columns \code{sex} (\code{"Combined"},
#'   \code{"Male"}, \code{"Female"}), \code{age_bin}, \code{life_expectancy}.
#' @keywords internal
load_india_life_expectancy <- function(le_path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to load the life expectancy table.")
  }
  if (!file.exists(le_path)) {
    stop(sprintf("Life expectancy file not found: %s", le_path))
  }

  raw <- readxl::read_excel(le_path, col_names = FALSE, .name_repair = "minimal")

  # Table layout (1-indexed rows) for life_expectancy_all.xlsx:
  #   Row  1       : title  ("LIFE EXPECTANCY 2019-2023 INDIA")
  #   Row  2       : header (Age category, India, state names …)
  #   Rows  3 – 21 : Combined, 19 age bins
  #   Row  22      : "Male" section header
  #   Rows 23 – 41 : Male, 19 age bins
  #   Row  42      : "Female" section header
  #   Rows 43 – 61 : Female, 19 age bins
  # India life expectancy is always in column 2.

  age_bins <- c(
    "0-1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30",
    "30-35", "35-40", "40-45", "45-50", "50-55", "55-60",
    "60-65", "65-70", "70-75", "75-80", "80-85", "85+"
  )

  india_col <- 2L
  parse_le <- function(rows) suppressWarnings(as.numeric(raw[[india_col]][rows]))

  rbind(
    data.frame(
      sex = "Combined", age_bin = age_bins,
      life_expectancy = parse_le(3:21), stringsAsFactors = FALSE
    ),
    data.frame(
      sex = "Male", age_bin = age_bins,
      life_expectancy = parse_le(23:41), stringsAsFactors = FALSE
    ),
    data.frame(
      sex = "Female", age_bin = age_bins,
      life_expectancy = parse_le(43:61), stringsAsFactors = FALSE
    )
  )
}


#' Compute YLL Associated with AMR (Patient-Level, Facility-Direct)
#'
#' Computes years of life lost (YLL) associated with AMR directly from
#' patient-level facility records, without requiring population-level
#' P_LK (syndrome fractions) or R_kd (resistance profile scalars).
#'
#' For every fatal patient with pathogen k (and optionally syndrome L) the
#' individual YLL contribution is:
#' \deqn{\text{YLL}_{r,k} = \text{LE}(\text{age\_bin}_r, \text{sex}_r)
#'   \times w_{r,k}}
#' where \eqn{w_{r,k}} is the polymicrobial death weight for patient r and
#' pathogen k (= 1 for monomicrobial, 0–1 for polymicrobial episodes).
#' Total YLL associated:
#' \deqn{\text{YLL}_{\text{associated}} = \sum_{r,k} \text{YLL}_{r,k}}
#'
#' \strong{Polymicrobial weights} are computed via
#' \code{flag_polymicrobial()} + \code{compute_polymicrobial_weight()}
#' from \code{weight.R}.  When \code{facility_col} is provided the weights
#' are derived per facility (reflecting local organism distributions), with
#' automatic fallback to globally-pooled proportions for facilities whose
#' monomicrobial reference pool is smaller than \code{min_mono_per_facility}.
#' If \code{date_culture_col} is \code{NULL} polymicrobial flagging is
#' skipped and all weights default to 1.
#'
#' @param data                  Data frame of patient-level microbiology records
#'   (one row per patient × antibiotic test or per patient × pathogen).
#' @param outcome_col           Character.  Final outcome column.
#' @param death_value           Character.  Value(s) indicating a fatal outcome.
#'   Default \code{"Death"}.  Pass a vector to match multiple labels
#'   (e.g. \code{c("Death","Died")}).
#' @param pathogen_col          Character.  Pathogen / organism column (k).
#' @param patient_col           Character.  Unique patient identifier column.
#' @param age_bin_col           Character.  Column containing GBD-standard age
#'   bin labels (e.g. \code{"0-1"}, \code{"1-5"}, …, \code{"85+"}).  Use
#'   \code{age_bin_map} to recode non-standard labels.
#' @param sex_col               Character.  Column containing patient sex.
#' @param facility_col          Character or \code{NULL}.  Facility identifier
#'   column.  When supplied, polymicrobial weights are computed per facility
#'   (with global fallback) and results include a \code{by_facility} breakdown.
#' @param syndrome_col          Character or \code{NULL}.  Syndrome column.
#'   When \code{NULL} all syndromes are pooled; when supplied the
#'   \code{by_syndrome} and \code{by_syndrome_pathogen} outputs are populated.
#' @param syndrome_name         Character or \code{NULL}.  If supplied, data
#'   are filtered to this syndrome before computation.  \code{NULL} = all.
#' @param date_culture_col      Character or \code{NULL}.  Culture date column
#'   used to define polymicrobial episodes (\code{flag_polymicrobial()}).
#'   Set \code{NULL} to skip polymicrobial flagging (all weights = 1).
#' @param specimen_col          Character or \code{NULL}.  Specimen type column
#'   used together with \code{date_culture_col} for episode definition.
#' @param poly_weight_method    Character.  Passed to
#'   \code{compute_polymicrobial_weight()}: \code{"monomicrobial_proportion"}
#'   (default), \code{"equal"}, or \code{"manual"}.
#' @param min_mono_per_facility Integer.  Minimum monomicrobial isolates a
#'   facility must have for per-facility weights.  Facilities below this
#'   threshold fall back to globally-pooled proportions.  Default \code{30L}.
#' @param gap_days              Numeric.  Maximum days between cultures
#'   considered the same episode.  Default \code{14}.
#' @param le_path               Character.  Path to the India life expectancy
#'   xlsx file.  Defaults to the bundled \code{inst/extdata} copy.
#' @param male_value            Character.  Value in \code{sex_col} for males.
#'   Default \code{"Male"}.
#' @param female_value          Character.  Value in \code{sex_col} for females.
#'   Default \code{"Female"}.  All other values use the combined LE.
#' @param age_bin_map           Named character vector remapping non-standard
#'   age bin labels to LE-table labels.  Default \code{c("<1" = "0-1")}.
#' @param stratify_by           Character vector or \code{NULL}.  Additional
#'   column names to include as stratification dimensions in the
#'   \code{stratified} output (e.g. \code{c("location", "Age_bin", "sex")}).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{total}}{Scalar: total YLL associated across all pathogens and
#'     facilities.}
#'   \item{\code{per_pathogen}}{Data frame: YLL summed per pathogen k, pooled
#'     across facilities.  Columns: \code{pathogen_col}, \code{n_patients},
#'     \code{YLL_associated_k}.}
#'   \item{\code{by_age_sex}}{Data frame: YLL by \code{age_bin_col} × sex.}
#'   \item{\code{by_pathogen_age_sex}}{Data frame: YLL by pathogen ×
#'     \code{age_bin_col} × sex.}
#'   \item{\code{by_facility}}{Data frame (only when \code{facility_col} is
#'     supplied): per-facility YLL, one row per facility × pathogen.}
#'   \item{\code{by_syndrome}}{Data frame (only when \code{syndrome_col} is
#'     supplied and \code{syndrome_name} is \code{NULL}): YLL by syndrome.}
#'   \item{\code{by_syndrome_pathogen}}{Data frame (only when
#'     \code{syndrome_col} supplied): YLL by syndrome × pathogen.}
#'   \item{\code{stratified}}{Data frame (only when \code{stratify_by} is
#'     supplied): YLL aggregated by the requested columns.}
#'   \item{\code{patient_data}}{The death-cohort data frame used for
#'     computation, with \code{polymicrobial_weight}, \code{life_expectancy},
#'     and \code{yll_contribution} columns attached.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' yll <- compute_yll_associated(
#'   data             = cohort_df,
#'   outcome_col      = "final_outcome",
#'   death_value      = "Death",
#'   pathogen_col     = "organism_name",
#'   patient_col      = "PatientInformation_id",
#'   age_bin_col      = "Age_bin",
#'   sex_col          = "gender",
#'   facility_col     = "center_name",
#'   syndrome_col     = "infectious_syndrome",
#'   date_culture_col = "culture_date",
#'   specimen_col     = "sample_type",
#'   stratify_by      = c("location", "infectious_syndrome")
#' )
#' yll$total
#' yll$per_pathogen
#' yll$by_age_sex
#' yll$stratified
#' }
compute_yll_associated <- function(
  data,
  outcome_col,
  death_value = "Death",
  pathogen_col,
  patient_col,
  age_bin_col,
  sex_col,
  facility_col = NULL,
  syndrome_col = NULL,
  syndrome_name = NULL,
  date_culture_col = NULL,
  specimen_col = NULL,
  poly_weight_method = "monomicrobial_proportion",
  min_mono_per_facility = 30L,
  gap_days = 14,
  le_path = here::here("inst", "extdata", "life_expectancy_all.xlsx"),
  male_value = "Male",
  female_value = "Female",
  age_bin_map = c("<1" = "0-1"),
  stratify_by = NULL
) {
  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(
    outcome_col, pathogen_col, patient_col,
    age_bin_col, sex_col
  )
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)
  if (!is.null(date_culture_col)) required_cols <- c(required_cols, date_culture_col)
  if (!is.null(specimen_col)) required_cols <- c(required_cols, specimen_col)
  if (!is.null(stratify_by)) required_cols <- c(required_cols, stratify_by)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Missing column(s) in data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.null(syndrome_name) && is.null(syndrome_col)) {
    stop("syndrome_col must be supplied when syndrome_name is specified.")
  }

  # ── Step 1: Load life expectancy lookup ───────────────────────────────────
  le_lookup <- load_india_life_expectancy(le_path)

  # ── Step 2: Polymicrobial weights (per facility; global fallback) ─────────
  # Computed on the FULL data before death filter so the monomicrobial
  # reference pool is not restricted to fatal cases only.


  .compute_poly_weights <- function(df_sub, label) {
    if (is.null(date_culture_col)) {
      df_sub$polymicrobial_weight <- 1.0
      df_sub$weight_method <- "none"
      df_sub$weight_confidence <- "n/a"
      return(df_sub)
    }

    df_flagged <- tryCatch(
      flag_polymicrobial(
        df_sub,
        patient_col  = patient_col,
        specimen_col = if (!is.null(specimen_col)) specimen_col else "sample_type",
        date_col     = date_culture_col,
        organism_col = pathogen_col,
        gap_days     = gap_days
      ),
      error = function(e) {
        message(sprintf(
          "  [%s] flag_polymicrobial error: %s — weights set to 1.",
          label, conditionMessage(e)
        ))
        df_sub$polymicrobial_weight <- 1
        df_sub$is_polymicrobial <- 0L
        df_sub$episode_id <- seq_len(nrow(df_sub))
        return(df_sub)
      }
    )

    # Ensure required columns exist
    if (!"is_polymicrobial" %in% names(df_flagged)) {
      message(sprintf("[%s] 'is_polymicrobial' missing — assuming monomicrobial.", label))
      df_flagged$is_polymicrobial <- 0L
    }

    if (!"episode_id" %in% names(df_flagged)) {
      df_flagged$episode_id <- seq_len(nrow(df_flagged))
    }

    n_mono <- sum(df_flagged$is_polymicrobial == 0L, na.rm = TRUE)

    use_method <- if (n_mono >= min_mono_per_facility) {
      poly_weight_method
    } else {
      "equal"
    }

    if (n_mono < min_mono_per_facility) {
      message(sprintf(
        "  [%s] n_mono=%d < threshold %d — falling back to equal weights.",
        label, n_mono, min_mono_per_facility
      ))
    }

    df_weighted <- compute_polymicrobial_weight(
      df_flagged,
      episode_col       = "episode_id",
      organism_col      = pathogen_col,
      polymicrobial_col = "is_polymicrobial",
      method            = use_method
    )

    return(df_weighted)
  }

  if (!is.null(facility_col)) {
    facilities <- sort(unique(data[[facility_col]]))
    message(sprintf(
      "Computing polymicrobial weights for %d facility/facilities...",
      length(facilities)
    ))
    data_weighted <- dplyr::bind_rows(lapply(facilities, function(fac) {
      df_fac <- data[data[[facility_col]] == fac, , drop = FALSE]
      .compute_poly_weights(df_fac, label = fac)
    }))
  } else {
    message("Computing polymicrobial weights (global)...")
    data_weighted <- .compute_poly_weights(data, label = "global")
  }

  # ── Step 3: Filter to death cohort ────────────────────────────────────────
  df <- data_weighted %>%
    dplyr::filter(
      .data[[outcome_col]] %in% death_value,
      !is.na(.data[[pathogen_col]])
    )

  if (!is.null(syndrome_name)) {
    df <- df %>% dplyr::filter(.data[[syndrome_col]] == syndrome_name)
  }

  if (nrow(df) == 0L) {
    stop("No fatal records found after applying filters. Check outcome_col, death_value, and syndrome_name.")
  }

  message(sprintf(
    "Death cohort: %d rows | %d patients | %d pathogens%s",
    nrow(df),
    dplyr::n_distinct(df[[patient_col]]),
    dplyr::n_distinct(df[[pathogen_col]]),
    if (!is.null(facility_col)) {
      sprintf(" | %d facilities", dplyr::n_distinct(df[[facility_col]]))
    } else {
      ""
    }
  ))

  # ── Step 4: Recode age bins; normalise sex ────────────────────────────────
  if (length(age_bin_map) > 0L) {
    df[[age_bin_col]] <- dplyr::recode(
      as.character(df[[age_bin_col]]),
      !!!age_bin_map
    )
  }

  df <- df %>%
    dplyr::mutate(
      .sex_norm = dplyr::case_when(
        .data[[sex_col]] == male_value ~ "Male",
        .data[[sex_col]] == female_value ~ "Female",
        TRUE ~ "Combined"
      )
    )

  # ── Step 5: Join life expectancy (age bin × sex) ──────────────────────────
  join_by <- stats::setNames(c("age_bin", "sex"), c(age_bin_col, ".sex_norm"))
  df <- dplyr::left_join(df, le_lookup, by = join_by)

  n_missing_le <- sum(is.na(df$life_expectancy))
  if (n_missing_le > 0L) {
    warning(sprintf(
      "%d row(s) had no life expectancy match (unrecognised age_bin or sex); YLL set to NA.",
      n_missing_le
    ))
  }

  # ── Step 6: YLL per patient-pathogen row ──────────────────────────────────
  # death_weight = polymicrobial_weight from weight.R
  #   mono patient  → weight = 1.0  (full LE attributed to this pathogen)
  #   poly patient  → weight = 0–1  (fractional LE per pathogen)
  df <- df %>%
    dplyr::mutate(
      death_weight     = dplyr::coalesce(polymicrobial_weight, 1.0),
      yll_contribution = life_expectancy * death_weight
    )

  # ── Step 7: Aggregate ─────────────────────────────────────────────────────

  .agg <- function(df, grp_cols) {
    df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) %>%
      dplyr::summarise(
        n_patients       = dplyr::n_distinct(.data[[patient_col]]),
        YLL_associated_k = sum(yll_contribution, na.rm = TRUE),
        .groups          = "drop"
      )
  }

  # Per pathogen (+ facility if supplied)
  path_grp <- c(if (!is.null(facility_col)) facility_col, pathogen_col)
  by_path_fac <- .agg(df, path_grp)

  # Pool to per-pathogen (across facilities)
  per_pathogen <- if (!is.null(facility_col)) {
    by_path_fac %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        n_patients       = sum(n_patients, na.rm = TRUE),
        YLL_associated_k = sum(YLL_associated_k, na.rm = TRUE),
        .groups          = "drop"
      )
  } else {
    by_path_fac
  }

  # By age_bin × sex
  by_age_sex <- .agg(df, c(age_bin_col, ".sex_norm")) %>%
    dplyr::rename(sex = ".sex_norm")

  # By pathogen × age_bin × sex
  by_pathogen_age_sex <- .agg(df, c(pathogen_col, age_bin_col, ".sex_norm")) %>%
    dplyr::rename(sex = ".sex_norm")

  # By facility (if supplied)
  by_facility <- if (!is.null(facility_col)) {
    .agg(df, c(facility_col, pathogen_col))
  } else {
    NULL
  }

  # By syndrome (if syndrome_col supplied and not filtered to one syndrome)
  by_syndrome <- if (!is.null(syndrome_col) && is.null(syndrome_name)) {
    .agg(df, syndrome_col)
  } else {
    NULL
  }

  by_syndrome_pathogen <- if (!is.null(syndrome_col)) {
    .agg(df, c(syndrome_col, pathogen_col))
  } else {
    NULL
  }

  # User-defined stratification
  stratified <- if (!is.null(stratify_by)) {
    .agg(df, unique(c(stratify_by, pathogen_col)))
  } else {
    NULL
  }

  # ── Step 8: Total and summary message ─────────────────────────────────────
  YLL_total <- sum(per_pathogen$YLL_associated_k, na.rm = TRUE)

  message(sprintf(
    "YLL associated: total = %.2f years | %d pathogens | %d patients",
    YLL_total,
    nrow(per_pathogen),
    sum(per_pathogen$n_patients, na.rm = TRUE)
  ))

  # ── Return ─────────────────────────────────────────────────────────────────
  out <- list(
    total               = YLL_total,
    per_pathogen        = per_pathogen,
    by_age_sex          = by_age_sex,
    by_pathogen_age_sex = by_pathogen_age_sex,
    patient_data        = df
  )
  if (!is.null(facility_col)) out$by_facility <- by_facility
  if (!is.null(syndrome_col)) out$by_syndrome <- by_syndrome
  if (!is.null(syndrome_col)) out$by_syndrome_pathogen <- by_syndrome_pathogen
  if (!is.null(stratify_by)) out$stratified <- stratified

  out
}


# ── YLL Attributable ──────────────────────────────────────────────────────────

#' Compute YLL Attributable to AMR
#'
#' Takes the per-patient YLL data produced by \code{compute_yll_associated()}
#' (which already contains \code{life_expectancy}, \code{death_weight}, and
#' \code{yll_contribution}) and multiplies by the mortality PAF from
#' \code{compute_paf_rr_mortality()} to produce AMR-attributable YLL.
#'
#' \strong{Two PAF modes:}
#' \describe{
#'   \item{PAF_k scalar mode (default)}{When \code{resistance_profile_col} is
#'     \code{NULL}, the overall mortality PAF per pathogen k (\code{PAF_k_mort})
#'     from \code{paf_mort} is used as a scalar multiplier.}
#'   \item{Per-profile mode}{When \code{resistance_profile_col} is supplied,
#'     each patient row is matched to its resistance profile delta and the
#'     profile-specific \code{PAF_mortality} is applied.  Patients whose
#'     profile does not appear in \code{paf_mort} receive \code{NA} and a
#'     warning is issued.}
#' }
#'
#' \deqn{\text{YLL}^{\text{attr}}_{i,k} =
#'   \text{yll\_contribution}_{i,k} \times \text{PAF}_{k(,\delta)}}
#'
#' @param yll_patient_data Data frame.  The \code{$patient_data} element from
#'   \code{compute_yll_associated()}.  Must contain \code{yll_contribution},
#'   \code{life_expectancy}, \code{death_weight}, plus \code{pathogen_col},
#'   \code{patient_col}, \code{age_bin_col}, and \code{sex_col}.
#' @param paf_mort Named list returned by \code{compute_paf_rr_mortality()}.
#'   Each entry corresponds to one pathogen k and must contain
#'   \code{$PAF_k_mort} (scalar) and, for per-profile mode,
#'   \code{$per_profile} (data frame with \code{profile_col} and
#'   \code{PAF_mortality}).
#' @param pathogen_col  Character.  Pathogen column in \code{yll_patient_data}.
#' @param patient_col   Character.  Patient identifier column.
#' @param age_bin_col   Character.  Age bin column.
#' @param sex_col       Character.  Normalised sex column in
#'   \code{yll_patient_data}.  \code{compute_yll_associated()} stores the
#'   normalised values in \code{".sex_norm"}; pass that here.
#'   Default \code{".sex_norm"}.
#' @param resistance_profile_col Character or \code{NULL}.  Column in
#'   \code{yll_patient_data} containing the resistance profile delta label.
#'   When \code{NULL} (default), scalar \code{PAF_k_mort} is used.
#' @param profile_col   Character.  Column name in \code{paf_mort$per_profile}
#'   holding the profile label.  Default \code{"profile"}.
#' @param facility_col  Character or \code{NULL}.  Facility column.
#' @param syndrome_col  Character or \code{NULL}.  Syndrome column.
#' @param stratify_by   Character vector or \code{NULL}.  Additional
#'   stratification columns (always includes \code{pathogen_col}).
#'
#' @return A list:
#'   \describe{
#'     \item{total}{Scalar: total AMR-attributable YLL.}
#'     \item{per_pathogen}{One row per pathogen: \code{pathogen_col},
#'       \code{n_patients}, \code{YLL_associated_k}, \code{PAF_k_mort},
#'       \code{YLL_attributable_k}.}
#'     \item{by_age_sex}{YLL attributable stratified by age_bin × sex.}
#'     \item{by_pathogen_age_sex}{YLL attributable by pathogen × age_bin × sex.}
#'     \item{by_facility}{Per-facility breakdown (when \code{facility_col}
#'       supplied).}
#'     \item{by_syndrome}{Per-syndrome breakdown (when \code{syndrome_col}
#'       supplied).}
#'     \item{by_syndrome_pathogen}{Per-syndrome × pathogen (when
#'       \code{syndrome_col} supplied).}
#'     \item{stratified}{User-defined stratification (when \code{stratify_by}
#'       supplied).}
#'     \item{patient_data}{Patient-level data augmented with \code{PAF_kd}
#'       and \code{YLL_attributable_contribution}.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' yll_assoc <- compute_yll_associated(data = cohort, ...)
#' mort_or <- fit_mortality_rr_logistic(data = rr_data, ...)
#' profiles_or <- assign_rr_to_profiles(profiles_out,
#'   rr_table = mort_or,
#'   rr_col = "OR_death"
#' )
#' paf_mort <- compute_paf_rr_mortality(profiles_or)
#'
#' yll_attr <- compute_yll_attributable(
#'   yll_patient_data = yll_assoc$patient_data,
#'   paf_mort         = paf_mort,
#'   pathogen_col     = "organism_name",
#'   patient_col      = "PatientInformation_id",
#'   age_bin_col      = "Age_bin",
#'   sex_col          = ".sex_norm",
#'   facility_col     = "center_name",
#'   syndrome_col     = "infectious_syndrome",
#'   stratify_by      = c("location", "infectious_syndrome")
#' )
#' }
compute_yll_attributable <- function(
  yll_patient_data,
  paf_mort,
  pathogen_col,
  patient_col,
  age_bin_col,
  sex_col = ".sex_norm",
  resistance_profile_col = NULL,
  profile_col = "profile",
  facility_col = NULL,
  syndrome_col = NULL,
  stratify_by = NULL
) {
  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(
    "yll_contribution", "life_expectancy", "death_weight",
    pathogen_col, patient_col, age_bin_col, sex_col
  )
  if (!is.null(facility_col)) required_cols <- c(required_cols, facility_col)
  if (!is.null(syndrome_col)) required_cols <- c(required_cols, syndrome_col)
  if (!is.null(stratify_by)) required_cols <- c(required_cols, stratify_by)
  if (!is.null(resistance_profile_col)) required_cols <- c(required_cols, resistance_profile_col)

  missing_cols <- setdiff(required_cols, names(yll_patient_data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "Missing column(s) in yll_patient_data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.list(paf_mort) || length(paf_mort) == 0L) {
    stop("paf_mort must be the non-empty named list from compute_paf_rr_mortality().")
  }

  use_profile_mode <- !is.null(resistance_profile_col)
  message(sprintf(
    "PAF mode: %s",
    if (use_profile_mode) "per-profile (PAF_kd)" else "PAF_k scalar"
  ))

  df <- yll_patient_data

  # ── Step 1: Build PAF lookup and join ─────────────────────────────────────

  # Always build scalar PAF_k lookup (used for per_pathogen output table).
  paf_scalar_df <- do.call(rbind, lapply(names(paf_mort), function(k) {
    data.frame(
      .pathogen = k,
      PAF_k_mort = paf_mort[[k]]$PAF_k_mort,
      stringsAsFactors = FALSE
    )
  }))
  names(paf_scalar_df)[1L] <- pathogen_col

  if (!use_profile_mode) {
    # Warn on pathogen mismatches
    pats_data <- unique(df[[pathogen_col]])
    pats_paf <- names(paf_mort)
    missing_paf <- setdiff(pats_data, pats_paf)
    missing_data <- setdiff(pats_paf, pats_data)
    if (length(missing_paf) > 0L) {
      warning(sprintf(
        "Pathogen(s) in yll_patient_data have no PAF in paf_mort: %s. YLL_attributable set to NA.",
        paste(missing_paf, collapse = ", ")
      ))
    }
    if (length(missing_data) > 0L) {
      warning(sprintf(
        "Pathogen(s) in paf_mort have no records in yll_patient_data: %s.",
        paste(missing_data, collapse = ", ")
      ))
    }

    df <- dplyr::left_join(df, paf_scalar_df, by = pathogen_col)
    df <- df %>%
      dplyr::mutate(
        PAF_kd                        = PAF_k_mort,
        YLL_attributable_contribution = yll_contribution * PAF_kd
      )
  } else {
    # Per-profile: build (pathogen, profile_col) → PAF_mortality lookup
    paf_profile_df <- do.call(rbind, lapply(names(paf_mort), function(k) {
      pp <- paf_mort[[k]]$per_profile
      if (is.null(pp) || !is.data.frame(pp) || !profile_col %in% names(pp)) {
        return(NULL)
      }
      df_k <- pp[, c(profile_col, "PAF_mortality"), drop = FALSE]
      df_k[[pathogen_col]] <- k
      df_k
    }))
    if (is.null(paf_profile_df) || nrow(paf_profile_df) == 0L) {
      stop(paste0(
        "No per-profile PAF data found in paf_mort. ",
        "Check paf_mort$per_profile and the profile_col argument."
      ))
    }

    # Join: left(pathogen_col, resistance_profile_col) → right(pathogen_col, profile_col)
    join_keys <- stats::setNames(
      c(pathogen_col, profile_col),
      c(pathogen_col, resistance_profile_col)
    )
    df <- dplyr::left_join(df, paf_profile_df, by = join_keys)

    n_missing_paf <- sum(is.na(df$PAF_mortality))
    if (n_missing_paf > 0L) {
      warning(sprintf(
        "%d row(s) had no PAF_mortality match (unrecognised pathogen/profile combination); YLL_attributable set to NA.",
        n_missing_paf
      ))
    }

    df <- df %>%
      dplyr::mutate(
        PAF_kd                        = PAF_mortality,
        YLL_attributable_contribution = yll_contribution * PAF_kd
      )
  }

  # ── Step 2: Aggregate ─────────────────────────────────────────────────────

  .agg_attr <- function(df, grp_cols) {
    df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) %>%
      dplyr::summarise(
        n_patients         = dplyr::n_distinct(.data[[patient_col]]),
        YLL_associated_k   = sum(yll_contribution, na.rm = TRUE),
        YLL_attributable_k = sum(YLL_attributable_contribution, na.rm = TRUE),
        .groups            = "drop"
      )
  }

  # Rename sex_col → "sex" in output tables
  .rename_sex <- function(tbl) {
    if (sex_col != "sex" && sex_col %in% names(tbl)) {
      names(tbl)[names(tbl) == sex_col] <- "sex"
    }
    tbl
  }

  # Per pathogen (+ facility if supplied)
  path_grp <- c(if (!is.null(facility_col)) facility_col, pathogen_col)
  by_path_fac <- .agg_attr(df, path_grp)

  # Pool to per-pathogen across facilities; attach PAF_k_mort
  per_pathogen <- if (!is.null(facility_col)) {
    by_path_fac %>%
      dplyr::group_by(.data[[pathogen_col]]) %>%
      dplyr::summarise(
        n_patients         = sum(n_patients, na.rm = TRUE),
        YLL_associated_k   = sum(YLL_associated_k, na.rm = TRUE),
        YLL_attributable_k = sum(YLL_attributable_k, na.rm = TRUE),
        .groups            = "drop"
      )
  } else {
    by_path_fac
  }
  per_pathogen <- dplyr::left_join(per_pathogen, paf_scalar_df, by = pathogen_col)

  # By age_bin × sex
  by_age_sex <- .agg_attr(df, c(age_bin_col, sex_col)) %>%
    .rename_sex()

  # By pathogen × age_bin × sex
  by_pathogen_age_sex <- .agg_attr(df, c(pathogen_col, age_bin_col, sex_col)) %>%
    .rename_sex()

  # By facility
  by_facility <- if (!is.null(facility_col)) {
    .agg_attr(df, c(facility_col, pathogen_col))
  } else {
    NULL
  }

  # By syndrome
  by_syndrome <- if (!is.null(syndrome_col)) {
    .agg_attr(df, syndrome_col)
  } else {
    NULL
  }

  by_syndrome_pathogen <- if (!is.null(syndrome_col)) {
    .agg_attr(df, c(syndrome_col, pathogen_col))
  } else {
    NULL
  }

  # User-defined stratification
  stratified <- if (!is.null(stratify_by)) {
    .agg_attr(df, unique(c(stratify_by, pathogen_col)))
  } else {
    NULL
  }

  # ── Step 3: Total and summary ──────────────────────────────────────────────
  YLL_attributable_total <- sum(per_pathogen$YLL_attributable_k, na.rm = TRUE)
  YLL_associated_total <- sum(per_pathogen$YLL_associated_k, na.rm = TRUE)

  message(sprintf(
    "YLL attributable: total = %.2f years | associated = %.2f years | %d pathogens | %d patients",
    YLL_attributable_total,
    YLL_associated_total,
    nrow(per_pathogen),
    sum(per_pathogen$n_patients, na.rm = TRUE)
  ))

  # ── Return ─────────────────────────────────────────────────────────────────
  out <- list(
    total               = YLL_attributable_total,
    per_pathogen        = per_pathogen,
    by_age_sex          = by_age_sex,
    by_pathogen_age_sex = by_pathogen_age_sex,
    patient_data        = df
  )
  if (!is.null(facility_col)) out$by_facility <- by_facility
  if (!is.null(syndrome_col)) out$by_syndrome <- by_syndrome
  if (!is.null(syndrome_col)) out$by_syndrome_pathogen <- by_syndrome_pathogen
  if (!is.null(stratify_by)) out$stratified <- stratified

  out
}


# ── PAF Mortality ─────────────────────────────────────────────────────────────

#' Compute Mortality Population Attributable Fraction per Resistance Profile
#'
#' Computes PAF_kd_mortality for each pathogen k and resistance profile delta
#' using the GBD multi-exposure Levin formula, substituting the mortality
#' odds ratio (OR_death) from \code{fit_mortality_rr_logistic()} in place of
#' the LOS relative risk used by \code{compute_paf_los()}:
#'
#' \deqn{\text{PAF}_{kd,\text{mort}} =
#'   \frac{R'_{K\delta}\,(\text{OR}_{K\delta} - 1)}
#'        {1 + \sum_\delta R'_{K\delta}\,(\text{OR}_{K\delta} - 1)}}
#'
#' The all-susceptible profile carries OR = 1 and contributes 0.
#' The denominator equals
#' \eqn{E[\text{OR}_k] = \sum_\delta R'_{K\delta} \cdot \text{OR}_{K\delta}},
#' numerically identical to the denominator produced by \code{compute_paf_los()}.
#'
#' Overall mortality PAF for pathogen k:
#'
#' \deqn{\text{PAF}_{k,\text{mort}} = \sum_\delta \text{PAF}_{kd,\text{mort}}
#'   = \frac{\sum_\delta R'_{K\delta}(\text{OR}_{K\delta}-1)}
#'          {1 + \sum_\delta R'_{K\delta}(\text{OR}_{K\delta}-1)}}
#'
#' \strong{Usage pipeline:}
#' \preformatted{
#'   # 1. Fit mortality OR per class
#'   mort_or <- fit_mortality_rr_logistic(data, ...)
#'
#'   # 2. Assign OR to profiles via max rule (rr_col = "OR_death")
#'   profiles_with_or <- assign_rr_to_profiles(
#'       profiles_output,
#'       rr_table = mort_or,
#'       rr_col   = "OR_death"
#'   )
#'
#'   # 3. Compute per-profile and overall mortality PAF
#'   paf_mort <- compute_paf_rr_mortality(profiles_with_or)
#' }
#'
#' @param profiles_with_rr Named list from \code{assign_rr_to_profiles()}
#'   called with \code{rr_col = "OR_death"}.  Each element is a profile
#'   data frame for one pathogen.
#' @param probability_col  Character.  Profile probability column name.
#'   Default \code{"probability"}.
#' @param rr_profile_col   Character.  Profile-level OR column as produced by
#'   \code{assign_rr_to_profiles()}.  Default \code{"RR_LOS_profile"}.
#' @param profile_col      Character.  Profile label column.
#'   Default \code{"profile"}.
#'
#' @return Named list (one entry per pathogen) containing:
#'   \itemize{
#'     \item \code{per_profile}: profile data frame augmented with
#'       \code{numerator_mort} (= \eqn{R'_{K\delta}(\text{OR}_{K\delta}-1)}),
#'       \code{PAF_mortality} (= numerator / denominator), and
#'       \code{denominator_mort} (= \eqn{1 + \sum_\delta} numerator).
#'     \item \code{PAF_k_mort}: overall mortality PAF for pathogen k.
#'     \item \code{denominator_mort}: shared denominator \eqn{E[\text{OR}_k]}.
#'   }
#' @export
compute_paf_rr_mortality <- function(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
) {
  if (!is.list(profiles_with_rr)) {
    stop("profiles_with_rr must be the list returned by assign_rr_to_profiles().")
  }

  out <- list()

  for (path in names(profiles_with_rr)) {
    df <- profiles_with_rr[[path]]

    for (col in c(probability_col, rr_profile_col, profile_col)) {
      if (!col %in% names(df)) {
        stop(sprintf(
          "Column '%s' not found in profiles for '%s'.",
          col, path
        ))
      }
    }

    p <- df[[probability_col]] # R'_kd  (sums to 1)
    or <- df[[rr_profile_col]] # OR_death_kd
    numerator_vec <- p * (or - 1.0) # R'_kd * (OR_kd - 1)
    denom <- 1.0 + sum(numerator_vec, na.rm = TRUE)

    if (!is.finite(denom) || denom <= 0) {
      warning(sprintf(
        "'%s': PAF_mortality denominator = %.6g (must be > 0) — all mortality ORs may be <= 1 or NA; skipping.",
        path, denom
      ))
      next
    }

    paf_vec <- numerator_vec / denom

    df$numerator_mort <- round(numerator_vec, 6L)
    df$PAF_mortality <- round(paf_vec, 6L)
    df$denominator_mort <- round(denom, 6L)

    paf_k_mort <- sum(paf_vec, na.rm = TRUE)

    out[[path]] <- list(
      per_profile      = df,
      PAF_k_mort       = round(paf_k_mort, 6L),
      denominator_mort = round(denom, 6L)
    )

    message(sprintf(
      "'%s': PAF_k_mortality = %.4f | E[OR_death] (denominator) = %.4f | %d profiles.",
      path, paf_k_mort, denom, nrow(df)
    ))
  }

  return(out)
}
