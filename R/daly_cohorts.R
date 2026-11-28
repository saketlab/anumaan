# daly_cohorts.R
# Cohort-level death and incidence calculations for DALY burden estimation

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


# D_J
#' Calculate Deaths by Underlying Cause (D_J)
#'
#' Computes D_J, the number of deaths for each underlying cause J from
#' population-level vital registration or mortality data. Optionally stratified
#' by grouping variables such as age group, sex, year, or location.
#'
#' This function requires population-level data. If only facility-level data
#' are available, use \code{daly_calc_deaths_by_syndrome()} with
#' \code{facility_data} instead, which directly counts deaths by syndrome.
#'
#' @param pop_data Data frame. Population-level vital registration or mortality
#'   data. Required -- this function does not accept facility-level data.
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
#' vital_reg <- data.frame(
#'   cause_of_death     = c("A41.9", "A41.9", "J18.9", "J18.9", "I21.9", "N39.0"),
#'   is_infection_death = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
#'   syndrome           = c("BSI", "BSI", "LRI", "LRI", NA, "UTI")
#' )
#' d_j <- daly_calc_deaths_by_cause(vital_reg)
#' d_j
daly_calc_deaths_by_cause <- function(pop_data,
                                      cause_col = "cause_of_death",
                                      deaths_col = NULL,
                                      groupby_cols = NULL) {
  if (missing(pop_data) || is.null(pop_data)) {
    stop(
      "pop_data is required for daly_calc_deaths_by_cause(). ",
      "If only facility-level data are available, use daly_calc_deaths_by_syndrome() ",
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


# S_J
#' Calculate Infection Fraction of Deaths by Cause (S_J)
#'
#' Computes S_J, the fraction of deaths for underlying cause J that were
#' related to infection. Values range from 0 (no infection involvement) to 1
#' (all deaths for cause J involved infection).
#'
#' This function requires population-level data where each death record carries
#' a binary flag indicating whether infection was involved. If only facility
#' data are available, use \code{daly_calc_deaths_by_syndrome()} with
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
#' vital_reg <- data.frame(
#'   cause_of_death     = c("A41.9", "A41.9", "J18.9", "J18.9", "I21.9", "N39.0"),
#'   is_infection_death = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
#'   syndrome           = c("BSI", "BSI", "LRI", "LRI", NA, "UTI")
#' )
#' s_j <- daly_calc_infection_fraction(vital_reg)
#' s_j
daly_calc_infection_fraction <- function(pop_data,
                                         cause_col = "cause_of_death",
                                         infection_flag_col = "is_infection_death",
                                         groupby_cols = NULL) {
  if (missing(pop_data) || is.null(pop_data)) {
    stop(
      "pop_data is required for daly_calc_infection_fraction(). ",
      "If only facility-level data are available, use daly_calc_deaths_by_syndrome() ",
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


# M_LJ
#' Calculate Infectious Syndrome Fraction (M_LJ)
#'
#' Computes M_LJ, the fraction of infection-related deaths for underlying cause
#' J that are attributed to infectious syndrome L. This distributes infection
#' deaths across clinical syndrome categories (e.g., bloodstream infection,
#' pneumonia, urinary tract infection).
#'
#' This function requires population-level data. If only facility data are
#' available, use \code{daly_calc_deaths_by_syndrome()} with \code{facility_data}
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
#' vital_reg <- data.frame(
#'   cause_of_death     = c("A41.9", "A41.9", "J18.9", "J18.9", "I21.9", "N39.0"),
#'   is_infection_death = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
#'   syndrome           = c("BSI", "BSI", "LRI", "LRI", NA, "UTI")
#' )
#' m_lj <- daly_calc_syndrome_fraction(vital_reg)
#' m_lj
daly_calc_syndrome_fraction <- function(pop_data,
                                        cause_col = "cause_of_death",
                                        syndrome_col = "syndrome",
                                        infection_flag_col = "is_infection_death",
                                        groupby_cols = NULL) {
  if (missing(pop_data) || is.null(pop_data)) {
    stop(
      "pop_data is required for daly_calc_syndrome_fraction(). ",
      "If only facility-level data are available, use daly_calc_deaths_by_syndrome() ",
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
      "No infection deaths found -- all values of '%s' are 0 or FALSE.",
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


# Deaths by Syndrome (D_L)
#' Calculate Deaths by Infectious Syndrome (D_L)
#'
#' Computes D_L, the number of deaths due to each infectious syndrome L.
#'
#' Three operating modes, evaluated in priority order:
#' \enumerate{
#'   \item \strong{Pre-computed components} (HIGH confidence): pass the outputs
#'     of \code{daly_calc_deaths_by_cause()}, \code{daly_calc_infection_fraction()},
#'     and \code{daly_calc_syndrome_fraction()} via \code{d_j}, \code{s_j}, and
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
#' @param d_j Data frame. Output of \code{daly_calc_deaths_by_cause()}.
#'   Supply together with \code{s_j} and \code{m_lj} for Mode 1.
#'   Default \code{NULL}.
#' @param s_j Data frame. Output of \code{daly_calc_infection_fraction()}.
#'   Default \code{NULL}.
#' @param m_lj Data frame. Output of \code{daly_calc_syndrome_fraction()}.
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
#' vital_reg <- data.frame(
#'   cause_of_death     = c("A41.9", "A41.9", "J18.9", "J18.9", "I21.9", "N39.0"),
#'   is_infection_death = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
#'   syndrome           = c("BSI", "BSI", "LRI", "LRI", NA, "UTI")
#' )
#' d_l <- daly_calc_deaths_by_syndrome(pop_data = vital_reg)
#' d_l
daly_calc_deaths_by_syndrome <- function(d_j = NULL,
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

  # Mode 1: Pre-computed components
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
        dplyr::filter(.name_key(!!rlang::sym(syndrome_col)) == .name_key(syndrome))
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

  # Mode 2: Raw population data -- compute all components internally
  if (has_pop_data) {
    message("Mode 2: Computing D_L from raw population data (all components computed internally)...")

    d_j_int <- daly_calc_deaths_by_cause(
      pop_data     = pop_data,
      cause_col    = cause_col,
      deaths_col   = deaths_col,
      groupby_cols = groupby_cols
    )
    s_j_int <- daly_calc_infection_fraction(
      pop_data           = pop_data,
      cause_col          = cause_col,
      infection_flag_col = infection_flag_col,
      groupby_cols       = groupby_cols
    )
    m_lj_int <- daly_calc_syndrome_fraction(
      pop_data           = pop_data,
      cause_col          = cause_col,
      syndrome_col       = syndrome_col,
      infection_flag_col = infection_flag_col,
      groupby_cols       = groupby_cols
    )

    # Recurse into Mode 1 with the computed components
    return(daly_calc_deaths_by_syndrome(
      d_j          = d_j_int,
      s_j          = s_j_int,
      m_lj         = m_lj_int,
      cause_col    = cause_col,
      syndrome_col = syndrome_col,
      syndrome     = syndrome,
      groupby_cols = groupby_cols
    ))
  }

  # Mode 3: Facility fallback -- count unique patients who died by syndrome
  message("Mode 3: No population data -- counting deaths by syndrome from facility data (LOW confidence)...")

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
      dplyr::filter(.name_key(!!rlang::sym(syndrome_col)) == .name_key(syndrome))
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

# Non-fatal incident cases by syndrome and hospital (direct count)
#' Count Non-Fatal Incident Cases by Syndrome and Hospital
#'
#' Counts the non-fatal incident cases of an infectious syndrome L at each
#' hospital h, \eqn{I_{h,L}}: the number of distinct \strong{discharged
#' admissions} with syndrome L. This is the direct count -- no CFR, no
#' \code{CR_L} adjustment and no pathogen weighting. It is the incidence used
#' in the YLD equations (see \code{daly_calc_yld_associated()}).
#'
#' Each admission is counted once per syndrome, however many isolates,
#' pathogens or antibiotic rows it has. Admissions with a missing or blank
#' syndrome are not counted under any syndrome; a message reports how many
#' discharged admissions that excludes. A warning is raised if an admission
#' carries more than one syndrome, because it is then counted under each.
#'
#' @param data Data frame with one or more rows per admission.
#' @param syndrome_col Character. Column holding the syndrome of each
#'   admission (e.g. \code{"chosen_syndrome"}).
#' @param syndrome_name Character vector or \code{NULL}. Syndrome(s) to
#'   count (matched ignoring case and spacing). \code{NULL} counts every
#'   syndrome present.
#' @param admission_col Character. Admission identifier column (e.g.
#'   \code{"readmission_id"}). Cases are distinct admissions, not rows or
#'   patients.
#' @param outcome_col Character. Final outcome column.
#' @param discharged_value Character. Value of \code{outcome_col} marking a
#'   discharged (non-fatal) admission. Default \code{"Discharged"}.
#' @param facility_col Character or \code{NULL}. Hospital column. When given,
#'   counts are returned per hospital; \code{NULL} pools all records.
#' @param facility_name Character or \code{NULL}. Restrict to one hospital.
#'
#' @return Data frame with \code{facility_col} (when given),
#'   \code{syndrome_col} and \code{N_NF_L} (distinct discharged admissions).
#'   Hospital x syndrome combinations with no discharged admission are not
#'   listed.
#' @export
#' @examples
#' adm <- data.frame(
#'   readmission_id = c("A1", "A1", "A2", "A3", "A4", "A5"),
#'   hospital = c("H1", "H1", "H1", "H1", "H2", "H2"),
#'   syndrome = c("BSI", "BSI", "BSI", "UTI", "BSI", ""),
#'   outcome = c("Discharged", "Discharged", "Death", "Discharged", "Discharged", "Discharged")
#' )
#' daly_count_non_fatal_incident_cases(adm, "syndrome", NULL, "readmission_id", "outcome",
#'   facility_col = "hospital"
#' )
daly_count_non_fatal_incident_cases <- function(data,
                                                syndrome_col,
                                                syndrome_name = NULL,
                                                admission_col,
                                                outcome_col,
                                                discharged_value = "Discharged",
                                                facility_col = NULL,
                                                facility_name = NULL) {
  # Input validation
  required_cols <- c(syndrome_col, admission_col, outcome_col, facility_col)
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

  # Step 1: discharged (non-fatal) records only
  df <- data[.name_key(data[[outcome_col]]) == .name_key(discharged_value) &
    !is.na(data[[outcome_col]]), , drop = FALSE]
  if (nrow(df) == 0) {
    warning(sprintf("No records where '%s' == '%s'.", outcome_col, discharged_value))
    return(data.frame())
  }

  # Step 2: optional single-hospital restriction
  if (!is.null(facility_name)) {
    df <- df[df[[facility_col]] == facility_name & !is.na(df[[facility_col]]), , drop = FALSE]
    if (nrow(df) == 0) {
      stop(sprintf("No discharged records found for facility '%s'.", facility_name))
    }
  }

  # Step 3: admissions without a syndrome are not counted under any syndrome
  no_syn <- is.na(df[[syndrome_col]]) | trimws(as.character(df[[syndrome_col]])) == ""
  n_no_syn <- length(setdiff(unique(df[[admission_col]][no_syn]), unique(df[[admission_col]][!no_syn])))
  if (n_no_syn > 0) {
    message(sprintf(
      "%d discharged admission(s) have no %s and are not counted under any syndrome.",
      n_no_syn, syndrome_col
    ))
  }
  df <- df[!no_syn, , drop = FALSE]

  # Step 4: optional syndrome filter
  if (!is.null(syndrome_name)) {
    df <- df[.name_key(df[[syndrome_col]]) %in% .name_key(syndrome_name), , drop = FALSE]
    if (nrow(df) == 0) {
      warning(sprintf(
        "No discharged records found for syndrome(s): %s",
        paste(syndrome_name, collapse = ", ")
      ))
      return(data.frame())
    }
  }

  # Syndromes are grouped ignoring case and spacing; the label reported is
  # the first spelling found in the data.
  df$.syn_key <- .name_key(df[[syndrome_col]])
  multi <- df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(admission_col, ".syn_key")))) %>%
    dplyr::count(.data[[admission_col]]) %>%
    dplyr::filter(.data[["n"]] > 1L)
  if (nrow(multi) > 0) {
    warning(sprintf(
      "%d admission(s) have more than one %s and are counted under each.",
      nrow(multi), syndrome_col
    ))
  }

  # Step 5: distinct discharged admissions per (hospital x) syndrome
  group_vars <- c(facility_col, ".syn_key")
  result <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      .syn_label = dplyr::first(as.character(.data[[syndrome_col]])),
      N_NF_L = dplyr::n_distinct(.data[[admission_col]]),
      .groups = "drop"
    ) %>%
    as.data.frame()
  result[[syndrome_col]] <- result$.syn_label
  result <- result[, c(facility_col, syndrome_col, "N_NF_L"), drop = FALSE]

  message(sprintf(
    "Non-fatal incident cases: %d discharged admission(s) across %d syndrome(s)%s.",
    sum(result$N_NF_L),
    dplyr::n_distinct(result[[syndrome_col]]),
    if (!is.null(facility_col)) sprintf(" and %d hospital(s)", dplyr::n_distinct(result[[facility_col]])) else ""
  ))

  result
}


# CR_L : CFR adjustment factor
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

daly_calc_cr_l <- function(data,
                           syndrome_col,
                           syndrome_name,
                           patient_col,
                           visit_type_col,
                           inpatient_value = "Inpatient",
                           outpatient_value = "Outpatient",
                           adjustment_ref,
                           facility_col = NULL,
                           facility_name = NULL) {
  # Input validation
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

  # Look up adjustment type
  adj_type <- adjustment_ref %>%
    dplyr::filter(.name_key(.data[["infectious_syndrome"]]) == .name_key(syndrome_name)) %>%
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
    message(sprintf("CR_L for '%s': adjustment type = None -> CR_L = 1.", syndrome_name))
    return(data.frame(
      syndrome        = syndrome_name,
      adjustment_type = "None",
      CR_L            = 1
    ))
  }

  # Filter to syndrome
  df <- data %>%
    dplyr::filter(.name_key(.data[[syndrome_col]]) == .name_key(syndrome_name))

  if (!is.null(facility_name)) {
    df <- df %>% dplyr::filter(.data[[facility_col]] == facility_name)
    if (nrow(df) == 0) {
      stop(sprintf("No records found for facility '%s'.", facility_name))
    }
  }

  fac_grp <- if (!is.null(facility_col)) facility_col else character(0)

  # Inpatient ratio
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

    # Outpatient to inpatient ratio
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


# Incidence (formula-based)
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
#' use \code{daly_count_non_fatal_incident_cases()} instead.
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
#'   \code{daly_calc_cr_l()}. Default \code{1} (no adjustment).
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

daly_calc_incidence_from_cfr <- function(deaths_L,
                                         cfr_lk_tbl,
                                         P_Lk_prime_tbl,
                                         CR_L = 1,
                                         pathogen_col = "pathogen",
                                         cfr_col = "CFR_LK",
                                         plk_col = "P_Lk_prime",
                                         facility_col = NULL,
                                         deaths_col = "deaths") {
  # Input validation
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

  # Join P'LK and CFR_LK
  join_cols <- if (!is.null(facility_col)) c(facility_col, pathogen_col) else pathogen_col

  merged <- dplyr::inner_join(P_Lk_prime_tbl, cfr_lk_tbl, by = join_cols)

  n_unmatched <- nrow(P_Lk_prime_tbl) - nrow(merged)
  if (n_unmatched > 0) {
    warning(sprintf(
      "%d pathogen(s) in P_Lk_prime_tbl have no matching CFR and are excluded from CFR_L.",
      n_unmatched
    ))
  }

  # CFR_L = sum_K( P'_LK x CFR_LK ) per facility
  grp <- if (!is.null(facility_col)) facility_col else character(0)

  cfr_L_tbl <- merged %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
    dplyr::summarise(
      CFR_L = sum(.data[[plk_col]] * .data[[cfr_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # Attach deaths and compute I_L
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
    "Incidence I_L computed: CFR_L = %.4f, CR_L = %.4f -> I_L = %.1f%s.",
    mean(result$CFR_L, na.rm = TRUE),
    mean(result$CR_L, na.rm = TRUE),
    mean(result$I_L, na.rm = TRUE),
    if (!is.null(facility_col)) sprintf(" (mean across %d facility/facilities)", nrow(result)) else ""
  ))
  return(result)
}
