# prep_ast_and_syndrome.R
# Mortality classification, MDR/XDR, resistance collapse/profile/selection,
# and polymicrobial flagging/weighting.
#
# Internal MDR/XDR reference helpers (not exported):
#   get_beta_lactam_hierarchy, get_magiorakos_thresholds,
#   get_antimicrobial_categories
#
# Moved to dedicated v2 files:
#   AST standardization    -> prep_ast.R
#   Contaminant flagging   -> prep_contaminants.R
#   HAI/CAI derivation     -> prep_hai_cai.R


# ---------------------------------------------------------------------------
# Internal helpers: MDR/XDR reference lookup tables (Magiorakos 2012)
# ---------------------------------------------------------------------------

get_beta_lactam_hierarchy <- function() {
  c(
    "Carbapenems",
    "Fourth-generation-cephalosporins",
    "Third-generation-cephalosporins",
    "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal",
    "Beta-lactam/beta-lactamase-inhibitor",
    "Aminopenicillins",
    "Penicillins"
  )
}

get_magiorakos_thresholds <- function() {
  tibble::tribble(
    ~organism_group,           ~mdr_threshold, ~xdr_threshold, ~total_categories,
    "Enterobacterales",         3,              "all_but_2",     9,
    "Pseudomonas aeruginosa",   3,              "all_but_2",    10,
    "Acinetobacter spp",        3,              "all_but_2",     9,
    "Staphylococcus aureus",    3,              "all_but_2",     9,
    "Enterococcus spp",         3,              "all_but_2",     7,
    "Streptococcus pneumoniae", 3,              "all_but_2",     5
  )
}

get_antimicrobial_categories <- function(organism_group = "Enterobacterales") {
  categories <- list(
    Enterobacterales = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (3rd gen)",
      "Cephalosporins (4th gen)",
      "Fluoroquinolones",
      "Monobactams",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Tigecycline"
    ),
    "Pseudomonas aeruginosa" = c(
      "Aminoglycosides",
      "Antipseudomonal carbapenems",
      "Antipseudomonal cephalosporins",
      "Antipseudomonal fluoroquinolones",
      "Antipseudomonal penicillins + beta-lactamase inhibitors",
      "Monobactams",
      "Phosphonic acids",
      "Polymyxins",
      "Ceftazidime-avibactam",
      "Ceftolozane-tazobactam"
    ),
    "Acinetobacter spp" = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (extended-spectrum)",
      "Fluoroquinolones",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Sulbactam",
      "Tetracyclines",
      "Tigecycline"
    ),
    "Staphylococcus aureus" = c(
      "Aminoglycosides",
      "Fluoroquinolones",
      "Folate pathway inhibitors",
      "Fusidic acid",
      "Glycopeptides",
      "Linezolid",
      "Mupirocin",
      "Rifampicin",
      "Tetracyclines"
    )
  )

  if (organism_group %in% names(categories)) {
    return(categories[[organism_group]])
  } else {
    warning(sprintf(
      "No category list for '%s', returning Enterobacterales default.",
      organism_group
    ))
    return(categories$Enterobacterales)
  }
}


#' Classify Infection-Related Mortality
#'
#' Determines if death was related to infection using date window logic.
#' Implements proxy logic when dates are missing.
#'
#' @param data Data frame
#' @param outcome_col Character. Outcome column. Default "final_outcome".
#' @param event_date_col Character. Event/culture date. Default "date_of_culture".
#' @param outcome_date_col Character. Outcome date. Default "date_of_final_outcome".
#' @param window Numeric. Days after event to classify death as infection-related.
#'   Default 14.
#'
#' @return Data frame with mortality_infection, mortality_method,
#'   mortality_confidence columns
#' @export
prep_classify_mortality <- function(data,
                                    outcome_col = "final_outcome",
                                    event_date_col = "date_of_culture",
                                    outcome_date_col = "date_of_final_outcome",
                                    window = 14) {
  # Check columns
  if (!outcome_col %in% names(data)) {
    stop(sprintf("Outcome column '%s' not found", outcome_col))
  }

  # Initialize
  data$mortality_infection <- "No"
  data$mortality_method <- NA_character_
  data$mortality_confidence <- NA_character_

  # Path A: Full date information available (HIGH CONFIDENCE)
  if (all(c(event_date_col, outcome_date_col) %in% names(data))) {
    data <- data %>%
      dplyr::mutate(
        gap_days = as.numeric(
          difftime(!!rlang::sym(outcome_date_col),
            !!rlang::sym(event_date_col),
            units = "days"
          )
        ),
        within_window = !is.na(gap_days) & gap_days >= 0 & gap_days <= window,
        mortality_infection = dplyr::case_when(
          !!rlang::sym(outcome_col) == "Died" & within_window ~ "Yes",
          !!rlang::sym(outcome_col) == "Died" & !within_window ~ "No",
          TRUE ~ "No"
        ),
        mortality_method = dplyr::if_else(
          !!rlang::sym(outcome_col) == "Died" & !is.na(gap_days),
          "date_calculated",
          NA_character_
        ),
        mortality_confidence = dplyr::case_when(
          mortality_method == "date_calculated" ~ "high",
          TRUE ~ NA_character_
        )
      )

    n_high_conf <- sum(data$mortality_confidence == "high", na.rm = TRUE)
    message(sprintf(
      "Classified mortality using dates (%d-day window): %d high confidence",
      window, n_high_conf
    ))
  } else {
    message("Date columns not available for mortality classification")
  }

  # Path B: PROXY - Only outcome available (LOW CONFIDENCE)
  data <- data %>%
    dplyr::mutate(
      mortality_infection = dplyr::case_when(
        # Keep high confidence classifications
        !is.na(mortality_confidence) ~ mortality_infection,

        # Proxy: Died but no dates -> mark as "Possible"
        !!rlang::sym(outcome_col) == "Died" ~ "Possible",
        TRUE ~ "No"
      ),
      mortality_method = dplyr::case_when(
        !is.na(mortality_method) ~ mortality_method,
        !!rlang::sym(outcome_col) == "Died" ~ "proxy_outcome_only",
        TRUE ~ NA_character_
      ),
      mortality_confidence = dplyr::case_when(
        !is.na(mortality_confidence) ~ mortality_confidence,
        mortality_method == "proxy_outcome_only" ~ "low",
        TRUE ~ NA_character_
      )
    )

  # Summary
  mortality_summary <- table(
    Method = data$mortality_method,
    Result = data$mortality_infection,
    useNA = "ifany"
  )

  message("\nMortality classification summary:")
  print(mortality_summary)

  n_proxy <- sum(data$mortality_method == "proxy_outcome_only", na.rm = TRUE)
  if (n_proxy > 0) {
    message(sprintf(
      "\n[!] Warning: %d deaths classified using PROXY (dates missing). Low confidence.",
      n_proxy
    ))
  }

  return(data)
}


#' Classify MDR and XDR
#'
#' Classifies isolates as MDR (Multidrug Resistant) and XDR (Extensively Drug
#' Resistant) in a single pass using Magiorakos 2012 criteria.
#'
#' MDR: resistant to at least one agent in three or more antimicrobial
#' categories. XDR: susceptible to at most two antimicrobial categories.
#'
#' Requires \code{class_result_event} column (run
#' \code{prep_collapse_class_level()} first, then rename \code{class_resistance}
#' to \code{class_result_event}).
#'
#' @param data Data frame with class-level resistance data.
#' @param definition Character. Classification criteria. Default "Magiorakos".
#' @param organism_group_col Character. Organism group column for
#'   pathogen-specific thresholds. Default "org_group".
#'
#' @return Data frame with mdr, mdr_confidence, mdr_method,
#'   n_resistant_categories, resistant_categories, xdr, xdr_confidence,
#'   and xdr_method columns added.
#' @export
#' @references
#' Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
prep_classify_mdr_xdr <- function(data,
                                   definition = "Magiorakos",
                                   organism_group_col = "org_group") {
  if (!"class_result_event" %in% names(data)) {
    stop("Must run prep_collapse_class_level() before MDR/XDR classification")
  }

  thresholds <- get_magiorakos_thresholds()
  n_total    <- dplyr::n_distinct(data$event_id)

  # --- MDR ---------------------------------------------------------------
  resistant_counts <- data %>%
    dplyr::filter(class_result_event == "R") %>%
    dplyr::group_by(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::summarise(
      n_resistant_categories = dplyr::n_distinct(antibiotic_class),
      resistant_categories   = paste(unique(antibiotic_class), collapse = "; "),
      .groups = "drop"
    )

  total_tested <- data %>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(
      n_total_categories = dplyr::n_distinct(antibiotic_class),
      .groups = "drop"
    )

  mdr_data <- resistant_counts %>%
    dplyr::left_join(total_tested, by = "event_id") %>%
    dplyr::left_join(thresholds, by = stats::setNames("organism_group", organism_group_col)) %>%
    dplyr::mutate(
      mdr_threshold  = dplyr::coalesce(mdr_threshold, 3L),
      mdr            = n_resistant_categories >= mdr_threshold,
      mdr_confidence = dplyr::case_when(
        n_total_categories >= 8 ~ "high",
        n_total_categories >= 5 ~ "medium",
        n_total_categories >= 3 ~ "low",
        TRUE ~ "insufficient_data"
      ),
      mdr_method = definition
    )

  data <- data %>%
    dplyr::left_join(
      mdr_data %>% dplyr::select(
        event_id, mdr, mdr_confidence, mdr_method,
        n_resistant_categories, resistant_categories
      ),
      by = "event_id"
    ) %>%
    dplyr::mutate(mdr = tidyr::replace_na(mdr, FALSE))

  n_mdr <- sum(data$mdr & data$mdr_confidence != "insufficient_data", na.rm = TRUE)
  message(sprintf("MDR classification (%s): %d/%d events (%.1f%%)",
                  definition, n_mdr, n_total, 100 * n_mdr / n_total))

  # --- XDR ---------------------------------------------------------------
  susceptible_counts <- data %>%
    dplyr::filter(class_result_event == "S") %>%
    dplyr::group_by(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::summarise(
      n_susceptible_categories = dplyr::n_distinct(antibiotic_class),
      .groups = "drop"
    )

  xdr_data <- data %>%
    dplyr::distinct(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::left_join(susceptible_counts, by = c("event_id", organism_group_col)) %>%
    dplyr::left_join(thresholds, by = stats::setNames("organism_group", organism_group_col)) %>%
    dplyr::mutate(
      n_susceptible_categories = tidyr::replace_na(n_susceptible_categories, 0),
      xdr            = n_susceptible_categories <= 2,
      xdr_confidence = dplyr::case_when(
        !is.na(total_categories) ~ "high",
        TRUE ~ "medium"
      ),
      xdr_method = definition
    )

  data <- data %>%
    dplyr::left_join(
      xdr_data %>% dplyr::select(event_id, xdr, xdr_confidence, xdr_method),
      by = "event_id"
    ) %>%
    dplyr::mutate(xdr = tidyr::replace_na(xdr, FALSE))

  n_xdr <- sum(data$xdr, na.rm = TRUE)
  message(sprintf("XDR classification (%s): %d/%d events (%.1f%%)",
                  definition, n_xdr, n_total, 100 * n_xdr / n_total))

  return(data)
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
#' class_level <- prep_collapse_class_level(data)
#' }
prep_collapse_class_level <- function(data,
                                      event_col = "event_id",
                                      organism_col = "organism_normalized",
                                      class_col = "antibiotic_class",
                                      susceptibility_col = "antibiotic_value",
                                      extra_cols = NULL) {
  required_cols <- c(event_col, organism_col, class_col, susceptibility_col)
  missing_cols  <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (!is.null(extra_cols)) {
    missing_extra <- setdiff(extra_cols, names(data))
    if (length(missing_extra) > 0) {
      stop(sprintf("Missing extra columns requested: %s", paste(missing_extra, collapse = ", ")))
    }
  }

  n_before   <- nrow(data)
  group_vars <- c(event_col, organism_col, class_col)

  message("Collapsing to antibiotic class level...")

  collapsed <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      class_resistance = dplyr::case_when(
        any(.data[[susceptibility_col]] == "R", na.rm = TRUE) ~ "R",
        any(.data[[susceptibility_col]] == "I", na.rm = TRUE) ~ "I",
        any(.data[[susceptibility_col]] == "S", na.rm = TRUE) ~ "S",
        TRUE ~ NA_character_
      ),
      n_drugs_in_class       = dplyr::n(),
      n_resistant            = sum(.data[[susceptibility_col]] == "R", na.rm = TRUE),
      pct_resistant_in_class = 100 * n_resistant / n_drugs_in_class,
      drugs_tested           = paste(
        sort(unique(.data[["antibiotic_normalized"]])),
        collapse = "; "
      ),
      .groups = "drop"
    )

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

  collapsed <- collapsed %>%
    dplyr::mutate(collapse_method = "class_any_R")

  n_after <- nrow(collapsed)
  message(sprintf("Collapsed: %d rows -> %d class-level rows", n_before, n_after))

  class_summary <- collapsed %>%
    dplyr::count(.data[[class_col]], class_resistance) %>%
    tidyr::pivot_wider(
      names_from  = class_resistance,
      values_from = n,
      values_fill = 0
    )

  message("\nClass-level resistance distribution:")
  print(class_summary)

  return(collapsed)
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
#' data_with_profile <- prep_create_resistance_profile(data)
#'
#' # Full S/R pattern
#' data_with_profile <- prep_create_resistance_profile(
#'   data,
#'   format = "full_pattern"
#' )
#'
#' # Resistant classes only
#' data_with_profile <- prep_create_resistance_profile(
#'   data,
#'   format = "class_summary"
#' )
#' }
prep_create_resistance_profile <- function(data,
                                           event_col = "event_id",
                                           antibiotic_col = "antibiotic_normalized",
                                           susceptibility_col = "antibiotic_value",
                                           format = "resistant_list",
                                           class_col = "antibiotic_class") {
  required_cols <- c(event_col, antibiotic_col, susceptibility_col)
  missing_cols  <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  if (format == "class_summary" && !class_col %in% names(data)) {
    stop(sprintf("Column '%s' required for format = 'class_summary'", class_col))
  }

  valid_formats <- c("resistant_list", "full_pattern", "class_summary")
  if (!format %in% valid_formats) {
    stop(sprintf("format must be one of: %s", paste(valid_formats, collapse = ", ")))
  }

  message(sprintf("Creating resistance profiles (format: %s)...", format))

  if (format == "resistant_list") {
    profiles <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R") %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(sort(unique(!!rlang::sym(antibiotic_col))), collapse = "; "),
        n_resistant        = dplyr::n_distinct(!!rlang::sym(antibiotic_col)),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(
        resistance_profile = dplyr::coalesce(resistance_profile, "None"),
        n_resistant        = dplyr::coalesce(n_resistant, 0L),
        profile_format     = "resistant_list"
      )

  } else if (format == "full_pattern") {
    profiles <- data %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::arrange(!!rlang::sym(antibiotic_col)) %>%
      dplyr::summarise(
        resistance_profile = paste(
          sprintf("%s:%s", !!rlang::sym(antibiotic_col), !!rlang::sym(susceptibility_col)),
          collapse = "; "
        ),
        n_resistant = sum(!!rlang::sym(susceptibility_col) == "R", na.rm = TRUE),
        n_tested    = dplyr::n(),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(profile_format = "full_pattern")

  } else if (format == "class_summary") {
    profiles <- data %>%
      dplyr::filter(!!rlang::sym(susceptibility_col) == "R") %>%
      dplyr::group_by(!!rlang::sym(event_col)) %>%
      dplyr::summarise(
        resistance_profile   = paste(sort(unique(!!rlang::sym(class_col))), collapse = "; "),
        n_resistant_classes  = dplyr::n_distinct(!!rlang::sym(class_col)),
        .groups = "drop"
      )

    data <- data %>%
      dplyr::left_join(profiles, by = event_col) %>%
      dplyr::mutate(
        resistance_profile  = dplyr::coalesce(resistance_profile, "None"),
        n_resistant_classes = dplyr::coalesce(n_resistant_classes, 0L),
        profile_format      = "class_summary"
      )
  }

  profile_freq <- data %>%
    dplyr::filter(!duplicated(!!rlang::sym(event_col))) %>%
    dplyr::count(resistance_profile, sort = TRUE) %>%
    utils::head(10)

  message("\nTop 10 resistance profiles:")
  print(profile_freq)

  return(data)
}
# flag_polymicrobial and compute_polymicrobial_weight have been moved to
# prep_polymicrobial.R. Backward-compatible aliases are defined there.
# compute_polymicrobial_weight has been moved to prep_polymicrobial.R.
