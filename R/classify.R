# classify.R
# Classification and enrichment functions for AMR preprocessing

#' Derive Age from Date of Birth
#'
#' Calculates age in years from date of birth and a reference date.
#' Implements proxy logic when Age is already present.
#'
#' @param data Data frame
#' @param dob_col Character. Name of DOB column. Default "DOB".
#' @param reference_date_col Character. Reference date column (usually culture date).
#'   Default "date_of_culture".
#' @param force Logical. If TRUE, recalculates even if Age present. Default FALSE.
#'
#' @return Data frame with Age column added/updated and Age_derived flag
#' @export
derive_age <- function(data, dob_col = "DOB",
                       reference_date_col = "date_of_culture",
                       force = FALSE) {
  # Check if Age already exists and force=FALSE
  if ("Age" %in% names(data) && !force) {
    n_present <- sum(!is.na(data$Age))
    if (n_present > 0) {
      message(sprintf(
        "Age already present for %d records. Use force=TRUE to recalculate.",
        n_present
      ))
      data$Age_derived <- FALSE
      return(data)
    }
  }

  # Check required columns
  if (!dob_col %in% names(data)) {
    stop(sprintf("DOB column '%s' not found", dob_col))
  }

  if (!reference_date_col %in% names(data)) {
    stop(sprintf("Reference date column '%s' not found", reference_date_col))
  }

  # Calculate age
  data$Age <- as.numeric(
    difftime(data[[reference_date_col]], data[[dob_col]], units = "days")
  ) / 365.25

  # Fix negative ages
  data$Age <- ifelse(data$Age < 0, 0, data$Age)

  # Add derivation flag
  data$Age_derived <- TRUE
  data$Age_method <- "calculated_from_dob"

  n_derived <- sum(!is.na(data$Age))
  message(sprintf("Derived Age for %d records from DOB", n_derived))

  return(data)
}


#' Assign Age Bins
#'
#' Categorizes age into bins for stratification.
#'
#' @param data Data frame with Age column
#' @param age_col Character. Name of age column. Default "Age".
#' @param bins Character or numeric vector. Either "GBD_standard", "pediatric",
#'   "geriatric", or custom bin breaks. Default "GBD_standard".
#'
#' @return Data frame with Age_bin column (factor)
#' @export
assign_age_bins <- function(data, age_col = "Age", bins = "GBD_standard") {
  if (!age_col %in% names(data)) {
    stop(sprintf("Age column '%s' not found", age_col))
  }

  # Get bin labels
  if (is.character(bins) && length(bins) == 1) {
    bin_labels <- get_age_bins(bins)
  } else {
    bin_labels <- bins
  }

  # Parse bin labels to get breaks
  breaks <- parse_age_bin_labels(bin_labels)

  # Assign bins
  data$Age_bin <- cut(
    data[[age_col]],
    breaks = breaks$breaks,
    labels = breaks$labels,
    right = FALSE,
    include.lowest = TRUE
  )

  n_binned <- sum(!is.na(data$Age_bin))
  n_unbinned <- sum(is.na(data$Age_bin) & !is.na(data[[age_col]]))

  message(sprintf(
    "Assigned age bins: %d binned, %d unbinned",
    n_binned, n_unbinned
  ))

  return(data)
}


#' Extract Genus from Organism Name
#'
#' Extracts bacterial genus (first word) from normalized organism name.
#'
#' @param data Data frame
#' @param organism_col Character. Organism column name.
#'   Default "organism_normalized".
#'
#' @return Data frame with org_genus column added
#' @export
extract_genus <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  data$org_genus <- stringr::str_extract(
    data[[organism_col]],
    "^[a-z]+"
  )

  n_extracted <- sum(!is.na(data$org_genus))
  n_unique <- length(unique(data$org_genus[!is.na(data$org_genus)]))

  message(sprintf(
    "Extracted genus: %d records, %d unique genera",
    n_extracted, n_unique
  ))

  return(data)
}


#' Extract Species from Organism Name
#'
#' Extracts bacterial species (second word) from normalized organism name.
#'
#' @param data Data frame
#' @param organism_col Character. Organism column name.
#'   Default "organism_normalized".
#'
#' @return Data frame with org_species column added
#' @export
extract_species <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  # Extract second word
  data$org_species <- stringr::str_extract(
    data[[organism_col]],
    "(?<=\\s)[a-z]+"
  )

  # Handle "spp" cases
  data$org_species <- ifelse(
    data$org_species == "spp",
    "species",
    data$org_species
  )

  n_extracted <- sum(!is.na(data$org_species))
  n_unique <- length(unique(data$org_species[!is.na(data$org_species)]))

  message(sprintf(
    "Extracted species: %d records, %d unique species",
    n_extracted, n_unique
  ))

  return(data)
}


#' Classify Organism Group
#'
#' Assigns organisms to taxonomic groups (Enterobacterales, Gram-positive, etc.)
#'
#' @param data Data frame
#' @param organism_col Character. Normalized organism column.
#'   Default "organism_normalized".
#'
#' @return Data frame with org_group column added
#' @export
classify_org_group <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  # Get taxonomy mapping
  taxonomy <- get_organism_taxonomy()

  # Join
  data <- data %>%
    dplyr::left_join(
      taxonomy,
      by = stats::setNames("organism_name", organism_col)
    )

  # Mark unmapped as "Other"
  data$org_group <- ifelse(
    is.na(data$org_group),
    "Other",
    data$org_group
  )

  # Summary
  group_counts <- table(data$org_group)
  message("Organism group distribution:")
  print(group_counts)

  return(data)
}


#' Classify Antibiotic to WHO Class
#'
#' Maps antibiotic names to WHO antibiotic classes.
#'
#' @param data Data frame
#' @param antibiotic_col Character. Normalized antibiotic column.
#'   Default "antibiotic_normalized".
#' @param who_table Data frame. WHO classification table. If NULL, uses
#'   built-in mapping.
#'
#' @return Data frame with antibiotic_class column added
#' @export
classify_antibiotic_class <- function(data,
                                      antibiotic_col = "antibiotic_normalized",
                                      who_table = NULL) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Antibiotic column '%s' not found", antibiotic_col))
  }

  # Use built-in if not provided
  if (is.null(who_table)) {
    # Placeholder - in production, load from inst/extdata/WHO_aware_class.csv
    message("Using built-in WHO class mapping")
    # For now, return with note
    data$antibiotic_class <- NA_character_
    data$class_source <- "needs_who_table"
    warning("WHO table not provided. Please supply or package with data.")
    return(data)
  }

  # Join with WHO table
  # Implementation when WHO table is provided
  data <- data %>%
    dplyr::left_join(
      who_table %>% dplyr::select(Antibiotic, Class),
      by = stats::setNames("Antibiotic", antibiotic_col)
    ) %>%
    dplyr::rename(antibiotic_class = Class)

  n_classified <- sum(!is.na(data$antibiotic_class))
  pct_classified <- 100 * n_classified / nrow(data)

  message(sprintf(
    "Classified antibiotics: %d/%d (%.1f%%)",
    n_classified, nrow(data), pct_classified
  ))

  return(data)
}


#' Classify AWaRe Category
#'
#' Assigns WHO AWaRe (Access, Watch, Reserve) categories to antibiotics.
#'
#' @param data Data frame
#' @param antibiotic_col Character. Antibiotic column.
#'   Default "antibiotic_normalized".
#' @param who_table Data frame. WHO AWaRe table. If NULL, uses built-in.
#'
#' @return Data frame with aware_category column added
#' @export
classify_aware <- function(data,
                           antibiotic_col = "antibiotic_normalized",
                           who_table = NULL) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Antibiotic column '%s' not found", antibiotic_col))
  }

  # Use built-in if not provided
  if (is.null(who_table)) {
    message("Using built-in AWaRe mapping")
    data$aware_category <- NA_character_
    data$aware_source <- "needs_who_table"
    warning("WHO AWaRe table not provided.")
    return(data)
  }

  # Join with WHO table
  data <- data %>%
    dplyr::left_join(
      who_table %>% dplyr::select(Antibiotic, Category),
      by = stats::setNames("Antibiotic", antibiotic_col)
    ) %>%
    dplyr::rename(aware_category = Category)

  n_classified <- sum(!is.na(data$aware_category))

  message(sprintf("Classified AWaRe: %d records", n_classified))

  aware_dist <- table(data$aware_category, useNA = "ifany")
  print(aware_dist)

  return(data)
}


#' Flag Contaminant Organisms
#'
#' Identifies likely contaminant organisms using multi-path logic.
#'
#' @param data Data frame
#' @param method Character. "auto" (try all methods), "device_based",
#'   "heuristic", "provided". Default "auto".
#' @param organism_col Character. Normalized organism column.
#' @param specimen_col Character. Specimen type column.
#'
#' @return Data frame with is_contaminant, contaminant_confidence,
#'   contaminant_method columns
#' @export
flag_contaminants <- function(data,
                              method = "auto",
                              organism_col = "organism_normalized",
                              specimen_col = "specimen_type") {
  # Check if required columns exist
  if (!organism_col %in% names(data)) {
    message(sprintf("[!] Column '%s' not found. Skipping contaminant flagging.", organism_col))
    data$is_contaminant <- FALSE
    data$contaminant_confidence <- "insufficient_data"
    data$contaminant_method <- "skipped"
    return(data)
  }

  if (!specimen_col %in% names(data)) {
    data[[specimen_col]] <- NA_character_
  }

  # Initialize
  data$is_contaminant <- FALSE
  data$contaminant_confidence <- "unknown"
  data$contaminant_method <- NA_character_

  # Path C: Use provided flag if exists
  if ("pathogen_contaminant" %in% names(data)) {
    data$is_contaminant <- data$pathogen_contaminant == 1
    data$contaminant_method <- "provided"
    data$contaminant_confidence <- "high"
    message("Using provided contaminant flags")
    return(data)
  }

  # Path A: Device-based (if device data available)
  if (all(c("device_inserted", "device_insertion_date", "date_of_culture") %in% names(data))) {
    data <- data %>%
      dplyr::mutate(
        days_since_device = as.numeric(
          difftime(date_of_culture, device_insertion_date, units = "days")
        )
      )

    # CoNS within 48h of line insertion = likely contaminant
    data <- data %>%
      dplyr::mutate(
        is_contaminant = dplyr::case_when(
          !!rlang::sym(organism_col) %in% c(
            "staphylococcus epidermidis",
            "staphylococcus haemolyticus"
          ) &
            !!rlang::sym(specimen_col) == "blood" &
            days_since_device <= 2 ~ TRUE,
          TRUE ~ is_contaminant
        ),
        contaminant_confidence = dplyr::if_else(
          !!rlang::sym(organism_col) %in% c(
            "staphylococcus epidermidis",
            "staphylococcus haemolyticus"
          ) &
            !!rlang::sym(specimen_col) == "blood",
          "high",
          contaminant_confidence
        ),
        contaminant_method = dplyr::if_else(
          is.na(contaminant_method),
          "device_based",
          contaminant_method
        )
      )

    n_device_based <- sum(data$contaminant_method == "device_based", na.rm = TRUE)
    message(sprintf("Device-based contaminants: %d", n_device_based))
  }

  # Path B: Heuristic (organism + specimen type)
  contam_blood <- get_contaminant_list(syndrome = "Bloodstream infections")
  contam_urine <- get_contaminant_list(syndrome = "Urinary tract infections / pyelonephritis")

  # Store previous contaminant flags
  prev_contaminant <- data$is_contaminant
  prev_method <- data$contaminant_method

  data <- data %>%
    dplyr::mutate(
      # Create temporary safe specimen column
      .specimen_lower = tolower(!!rlang::sym(specimen_col)),
      is_contaminant = dplyr::case_when(
        # Blood contaminants -- grepl handles "Blood-peripheral", "Blood-central catheter" etc.
        grepl("blood", .specimen_lower) &
          !!rlang::sym(organism_col) %in% contam_blood$names ~ TRUE,

        # Urine contaminants -- grepl handles "Urine culture", "Urine" etc.
        grepl("urine", .specimen_lower) &
          !!rlang::sym(organism_col) %in% contam_urine$names ~ TRUE,

        # Keep existing TRUE
        prev_contaminant == TRUE ~ TRUE,
        TRUE ~ FALSE
      ),
      contaminant_confidence = dplyr::case_when(
        is_contaminant & !is.na(prev_method) & prev_method == "device_based" ~ "high",
        is_contaminant ~ "medium",
        TRUE ~ "low"
      ),
      contaminant_method = dplyr::if_else(
        is.na(prev_method) | prev_method == "",
        "heuristic",
        prev_method
      ),
      # Remove temporary column
      .specimen_lower = NULL
    )

  # Summary
  n_contam <- sum(data$is_contaminant)
  n_high <- sum(data$contaminant_confidence == "high", na.rm = TRUE)
  n_medium <- sum(data$contaminant_confidence == "medium", na.rm = TRUE)

  message(sprintf(
    "Contaminants flagged: %d total (high: %d, medium: %d)",
    n_contam, n_high, n_medium
  ))

  return(data)
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
classify_mortality <- function(data,
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


#' Classify MDR (Multidrug Resistant)
#'
#' Classifies isolates as MDR using Magiorakos 2012 criteria.
#'
#' @param data Data frame
#' @param definition Character. "Magiorakos" or "WHO". Default "Magiorakos".
#' @param organism_group_col Character. Organism group column for
#'   pathogen-specific thresholds.
#'
#' @return Data frame with mdr, mdr_confidence, mdr_method,
#'   n_resistant_categories columns
#' @export
#' @references
#' Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
classify_mdr <- function(data,
                         definition = "Magiorakos",
                         organism_group_col = "org_group") {
  # Requires class-level resistance data
  if (!"class_result_event" %in% names(data)) {
    stop("Must run collapse_to_class_level() before MDR classification")
  }

  # Get thresholds
  thresholds <- get_magiorakos_thresholds()

  # Count resistant categories per event
  resistant_counts <- data %>%
    dplyr::filter(class_result_event == "R") %>%
    dplyr::group_by(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::summarise(
      n_resistant_categories = dplyr::n_distinct(antibiotic_class),
      resistant_categories = paste(unique(antibiotic_class), collapse = "; "),
      .groups = "drop"
    )

  # Total categories tested
  total_tested <- data %>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(
      n_total_categories = dplyr::n_distinct(antibiotic_class),
      .groups = "drop"
    )

  # Merge
  mdr_data <- resistant_counts %>%
    dplyr::left_join(total_tested, by = "event_id") %>%
    dplyr::left_join(thresholds, by = stats::setNames("organism_group", organism_group_col))

  # Apply MDR threshold
  mdr_data <- mdr_data %>%
    dplyr::mutate(
      mdr_threshold = dplyr::coalesce(mdr_threshold, 3), # Default to 3 if no match
      mdr = n_resistant_categories >= mdr_threshold,
      mdr_confidence = dplyr::case_when(
        n_total_categories >= 8 ~ "high",
        n_total_categories >= 5 ~ "medium",
        n_total_categories >= 3 ~ "low",
        TRUE ~ "insufficient_data"
      ),
      mdr_method = definition
    )

  # Join back to main data
  data <- data %>%
    dplyr::left_join(
      mdr_data %>% dplyr::select(
        event_id, mdr, mdr_confidence, mdr_method,
        n_resistant_categories, resistant_categories
      ),
      by = "event_id"
    ) %>%
    dplyr::mutate(mdr = tidyr::replace_na(mdr, FALSE))

  # Summary
  n_mdr <- sum(data$mdr & data$mdr_confidence != "insufficient_data", na.rm = TRUE)
  n_total <- dplyr::n_distinct(data$event_id)

  message(sprintf(
    "MDR classification (%s): %d/%d events (%.1f%%)",
    definition, n_mdr, n_total, 100 * n_mdr / n_total
  ))

  return(data)
}


#' Classify XDR (Extensively Drug Resistant)
#'
#' Classifies isolates as XDR using Magiorakos 2012 criteria.
#'
#' @param data Data frame
#' @param definition Character. "Magiorakos" or "WHO". Default "Magiorakos".
#' @param organism_group_col Character. Organism group column.
#'
#' @return Data frame with xdr, xdr_confidence, xdr_method columns
#' @export
#' @references
#' Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
classify_xdr <- function(data,
                         definition = "Magiorakos",
                         organism_group_col = "org_group") {
  # Requires MDR to be run first
  if (!"mdr" %in% names(data)) {
    message("Running MDR classification first...")
    data <- classify_mdr(data, definition, organism_group_col)
  }

  # Get thresholds
  thresholds <- get_magiorakos_thresholds()

  # Count susceptible categories
  susceptible_counts <- data %>%
    dplyr::filter(class_result_event == "S") %>%
    dplyr::group_by(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::summarise(
      n_susceptible_categories = dplyr::n_distinct(antibiotic_class),
      .groups = "drop"
    )

  # Get total categories per organism
  xdr_data <- data %>%
    dplyr::distinct(event_id, !!rlang::sym(organism_group_col)) %>%
    dplyr::left_join(susceptible_counts, by = c("event_id", organism_group_col)) %>%
    dplyr::left_join(thresholds, by = stats::setNames("organism_group", organism_group_col)) %>%
    dplyr::mutate(
      n_susceptible_categories = tidyr::replace_na(n_susceptible_categories, 0),
      # XDR: susceptible to <=2 categories
      xdr = n_susceptible_categories <= 2,
      xdr_confidence = dplyr::case_when(
        !is.na(total_categories) ~ "high", # Known pathogen
        TRUE ~ "medium"
      ),
      xdr_method = definition
    )

  # Join back
  data <- data %>%
    dplyr::left_join(
      xdr_data %>% dplyr::select(event_id, xdr, xdr_confidence, xdr_method),
      by = "event_id"
    ) %>%
    dplyr::mutate(xdr = tidyr::replace_na(xdr, FALSE))

  # Summary
  n_xdr <- sum(data$xdr, na.rm = TRUE)
  n_total <- dplyr::n_distinct(data$event_id)

  message(sprintf(
    "XDR classification (%s): %d/%d events (%.1f%%)",
    definition, n_xdr, n_total, 100 * n_xdr / n_total
  ))

  return(data)
}


#' Map Organism to RR Pathogen Category
#'
#' Maps normalized organism names to RR (Relative Risk) pathogen categories
#' used in burden estimation.
#'
#' @param data Data frame
#' @param organism_col Character. Normalized organism column.
#'
#' @return Data frame with rr_pathogen column added
#' @export
map_to_rr_pathogen <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Organism column '%s' not found", organism_col))
  }

  # Get mapping
  rr_map <- get_rr_pathogen_map()

  # Join
  data <- data %>%
    dplyr::left_join(
      rr_map,
      by = stats::setNames("organism_name", organism_col)
    )

  n_mapped <- sum(!is.na(data$rr_pathogen))
  pct_mapped <- 100 * n_mapped / nrow(data)

  message(sprintf(
    "Mapped to RR pathogen: %d/%d (%.1f%%)",
    n_mapped, nrow(data), pct_mapped
  ))

  unmapped_orgs <- unique(data[[organism_col]][is.na(data$rr_pathogen)])
  if (length(unmapped_orgs) > 0) {
    message("Unmapped organisms: ", paste(head(unmapped_orgs, 10), collapse = ", "))
  }

  return(data)
}


#' Map Antibiotic Class to RR Drug Category
#'
#' Maps WHO antibiotic classes to RR drug categories for burden estimation.
#'
#' @param data Data frame
#' @param class_col Character. Antibiotic class column. Default "antibiotic_class".
#'
#' @return Data frame with rr_drug column added
#' @export
map_class_to_rr <- function(data, class_col = "antibiotic_class") {
  if (!class_col %in% names(data)) {
    stop(sprintf("Class column '%s' not found", class_col))
  }

  # Get mapping
  class_map <- get_class_rr_map()

  # Join
  data <- data %>%
    dplyr::left_join(
      class_map,
      by = stats::setNames("Class", class_col)
    )

  n_mapped <- sum(!is.na(data$rr_drug))
  pct_mapped <- 100 * n_mapped / nrow(data)

  message(sprintf(
    "Mapped to RR drug: %d/%d (%.1f%%)",
    n_mapped, nrow(data), pct_mapped
  ))

  return(data)
}


#' Lookup Relative Risk Values
#'
#' Looks up RR values from table. Implements 3GC/4GC proxy logic.
#'
#' @param data Data frame with rr_pathogen and rr_drug columns
#' @param rr_table Data frame with RR values. If NULL, uses built-in.
#'
#' @return Data frame with rr_value column added
#' @export
lookup_rr <- function(data, rr_table = NULL) {
  if (!all(c("rr_pathogen", "rr_drug") %in% names(data))) {
    stop("Must have rr_pathogen and rr_drug columns. Run mapping functions first.")
  }

  if (is.null(rr_table)) {
    message("RR table not provided. Cannot lookup RR values.")
    data$rr_value <- NA_real_
    data$rr_source <- "no_table"
    return(data)
  }

  # Normalize for matching
  data_norm <- data %>%
    dplyr::mutate(
      rr_pathogen_norm = normalize_join(rr_pathogen),
      rr_drug_norm = normalize_join(rr_drug)
    )

  rr_norm <- rr_table %>%
    dplyr::mutate(
      rr_pathogen_norm = normalize_join(rr_pathogen),
      rr_drug_norm = normalize_join(rr_drug)
    )

  # Join
  data <- data_norm %>%
    dplyr::left_join(
      rr_norm %>% dplyr::select(rr_pathogen_norm, rr_drug_norm, RR),
      by = c("rr_pathogen_norm", "rr_drug_norm")
    ) %>%
    dplyr::rename(rr_value = RR) %>%
    dplyr::select(-rr_pathogen_norm, -rr_drug_norm)

  n_matched <- sum(!is.na(data$rr_value))
  pct_matched <- 100 * n_matched / nrow(data)

  message(sprintf(
    "RR values found: %d/%d (%.1f%%)",
    n_matched, nrow(data), pct_matched
  ))

  # TODO: Implement 3GC/4GC proxy logic if needed

  return(data)
}


#' Calculate Length of Stay
#'
#' Calculates hospital length of stay in days.
#'
#' @param data Data frame
#' @param admission_col Character. Admission date column.
#' @param outcome_col Character. Outcome/discharge date column.
#'
#' @return Data frame with Length_of_stay column
#' @export
calculate_los <- function(data,
                          admission_col = "date_of_admission",
                          outcome_col = "date_of_final_outcome") {
  if (!all(c(admission_col, outcome_col) %in% names(data))) {
    stop("Admission and outcome date columns not found")
  }

  data$Length_of_stay <- as.numeric(
    difftime(data[[outcome_col]], data[[admission_col]], units = "days")
  )

  # Flag suspicious values
  data$LOS_suspicious <- data$Length_of_stay < 0 | data$Length_of_stay > 365

  n_suspicious <- sum(data$LOS_suspicious, na.rm = TRUE)
  if (n_suspicious > 0) {
    warning(sprintf(
      "%d records have suspicious LOS (<0 or >365 days)",
      n_suspicious
    ))
  }

  message(sprintf(
    "Calculated LOS: median=%.1f days, mean=%.1f days",
    median(data$Length_of_stay, na.rm = TRUE),
    mean(data$Length_of_stay, na.rm = TRUE)
  ))

  return(data)
}


# ===== Helper Functions =====

#' Parse Age Bin Labels to Breaks
#' @keywords internal
parse_age_bin_labels <- function(labels) {
  breaks <- numeric(length(labels) + 1)
  clean_labels <- character(length(labels))

  for (i in seq_along(labels)) {
    label <- labels[i]

    if (grepl("\\+$", label)) {
      # Handle "85+" format
      lower <- as.numeric(gsub("\\+", "", label))
      breaks[i] <- lower
      breaks[i + 1] <- Inf
      clean_labels[i] <- label
    } else if (grepl("^<", label)) {
      # Handle "<1" format -- lower = -Inf to capture ages like -1
      upper <- as.numeric(gsub("^<", "", label))
      breaks[i] <- -Inf
      breaks[i + 1] <- upper
      clean_labels[i] <- label
    } else if (grepl("-", label)) {
      # Handle "1-5" format
      parts <- strsplit(label, "-")[[1]]
      lower <- as.numeric(parts[1])
      upper <- as.numeric(parts[2])
      breaks[i] <- lower
      breaks[i + 1] <- upper
      clean_labels[i] <- label
    } else {
      stop(sprintf("Cannot parse age bin label: %s", label))
    }
  }

  # Remove duplicate breaks
  breaks <- unique(breaks)

  return(list(breaks = breaks, labels = clean_labels))
}
