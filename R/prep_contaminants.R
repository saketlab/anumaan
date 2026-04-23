# prep_contaminants.R
# Layer 6c: Contaminant flagging
#
# Functions moved here from prep_ast_and_syndrome.R:
#   - prep_get_contaminant_list
#   - prep_is_contaminant
#   - prep_flag_contaminants  (extended with context parameter)
#
# New functions: prep_filter_fungal, prep_flag_opportunists


#' Get Contaminant List from Reference File
#'
#' Loads contaminant organisms from \code{common_commensals.csv} and filters
#' by syndrome and/or specimen type.
#'
#' @param syndrome Character. Syndrome name to filter (e.g., "Bloodstream infections").
#' @param specimen_type Character. Specimen type to filter (e.g., "Blood culture").
#' @param return_all Logical. If TRUE, returns all contaminants across syndromes.
#'   Default FALSE.
#'
#' @return A list with \code{names} (character vector) and \code{patterns}
#'   (list of pattern information for flexible matching).
#' @export
prep_get_contaminant_list <- function(syndrome      = NULL,
                                       specimen_type = NULL,
                                       return_all    = FALSE) {
  commensals_path <- find_extdata_file("common_commensals.csv")

  if (commensals_path == "" || !file.exists(commensals_path)) {
    warning("common_commensals.csv not found in inst/extdata/. Returning empty list.")
    return(list(names = character(0), patterns = list()))
  }

  commensals_data <- readr::read_csv(commensals_path, show_col_types = FALSE)

  if (!is.null(syndrome) && !return_all)
    commensals_data <- dplyr::filter(commensals_data, Syndrome == syndrome)
  if (!is.null(specimen_type) && !return_all)
    commensals_data <- dplyr::filter(commensals_data, `Type of culture/specimen` == specimen_type)

  contaminants <- commensals_data %>%
    dplyr::pull(`Common commensals`) %>%
    stringr::str_split(";\\s*") %>%
    unlist() %>%
    stringr::str_trim() %>%
    unique()
  contaminants <- contaminants[contaminants != ""]

  contaminant_patterns <- lapply(contaminants, function(name) {
    name_lower <- tolower(name)
    parts      <- strsplit(name_lower, "\\s+")[[1]]
    patterns   <- c(name_lower)

    if (length(parts) >= 1) {
      genus              <- parts[1]
      is_genus_level     <- length(parts) == 1 ||
        grepl("^sp(p\\.?|ecies)$", parts[2], ignore.case = TRUE)
      genus_abbrev_base  <- paste0("^", substr(genus, 1, 1), "\\.")
      genus_abbrev_alone <- paste0("^", substr(genus, 1, 1), "\\.\\s")

      if (is_genus_level) patterns <- c(patterns, genus, genus_abbrev_alone)
      else                patterns <- c(patterns, genus_abbrev_alone)

      if (length(parts) >= 2) {
        species         <- parts[2]
        is_real_species <- !grepl("^sp(p\\.?|ecies)$", species, ignore.case = TRUE) &&
          !grepl("-", genus)
        if (is_real_species) patterns <- c(patterns, species)
        patterns <- c(patterns,
                      paste0(genus_abbrev_base, "\\s*", species),
                      paste(genus, species))
      }

      if (grepl("staphylococcus", genus)) {
        if (is_genus_level)
          patterns <- c(patterns, "staph", "coag.*neg", "coagulase.*negative")
        else
          patterns <- c(patterns, paste0("staph.*", parts[2]), "coag.*neg", "coagulase.*negative")
      }
      if (grepl("streptococcus",    genus)) patterns <- c(patterns, "strep")
      if (grepl("escherichia",      genus)) patterns <- c(patterns, "e\\.?\\s*coli")
      if (grepl("pseudomonas",      genus)) patterns <- c(patterns, "pseudo")
      if (grepl("klebsiella",       genus)) patterns <- c(patterns, "kleb")
      if (grepl("acinetobacter",    genus)) patterns <- c(patterns, "acin")
      if (grepl("enterococcus",     genus)) patterns <- c(patterns, "entero")
      if (grepl("corynebacterium",  genus)) patterns <- c(patterns, "coryno", "diphtheroids?")
      if (grepl("bacillus",         genus)) patterns <- c(patterns, "bacil")
      if (grepl("micrococcus",      genus)) patterns <- c(patterns, "micro")
      if (grepl("cutibacterium|propionibacterium", genus))
        patterns <- c(patterns, "propioni", "cuti", "p\\.?\\s*acnes")
      if (grepl("lactobacillus",    genus)) patterns <- c(patterns, "lacto")
    }

    list(original = name, patterns = unique(patterns))
  })

  list(names = contaminants, patterns = contaminant_patterns)
}


#' Check if Organism is a Contaminant
#'
#' Checks whether an organism name matches any known contaminant for a specific
#' syndrome or specimen type using flexible pattern matching.
#'
#' @param organism_name Character vector of organism name(s) to check.
#' @param syndrome Character. Optional syndrome filter.
#' @param specimen_type Character. Optional specimen type filter.
#'
#' @return Logical vector (same length as \code{organism_name}).
#' @export
prep_is_contaminant <- function(organism_name,
                                 syndrome      = NULL,
                                 specimen_type = NULL) {
  contaminant_data <- prep_get_contaminant_list(syndrome      = syndrome,
                                                 specimen_type = specimen_type)

  sapply(organism_name, function(org) {
    if (is.na(org) || org == "") return(FALSE)
    org_lower <- tolower(trimws(org))

    for (contam_info in contaminant_data$patterns)
      for (pattern in contam_info$patterns)
        if (grepl(pattern, org_lower, perl = TRUE)) return(TRUE)

    FALSE
  }, USE.NAMES = FALSE)
}


#' Flag Contaminant Organisms
#'
#' Identifies likely contaminant organisms using multi-path logic (provided flags,
#' device-based context, or heuristic organism+specimen matching).
#'
#' When \code{context = "icu_device"} (AIIMS ICU BSI data), CoNS in blood with
#' a central line present AND line in place > 2 days are classified as CLABSI
#' pathogens rather than contaminants.
#'
#' @param data Data frame.
#' @param method Character. "auto" (try all methods), "device_based",
#'   "heuristic", or "provided". Default "auto".
#' @param context Character. "general" or "icu_device". Default "general".
#' @param organism_col Character. Normalized organism column.
#'   Default "organism_normalized".
#' @param specimen_col Character. Specimen type column. Default "specimen_type".
#'
#' @return Data frame with \code{is_contaminant}, \code{contaminant_confidence},
#'   and \code{contaminant_method} columns.
#' @export
prep_flag_contaminants <- function(data,
                                    method       = "auto",
                                    context      = c("general", "icu_device"),
                                    organism_col = "organism_normalized",
                                    specimen_col = "specimen_type") {
  context <- match.arg(context)

  if (!organism_col %in% names(data)) {
    message(sprintf("[!] Column '%s' not found. Skipping contaminant flagging.", organism_col))
    data$is_contaminant          <- FALSE
    data$contaminant_confidence  <- "insufficient_data"
    data$contaminant_method      <- "skipped"
    return(data)
  }

  if (!specimen_col %in% names(data)) data[[specimen_col]] <- NA_character_

  data$is_contaminant         <- FALSE
  data$contaminant_confidence <- "unknown"
  data$contaminant_method     <- NA_character_

  # Use provided flag if it exists
  if ("pathogen_contaminant" %in% names(data)) {
    data$is_contaminant         <- data$pathogen_contaminant == 1
    data$contaminant_method     <- "provided"
    data$contaminant_confidence <- "high"
    message("Using provided contaminant flags")
    return(data)
  }

  # ICU device context: CoNS with central line in place > 2 days = CLABSI, NOT contaminant
  cons_organisms <- c("staphylococcus epidermidis", "staphylococcus haemolyticus",
                       "staphylococcus hominis", "staphylococcus capitis",
                       "coagulase-negative staphylococci")

  if (context == "icu_device" &&
      all(c("central_line_present", "central_line_days") %in% names(data))) {
    is_clabsi_candidate <- !is.na(data[[organism_col]]) &
      tolower(data[[organism_col]]) %in% cons_organisms &
      !is.na(data[[specimen_col]]) & grepl("blood", tolower(data[[specimen_col]])) &
      !is.na(data$central_line_present) & data$central_line_present == 1 &
      !is.na(data$central_line_days) & data$central_line_days > 2

    # Mark these as CLABSI (not contaminant)
    if (any(is_clabsi_candidate, na.rm = TRUE)) {
      data$is_contaminant[is_clabsi_candidate]         <- FALSE
      data$contaminant_confidence[is_clabsi_candidate] <- "high"
      data$contaminant_method[is_clabsi_candidate]     <- "icu_device_clabsi"
      message(sprintf("[prep_flag_contaminants] %d CoNS row(s) classified as CLABSI (not contaminant).",
                      sum(is_clabsi_candidate, na.rm = TRUE)))
    }
  }

  prev_contaminant <- data$is_contaminant
  prev_method      <- data$contaminant_method

  # Device-based path (general)
  if (method %in% c("auto", "device_based") &&
      all(c("device_inserted", "device_insertion_date", "date_of_culture") %in% names(data))) {
    data <- data %>%
      dplyr::mutate(
        .days_device = as.numeric(difftime(date_of_culture, device_insertion_date, units = "days")),
        is_contaminant = dplyr::case_when(
          tolower(!!rlang::sym(organism_col)) %in% cons_organisms &
            grepl("blood", tolower(!!rlang::sym(specimen_col))) &
            !is.na(.days_device) & .days_device <= 2   ~ TRUE,
          TRUE                                           ~ is_contaminant
        ),
        contaminant_confidence = dplyr::if_else(
          tolower(!!rlang::sym(organism_col)) %in% cons_organisms &
            grepl("blood", tolower(!!rlang::sym(specimen_col))),
          "high", contaminant_confidence
        ),
        contaminant_method = dplyr::if_else(
          is.na(contaminant_method), "device_based", contaminant_method
        ),
        .days_device = NULL
      )

    message(sprintf("Device-based contaminants: %d",
                    sum(data$contaminant_method == "device_based", na.rm = TRUE)))
    prev_contaminant <- data$is_contaminant
    prev_method      <- data$contaminant_method
  }

  # Heuristic path: organism + specimen type
  if (method %in% c("auto", "heuristic")) {
    contam_blood <- prep_get_contaminant_list(syndrome = "Bloodstream infections")
    contam_urine <- prep_get_contaminant_list(syndrome = "Urinary tract infections / pyelonephritis")

    data <- data %>%
      dplyr::mutate(
        .spec_lower = tolower(!!rlang::sym(specimen_col)),
        is_contaminant = dplyr::case_when(
          grepl("blood", .spec_lower) &
            !!rlang::sym(organism_col) %in% contam_blood$names       ~ TRUE,
          grepl("urine", .spec_lower) &
            !!rlang::sym(organism_col) %in% contam_urine$names       ~ TRUE,
          prev_contaminant == TRUE                                     ~ TRUE,
          TRUE                                                          ~ FALSE
        ),
        contaminant_confidence = dplyr::case_when(
          is_contaminant & !is.na(prev_method) & prev_method == "device_based" ~ "high",
          is_contaminant ~ "medium",
          TRUE           ~ "low"
        ),
        contaminant_method = dplyr::if_else(
          is.na(prev_method) | prev_method == "",
          "heuristic", prev_method
        ),
        .spec_lower = NULL
      )
  }

  n_contam  <- sum(data$is_contaminant, na.rm = TRUE)
  n_high    <- sum(data$contaminant_confidence == "high",   na.rm = TRUE)
  n_medium  <- sum(data$contaminant_confidence == "medium", na.rm = TRUE)
  message(sprintf("Contaminants flagged: %d total (high: %d, medium: %d)",
                  n_contam, n_high, n_medium))

  return(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 6c)
# ---------------------------------------------------------------------------

#' Filter Fungal Organisms
#'
#' Removes rows where the organism belongs to a fungal group. Logs the count
#' of rows removed.
#'
#' @param data Data frame.
#' @param group_col Character. Organism group column. Default "organism_group".
#' @param organism_col Character. Normalized organism column used as fallback
#'   when \code{group_col} is absent. Default "organism_normalized".
#'
#' @return Data frame with fungal rows removed.
#' @export
prep_filter_fungal <- function(data,
                                group_col    = "organism_group",
                                organism_col = "organism_normalized") {
  n_before <- nrow(data)

  if (group_col %in% names(data)) {
    is_fungal <- !is.na(data[[group_col]]) &
      grepl("fung|yeast|candida|asperg|cryptococ|mucor|fung",
            tolower(data[[group_col]]), ignore.case = TRUE)
  } else if (organism_col %in% names(data)) {
    is_fungal <- !is.na(data[[organism_col]]) &
      grepl("candida|aspergillus|cryptococcus|mucor|trichosporon|fusarium|rhizopus|pichia",
            tolower(data[[organism_col]]), ignore.case = TRUE)
    message(sprintf("[prep_filter_fungal] '%s' not found; using organism name patterns.", group_col))
  } else {
    warning("[prep_filter_fungal] Neither group nor organism column found. Returning data unchanged.")
    return(data)
  }

  n_fungal <- sum(is_fungal, na.rm = TRUE)
  if (n_fungal > 0) {
    message(sprintf("[prep_filter_fungal] Removing %d fungal row(s) (%.1f%%).",
                    n_fungal, 100 * n_fungal / n_before))
    data <- data[!is_fungal, ]
  } else {
    message("[prep_filter_fungal] No fungal organisms detected.")
  }

  return(data)
}


