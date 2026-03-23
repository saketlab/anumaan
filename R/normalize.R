# normalize.R
# Phase 1: Normalization functions for AMR data preprocessing


# Store package root when normalize.R is sourced (set at source time)
.amrburden_pkg_root <- NULL

# Try to detect package root at source time
local({
  # Get the path of this source file
  source_path <- tryCatch(
    {
      # Check various ways to get the source file path
      if (exists("ofile", envir = parent.frame(2))) {
        get("ofile", envir = parent.frame(2))
      } else {
        # Try to get from sys.frames
        for (i in seq_len(sys.nframe())) {
          env <- sys.frame(i)
          if (exists("ofile", envir = env)) {
            return(get("ofile", envir = env))
          }
        }
        NULL
      }
    },
    error = function(e) NULL
  )

  if (!is.null(source_path) && file.exists(source_path)) {
    # This file is in R/, so package root is parent of R/
    r_dir <- dirname(normalizePath(source_path))
    pkg_root <- dirname(r_dir)

    # Verify it's the anumaan package
    desc_path <- file.path(pkg_root, "DESCRIPTION")
    if (file.exists(desc_path)) {
      .amrburden_pkg_root <<- pkg_root
    }
  }
})


#' Find Package Data File
#'
#' Helper function to locate data files in inst/extdata.
#' Works both when the package is installed and during development.
#'
#' @param filename Character. Name of the file to find (e.g., "organisms.csv").
#' @return Character. Full path to the file, or empty string if not found.
#' @keywords internal
find_extdata_file <- function(filename) {
  # First try system.file (works when package is installed)
  file_path <- system.file("extdata", filename, package = "anumaan")

  if (file_path != "" && file.exists(file_path)) {
    return(file_path)
  }

  # Strategy 1: Use package root detected at source time
  if (!is.null(.amrburden_pkg_root)) {
    extdata_path <- file.path(.amrburden_pkg_root, "inst", "extdata", filename)
    if (file.exists(extdata_path)) {
      return(normalizePath(extdata_path))
    }
  }

  # Strategy 2: Check relative to current working directory
  # Cover common directory structures (including deep nesting like analysis/X/Y/)
  possible_roots <- c(
    ".", # Working dir is package root
    "..", # Working dir is R/ or subdirectory
    "../..", # Working dir is 2 levels deep
    "../../..", # Working dir is 3 levels deep (e.g., analysis/X/Y/)
    "../../../..", # Working dir is 4 levels deep
    "anumaan", # Working dir is parent of package
    "../anumaan", # Working dir is sibling of package
    "../../anumaan", # Working dir is nested sibling
    "../../../anumaan" # Working dir is deeply nested
  )

  for (root in possible_roots) {
    extdata_path <- file.path(root, "inst", "extdata", filename)
    if (file.exists(extdata_path)) {
      return(normalizePath(extdata_path))
    }
  }

  # File not found
  return("")
}


#' Standardize Column Names to Package Convention
#'
#' Maps incoming dataset column names to standardized names used throughout
#' the package. Supports exact matching and optional fuzzy matching for
#' unmatched columns.
#'
#' @param data A data frame with raw column names
#' @param mapping Named list where names are standard column names and values
#'   are character vectors of acceptable aliases. Default uses
#'   \code{default_column_mappings}.
#' @param fuzzy_match Logical. If TRUE, attempts fuzzy matching for unmapped
#'   columns using string distance. Default TRUE.
#' @param fuzzy_threshold Numeric. Maximum string distance (0-1) for fuzzy
#'   matching. Lower values require closer matches. Default 0.3.
#' @param interactive Logical. If TRUE and fuzzy matches found, prompts user
#'   for confirmation. Default FALSE (auto-accept).
#'
#' @return A list with components:
#'   \itemize{
#'     \item data: Data frame with standardized column names
#'     \item mapping_log: List documenting which columns were mapped and how
#'     \item unmapped: Character vector of columns that couldn't be mapped
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' raw_data <- data.frame(
#'   PatientID = 1:10,
#'   Organism = rep("E. coli", 10),
#'   Drug = rep("Ampicillin", 10)
#' )
#' result <- standardize_column_names(raw_data)
#' clean_data <- result$data
#' }
standardize_column_names <- function(data,
                                     mapping = default_column_mappings,
                                     fuzzy_match = TRUE,
                                     fuzzy_threshold = 0.3,
                                     interactive = FALSE) {
  original_names <- names(data)
  new_names <- names(data)
  mapping_log <- list()

  # Phase 1: Exact matching
  for (std_name in names(mapping)) {
    aliases <- mapping[[std_name]]
    matches <- which(original_names %in% aliases)

    if (length(matches) > 1) {
      warning(sprintf(
        "Multiple columns match '%s': %s. Using first match.",
        std_name,
        paste(original_names[matches], collapse = ", ")
      ))
      matches <- matches[1]
    }

    if (length(matches) == 1) {
      new_names[matches] <- std_name
      mapping_log[[std_name]] <- list(
        original = original_names[matches],
        method = "exact_match"
      )
    }
  }

  # Phase 2: Fuzzy matching (optional)
  if (fuzzy_match) {
    unmapped_idx <- which(new_names == original_names)
    unmapped <- original_names[unmapped_idx]

    for (std_name in names(mapping)) {
      if (std_name %in% new_names) next # Already mapped

      # Calculate string distances
      distances <- stringdist::stringdist(
        tolower(std_name),
        tolower(unmapped),
        method = "jw" # Jaro-Winkler distance
      )

      best_match_idx <- which.min(distances)

      if (length(best_match_idx) > 0 && distances[best_match_idx] < fuzzy_threshold) {
        if (interactive) {
          message(sprintf(
            "Fuzzy match: '%s' -> '%s' (distance: %.2f)",
            unmapped[best_match_idx],
            std_name,
            distances[best_match_idx]
          ))
          confirm <- readline(prompt = "Accept this mapping? (y/n): ")
          accept <- tolower(confirm) == "y"
        } else {
          accept <- TRUE
          message(sprintf(
            "Auto-accepted fuzzy match: '%s' -> '%s' (distance: %.2f)",
            unmapped[best_match_idx],
            std_name,
            distances[best_match_idx]
          ))
        }

        if (accept) {
          new_names[unmapped_idx[best_match_idx]] <- std_name
          mapping_log[[std_name]] <- list(
            original = unmapped[best_match_idx],
            method = "fuzzy_match",
            distance = distances[best_match_idx]
          )
        }
      }
    }
  }

  # Apply new names
  names(data) <- new_names

  # Report unmapped columns
  unmapped_final <- setdiff(new_names, names(mapping))
  if (length(unmapped_final) > 0) {
    message(
      "Columns not mapped to standard names: ",
      paste(unmapped_final, collapse = ", ")
    )
  }

  return(list(
    data = data,
    mapping_log = mapping_log,
    unmapped = unmapped_final
  ))
}


#' Normalize Antibiotic Names
#'
#' Standardizes antibiotic/drug names by removing punctuation, spaces, and
#' applying consistent lowercase formatting.
#'
#' @param data Data frame containing antibiotic names
#' @param col Character. Name of the column containing antibiotic names.
#'   Default "antibiotic_name".
#' @param custom_map Optional named character vector of additional mappings.
#'
#' @return Data frame with added column \code{antibiotic_normalized}
#'
#' @export
normalize_antibiotic <- function(data, col = "antibiotic_name", custom_map = NULL) {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  # Apply normalization
  data$antibiotic_normalized <- data[[col]] %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[[:punct:][:space:]]+", "") %>%
    stringr::str_trim()

  # Apply any custom mappings
  if (!is.null(custom_map)) {
    for (variant in names(custom_map)) {
      data$antibiotic_normalized[data$antibiotic_normalized == variant] <-
        custom_map[variant]
    }
  }

  message(sprintf(
    "Normalized %d antibiotic names (%d unique values)",
    nrow(data),
    dplyr::n_distinct(data$antibiotic_normalized)
  ))

  return(data)
}


#' Normalize Organism Names
#'
#' Normalizes organism names using organisms.csv reference file.
#' Automatically handles abbreviations (E. coli), case variations, and typos.
#'
#' @param data Data frame with organism column
#' @param organism_col Character. Organism column name. Default "organism_name".
#' @param add_organism_group Logical. Add organism_group column from CSV. Default TRUE.
#' @param add_resistance_flags Logical. Add resistance flag columns (0/1). Default TRUE.
#'   Creates: is_MRSA, is_MSSA, is_MRCONS, is_MSCONS
#'
#' @return Data frame with organism_normalized, organism_group, and optionally resistance flag columns
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(organism = c("E. coli", "S. aureus", "MRSA"))
#' result <- normalize_organism(data, organism_col = "organism")
#' }
normalize_organism <- function(data,
                               organism_col = "organism_name",
                               add_organism_group = TRUE,
                               add_resistance_flags = TRUE) {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", organism_col))
  }

  # Load organisms.csv reference
  csv_path <- find_extdata_file("organisms.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    warning("organisms.csv not found in inst/extdata/. Organism normalization skipped.")
    data$organism_normalized <- data[[organism_col]]
    if (add_organism_group) data$organism_group <- NA_character_
    return(data)
  }

  # Load reference organisms
  org_ref <- readr::read_csv(csv_path, show_col_types = FALSE)

  # Helper function: extract keywords from organism name
  extract_keywords <- function(name) {
    if (is.na(name) || name == "") {
      return(character(0))
    }

    # Remove special characters, lowercase, split by spaces
    cleaned <- tolower(gsub("[^a-z0-9\\s]", " ", name))
    words <- unlist(strsplit(cleaned, "\\s+"))

    # Remove common filler words
    filler <- c(
      "spp", "sp", "species", "not", "specified", "other", "resistant",
      "positive", "negative", "gen", "generation"
    )
    words <- words[!words %in% filler & nchar(words) > 1]

    return(unique(words))
  }

  # Create lookup: extract keywords from reference organisms
  org_ref$keywords <- lapply(org_ref$organism_name, extract_keywords)
  org_ref$ref_lower <- tolower(trimws(org_ref$organism_name))

  # Process each unique input organism
  data$temp_org_input <- trimws(as.character(data[[organism_col]]))

  # ============================================================================
  # STEP 1: Detect and flag methicillin resistance patterns (MRSA, MRCONS, etc.)
  # ============================================================================
  if (add_resistance_flags) {
    data$is_MRSA <- 0L
    data$is_MSSA <- 0L
    data$is_MRCONS <- 0L
    data$is_MSCONS <- 0L

    # Detect MRSA patterns (Methicillin-resistant Staphylococcus aureus)
    # Matches: "MRSA", "S. aureus (MRSA)", "Staphylococcus aureus MRSA", "methicillin resistant S. aureus"
    mrsa_pattern <- "\\bmrsa\\b|\\(mrsa\\)|methicillin[- ]*resist[a-z]*[- ]*s[a-z]*[- ]*aureus|s\\.?\\s*aureus.*mrsa|staphylococcus\\s+aureus.*mrsa"
    data$is_MRSA[grepl(mrsa_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    # Detect MSSA patterns (Methicillin-sensitive Staphylococcus aureus)
    # Matches: "MSSA", "S. aureus (MSSA)", "methicillin sensitive S. aureus"
    mssa_pattern <- "\\bmssa\\b|\\(mssa\\)|methicillin[- ]*sens[a-z]*[- ]*s[a-z]*[- ]*aureus|methicillin[- ]*suscept[a-z]*[- ]*s[a-z]*[- ]*aureus|s\\.?\\s*aureus.*mssa|staphylococcus\\s+aureus.*mssa"
    data$is_MSSA[grepl(mssa_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    # Detect MRCONS/MR-CONS patterns (Methicillin-resistant Coagulase-negative staphylococci)
    # Matches: "MRCONS", "MR-CONS", "MR CONS", "(MRCONS)", "methicillin resistant coagulase negative"
    mrcons_pattern <- "\\bmr[- ]?cons\\b|\\bmrcos\\b|\\(mr[- ]?cons\\)|methicillin[- ]*resist[a-z]*[- ]*coagulase[- ]*neg"
    data$is_MRCONS[grepl(mrcons_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    # Detect MSCONS patterns (Methicillin-sensitive Coagulase-negative staphylococci)
    # Matches: "MSCONS", "MS-CONS", "methicillin sensitive coagulase negative"
    mscons_pattern <- "\\bms[- ]?cons\\b|\\(ms[- ]?cons\\)|methicillin[- ]*sens[a-z]*[- ]*coagulase[- ]*neg|methicillin[- ]*suscept[a-z]*[- ]*coagulase[- ]*neg"
    data$is_MSCONS[grepl(mscons_pattern, data$temp_org_input, ignore.case = TRUE)] <- 1L

    # Map MRSA/MSSA to Staphylococcus aureus
    data$temp_org_input[data$is_MRSA == 1 | data$is_MSSA == 1] <- "Staphylococcus aureus"

    # Map MRCONS/MSCONS/CONS variants to Coagulase-negative staphylococci
    # Matches: "CONS", "CNS", "coagulase negative", "coag-neg"
    cons_pattern <- "\\bcons\\b|\\bcns\\b|\\(cons\\)|\\(cns\\)|coagulase[- ]*neg|\\bco[a]?g[- ]*neg"
    data$temp_org_input[data$is_MRCONS == 1 | data$is_MSCONS == 1 | grepl(cons_pattern, data$temp_org_input, ignore.case = TRUE)] <- "Coagulase-negative staphylococci"
  } else {
    # Still normalize CONS/MRSA/MSSA variants even without flags
    cons_pattern <- "\\bcons\\b|\\bcns\\b|\\(cons\\)|\\(cns\\)|coagulase[- ]*neg|\\bco[a]?g[- ]*neg|\\bmr[- ]?cons\\b|\\bmrcos\\b|\\bms[- ]?cons\\b"
    data$temp_org_input[grepl(cons_pattern, data$temp_org_input, ignore.case = TRUE)] <- "Coagulase-negative staphylococci"

    mrsa_mssa_pattern <- "\\bmrsa\\b|\\(mrsa\\)|\\bmssa\\b|\\(mssa\\)"
    data$temp_org_input[grepl(mrsa_mssa_pattern, data$temp_org_input, ignore.case = TRUE)] <- "Staphylococcus aureus"
  }

  # ============================================================================
  # STEP 2: Standardize species abbreviations
  # ============================================================================
  # "sp." -> "spp." (plural form is standard)
  data$temp_org_input <- gsub("\\bsp\\.\\s*$", "spp.", data$temp_org_input, ignore.case = TRUE)
  data$temp_org_input <- gsub("\\bsp\\s+", "spp. ", data$temp_org_input, ignore.case = TRUE)
  # Also handle cases like "Genus sp -" -> "Genus spp. -"
  data$temp_org_input <- gsub("\\bsp\\b", "spp.", data$temp_org_input, ignore.case = TRUE)

  # Standardize "non fermenting" variations
  data$temp_org_input <- gsub("non[- ]?ferm[ea]nt[ia]ng", "Non-fermenting", data$temp_org_input, ignore.case = TRUE)

  data$organism_normalized <- NA_character_

  unique_inputs <- unique(data$temp_org_input[!is.na(data$temp_org_input) & data$temp_org_input != ""])

  # Create a mapping for each unique input
  organism_map <- setNames(rep(NA_character_, length(unique_inputs)), unique_inputs)

  for (input_org in unique_inputs) {
    input_keywords <- extract_keywords(input_org)

    if (length(input_keywords) == 0) {
      organism_map[input_org] <- tolower(input_org)
      next
    }

    # Score each reference organism by keyword overlap
    scores <- sapply(1:nrow(org_ref), function(i) {
      ref_keywords <- org_ref$keywords[[i]]

      if (length(ref_keywords) == 0) {
        return(0)
      }

      # Count overlapping keywords
      overlap <- sum(input_keywords %in% ref_keywords)

      # Exact match bonus
      if (tolower(input_org) == org_ref$ref_lower[i]) {
        return(1000)
      }

      # Partial string match bonus
      if (grepl(org_ref$ref_lower[i], tolower(input_org), fixed = TRUE) ||
        grepl(tolower(input_org), org_ref$ref_lower[i], fixed = TRUE)) {
        overlap <- overlap + 2
      }

      # Return score: overlap / total unique keywords (Jaccard similarity)
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) {
        return(0)
      }

      return(overlap / total_unique)
    })

    # Pick best match if score > 0.3 (at least 30% keyword overlap)
    max_score <- max(scores)
    if (max_score > 0.3) {
      # If multiple organisms have the same score, use string distance as tie-breaker
      top_matches <- which(scores == max_score)

      if (length(top_matches) == 1) {
        organism_map[input_org] <- org_ref$ref_lower[top_matches[1]]
      } else {
        # Tie-breaker: use string distance
        input_lower <- tolower(input_org)
        tie_distances <- sapply(top_matches, function(idx) {
          adist(input_lower, org_ref$ref_lower[idx], ignore.case = TRUE)[1, 1]
        })
        best_tie_idx <- top_matches[which.min(tie_distances)]
        organism_map[input_org] <- org_ref$ref_lower[best_tie_idx]
      }
    } else {
      # Fallback: try string distance matching for typos
      input_lower <- tolower(input_org)

      # Extract genus (first meaningful word, skip common prefixes)
      input_words <- strsplit(input_lower, "\\s+")[[1]]
      input_genus <- input_words[1]
      if (input_genus %in% c("non", "multi", "methicillin")) {
        input_genus <- if (length(input_words) > 1) input_words[2] else input_words[1]
      }

      # Calculate genus distance for all reference organisms
      genus_info <- lapply(org_ref$ref_lower, function(ref) {
        ref_words <- strsplit(ref, "\\s+")[[1]]
        ref_genus <- ref_words[1]
        if (ref_genus %in% c("non", "multi", "methicillin")) {
          ref_genus <- if (length(ref_words) > 1) ref_words[2] else ref_words[1]
        }
        genus_dist <- adist(input_genus, ref_genus, ignore.case = TRUE)[1, 1]
        list(ref_genus = ref_genus, genus_dist = genus_dist)
      })

      genus_distances <- sapply(genus_info, function(x) x$genus_dist)

      # Genus matching threshold
      genus_threshold <- if (nchar(input_genus) >= 4) 2 else 1

      # Filter to organisms with matching genus
      genus_matches <- which(genus_distances <= genus_threshold)

      if (length(genus_matches) > 0) {
        # Among genus matches, find the best full string match
        best_match_idx <- genus_matches[1]
        best_match_dist <- Inf

        for (idx in genus_matches) {
          full_dist <- adist(input_lower, org_ref$ref_lower[idx], ignore.case = TRUE)[1, 1]
          if (full_dist < best_match_dist) {
            best_match_dist <- full_dist
            best_match_idx <- idx
          }
        }

        organism_map[input_org] <- org_ref$ref_lower[best_match_idx]
      } else {
        # No genus match - keep original lowercased
        organism_map[input_org] <- tolower(input_org)
      }
    }
  }

  # Apply mapping to data
  data$organism_normalized <- organism_map[data$temp_org_input]
  data$organism_normalized[is.na(data$temp_org_input) | data$temp_org_input == ""] <- NA_character_

  # Add organism_group
  if (add_organism_group) {
    data <- data %>%
      dplyr::left_join(
        org_ref %>% dplyr::select(ref_lower, organism_group),
        by = c("organism_normalized" = "ref_lower")
      )
  }

  # Special case overrides: organisms not in reference CSV
  if (add_organism_group && "organism_group" %in% names(data)) {
    norm_lower <- tolower(data$organism_normalized)

    # Gram-negative bacilli overrides
    data$organism_group[!is.na(norm_lower) & grepl("chryseobacterium", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("vibrio", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("ralstonia", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("delftia", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("elizabethkingia", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("ochrobacterium", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("shewanella", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("non fermenter", norm_lower)] <- "Gram-negative bacilli"
    data$organism_group[!is.na(norm_lower) & grepl("non lactose fermenting", norm_lower)] <- "Gram-negative bacilli"

    # Enterobacterales overrides
    data$organism_group[!is.na(norm_lower) & grepl("pantoea", norm_lower)] <- "Enterobacterales"
    data$organism_group[!is.na(norm_lower) & grepl("kluyvera", norm_lower)] <- "Enterobacterales"
    data$organism_group[!is.na(norm_lower) & grepl("leclercia", norm_lower)] <- "Enterobacterales"
  }

  # Clean up
  data$temp_org_input <- NULL

  # Report
  n_matched <- sum(!is.na(data$organism_normalized))
  message(sprintf(
    "Normalized %d/%d organisms (%.1f%%)",
    n_matched, nrow(data), 100 * n_matched / nrow(data)
  ))

  n_unique <- dplyr::n_distinct(data$organism_normalized, na.rm = TRUE)
  message(sprintf("Result: %d unique organisms", n_unique))

  # Report resistance flags
  if (add_resistance_flags) {
    n_mrsa <- sum(data$is_MRSA == 1, na.rm = TRUE)
    n_mssa <- sum(data$is_MSSA == 1, na.rm = TRUE)
    n_mrcons <- sum(data$is_MRCONS == 1, na.rm = TRUE)
    n_mscons <- sum(data$is_MSCONS == 1, na.rm = TRUE)
    if (n_mrsa > 0 || n_mssa > 0 || n_mrcons > 0 || n_mscons > 0) {
      message(sprintf(
        "Resistance flags: MRSA=%d, MSSA=%d, MRCONS=%d, MSCONS=%d",
        n_mrsa, n_mssa, n_mrcons, n_mscons
      ))
    }
  }

  if (add_organism_group) {
    n_with_group <- sum(!is.na(data$organism_group))
    n_without_group <- sum(is.na(data$organism_group) & !is.na(data$organism_normalized))
    message(sprintf(
      "With organism_group: %d (%.1f%%)",
      n_with_group, 100 * n_with_group / nrow(data)
    ))
    if (n_without_group > 0) {
      message(sprintf("Warning: %d organisms matched but missing organism_group", n_without_group))
    }
  }

  return(data)
}


#' Parse Dates Safely
#'
#' Attempts to parse date columns using multiple common formats via lubridate.
#' Returns Date objects or NA for unparseable values.
#'
#' @param data Data frame containing date columns
#' @param date_columns Character vector of column names to parse as dates.
#'   Default includes common date fields.
#'
#' @return Data frame with date columns converted to Date class
#'
#' @export
parse_dates <- function(data, date_columns = c(
                          "date_of_admission",
                          "date_of_culture",
                          "date_of_final_outcome",
                          "DOB"
                        )) {
  for (col in date_columns) {
    if (!col %in% names(data)) {
      message(sprintf("Date column '%s' not found, skipping", col))
      next
    }

    # Skip if already Date class
    if (inherits(data[[col]], "Date")) {
      next
    }

    original_na <- sum(is.na(data[[col]]))

    # Parse using multiple format orders
    data[[col]] <- suppressWarnings(
      lubridate::parse_date_time(
        data[[col]],
        orders = c(
          "Ymd", "ymd", "dmy", "mdy", "Y-m-d", "d-m-Y", "d/m/Y",
          "Ymd HMS", "ymd HMS", "dmy HMS", "mdy HMS"
        ),
        tz = "UTC"
      )
    ) %>% lubridate::as_date()

    new_na <- sum(is.na(data[[col]]))
    parsed_na <- new_na - original_na

    if (parsed_na > 0) {
      warning(sprintf(
        "Column '%s': %d values could not be parsed as dates",
        col, parsed_na
      ))
    }

    message(sprintf(
      "Parsed '%s': %d dates (%d NA)",
      col, sum(!is.na(data[[col]])), new_na
    ))
  }

  return(data)
}


#' Standardize Gender Values
#'
#' Maps various gender representations to standard "M" or "F" values.
#'
#' @param data Data frame containing gender column
#' @param col Character. Name of gender column. Default "gender".
#'
#' @return Data frame with standardized gender values
#'
#' @export
standardize_gender <- function(data, col = "gender") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  original_values <- unique(data[[col]])

  data[[col]] <- dplyr::case_when(
    toupper(data[[col]]) %in% c("M", "MALE", "MAN") ~ "M",
    toupper(data[[col]]) %in% c("F", "FEMALE", "WOMAN") ~ "F",
    TRUE ~ NA_character_
  )

  na_count <- sum(is.na(data[[col]]))
  if (na_count > 0) {
    warning(sprintf(
      "Gender column: %d values could not be standardized to M/F",
      na_count
    ))
  }

  message(sprintf(
    "Standardized gender: M=%d, F=%d, NA=%d",
    sum(data[[col]] == "M", na.rm = TRUE),
    sum(data[[col]] == "F", na.rm = TRUE),
    na_count
  ))

  return(data)
}


#' Standardize Outcome Values
#'
#' Maps various outcome representations to standard values:
#' "Died", "Discharged", "LAMA", "Unknown"
#'
#' @param data Data frame containing outcome column
#' @param col Character. Name of outcome column. Default "final_outcome".
#'
#' @return Data frame with standardized outcome values
#'
#' @export
standardize_outcome <- function(data, col = "final_outcome") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  data[[col]] <- dplyr::case_when(
    toupper(data[[col]]) %in% c("DIED", "DEATH", "DECEASED", "EXPIRED") ~ "Died",
    toupper(data[[col]]) %in% c("DISCHARGED", "DISCHARGE", "ALIVE", "SURVIVED") ~ "Discharged",
    toupper(data[[col]]) %in% c("LAMA", "DAMA", "LEFT", "ABSCONDED") ~ "LAMA",
    TRUE ~ "Unknown"
  )

  # Summary
  outcome_counts <- table(data[[col]])
  message("Outcome distribution:")
  print(outcome_counts)

  return(data)
}


#' Standardize Susceptibility Values
#'
#' Maps various susceptibility/resistance representations to standard
#' "S" (Susceptible), "R" (Resistant), "I" (Intermediate) values.
#'
#' @param data Data frame containing susceptibility column
#' @param col Character. Name of susceptibility column. Default "antibiotic_value".
#'
#' @return Data frame with standardized susceptibility values and added column
#'   \code{antibiotic_value_std}
#'
#' @export
standardize_susceptibility <- function(data, col = "antibiotic_value") {
  if (!col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", col))
  }

  data$antibiotic_value_std <- dplyr::case_when(
    toupper(data[[col]]) %in% c("R", "RESISTANT", "RES") ~ "R",
    toupper(data[[col]]) %in% c("S", "SENSITIVE", "SUSCEPTIBLE", "SUS") ~ "S",
    toupper(data[[col]]) %in% c("I", "INTERMEDIATE", "INT") ~ "I",
    TRUE ~ NA_character_
  )

  na_count <- sum(is.na(data$antibiotic_value_std))
  if (na_count > 0) {
    warning(sprintf(
      "Susceptibility column: %d values could not be standardized to S/R/I",
      na_count
    ))
  }

  # Summary
  susc_counts <- table(data$antibiotic_value_std, useNA = "ifany")
  message("Susceptibility distribution:")
  print(susc_counts)

  return(data)
}


#' Load RR (Relative Risk) Reference Data
#'
#' Loads the rr_list_gbd.xlsx file containing pathogen-drug combinations
#' with relative risk values for burden estimation.
#'
#' @return Data frame with Pathogen, Drug, and RR columns
#' @export
#'
#' @examples
#' \dontrun{
#' rr_data <- load_rr_reference()
#' head(rr_data)
#' }
load_rr_reference <- function() {
  # Find Excel file
  xlsx_path <- find_extdata_file("rr_list_gbd.xlsx")

  if (xlsx_path == "" || !file.exists(xlsx_path)) {
    stop("rr_list_gbd.xlsx not found in inst/extdata/")
  }

  # Load Excel file
  rr_data <- readxl::read_excel(xlsx_path, sheet = "Sheet2")

  # Clean column names
  names(rr_data) <- c(
    "pathogen",
    "drug",
    "sample_size",
    "mean_rr",
    "lower_bound",
    "upper_bound"
  )

  # Remove rows with NA pathogen (like footer notes)
  rr_data <- rr_data %>%
    dplyr::filter(!is.na(pathogen), !is.na(drug))

  message(sprintf(
    "✓ Loaded %d pathogen-drug RR combinations",
    nrow(rr_data)
  ))

  return(rr_data)
}


#' Add RR Pathogen and Drug Mappings
#'
#' Adds rr_pathogen and rr_drug columns to data by matching against
#' the GBD RR reference list.
#'
#' @param data Data frame with organism and antibiotic columns
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_class".
#'
#' @return Data frame with rr_pathogen and rr_drug columns added
#' @export
#'
#' @examples
#' \dontrun{
#' data <- add_rr_mappings(data)
#' }
add_rr_mappings <- function(data,
                            organism_col = "organism_normalized",
                            antibiotic_col = "antibiotic_class") {
  if (!organism_col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping RR mapping.", organism_col))
    return(data)
  }

  # Load RR reference
  rr_ref <- load_rr_reference()

  # Create lowercase versions for matching
  data$temp_org_lower <- tolower(trimws(data[[organism_col]]))
  rr_ref$temp_path_lower <- tolower(trimws(rr_ref$pathogen))

  # Join to add rr_pathogen
  data <- data %>%
    dplyr::left_join(
      rr_ref %>%
        dplyr::select(temp_path_lower, pathogen) %>%
        dplyr::distinct(),
      by = c("temp_org_lower" = "temp_path_lower")
    ) %>%
    dplyr::rename(rr_pathogen = pathogen)

  # Join to add rr_drug (if antibiotic_col exists)
  if (antibiotic_col %in% names(data)) {
    data$temp_drug_lower <- tolower(trimws(data[[antibiotic_col]]))
    rr_ref$temp_drug_ref_lower <- tolower(trimws(rr_ref$drug))

    data <- data %>%
      dplyr::left_join(
        rr_ref %>%
          dplyr::select(temp_drug_ref_lower, drug) %>%
          dplyr::distinct(),
        by = c("temp_drug_lower" = "temp_drug_ref_lower")
      ) %>%
      dplyr::rename(rr_drug = drug)

    data$temp_drug_lower <- NULL
  }

  # Clean up
  data$temp_org_lower <- NULL

  # Report
  n_pathogen_mapped <- sum(!is.na(data$rr_pathogen))
  message(sprintf(
    "RR pathogen mapping: %d/%d (%.1f%%)",
    n_pathogen_mapped, nrow(data),
    100 * n_pathogen_mapped / nrow(data)
  ))

  if (antibiotic_col %in% names(data)) {
    n_drug_mapped <- sum(!is.na(data$rr_drug))
    message(sprintf(
      "RR drug mapping: %d/%d (%.1f%%)",
      n_drug_mapped, nrow(data),
      100 * n_drug_mapped / nrow(data)
    ))
  }

  return(data)
}


#' Normalize Specimen/Sample Type
#'
#' Normalizes specimen/sample type names and adds sample_category and
#' sterile_classification from the reference CSV file.
#'
#' @param data Data frame with specimen column
#' @param specimen_col Character. Specimen column name. Default "specimen_type".
#' @param add_categories Logical. Add sample_category and sterile_classification.
#'   Default TRUE.
#'
#' @return Data frame with specimen_normalized, sample_category, and
#'   sterile_classification columns
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(specimen = c("Blood", "Urine", "CSF"))
#' result <- normalize_specimen(data, specimen_col = "specimen")
#' }
normalize_specimen <- function(data,
                               specimen_col = "specimen_type",
                               add_categories = TRUE) {
  if (!specimen_col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping specimen normalization.", specimen_col))
    return(data)
  }

  # Load sample type classification CSV
  csv_path <- find_extdata_file("sample_type_classification.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    warning("sample_type_classification.csv not found in inst/extdata/. Specimen normalization skipped.")
    data$specimen_normalized <- data[[specimen_col]]
    if (add_categories) {
      data$sample_category <- NA_character_
      data$sterile_classification <- NA_character_
    }
    return(data)
  }

  # Load reference
  specimen_ref <- readr::read_csv(csv_path, show_col_types = FALSE)

  # Create lowercase working copy and clean messy input
  data$temp_spec_input <- tolower(trimws(data[[specimen_col]]))

  # Strip common suffixes, prefixes, and remarks
  data$temp_spec_clean <- data$temp_spec_input

  # Remove culture-related terms first
  data$temp_spec_clean <- gsub("\\s*culture\\s*/\\s*sensitivity.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*culture.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*c/s.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*sensitivity.*$", "", data$temp_spec_clean, ignore.case = TRUE)

  # Remove content in parentheses or brackets (e.g., "CSF (cerebrospinal fluid)" -> "CSF")
  data$temp_spec_clean <- gsub("\\s*\\(.*?\\)\\s*", " ", data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\s*\\[.*?\\]\\s*", " ", data$temp_spec_clean)

  # Remove descriptive phrases after "from" (e.g., "Urine from Catheter" -> "Urine")
  data$temp_spec_clean <- gsub("\\s*from\\s+.*$", "", data$temp_spec_clean, ignore.case = TRUE)

  # Remove catheter-related descriptors (after the specimen type)
  data$temp_spec_clean <- gsub("\\s*(catheter|foley|clean catch|voided).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*-\\s*(catheter|foley|midstream|peripheral|central).*$", "", data$temp_spec_clean, ignore.case = TRUE)

  # Remove common prefixes before specimen type (e.g., "Midstream Urine" -> "Urine", "Peripheral Blood" -> "Blood")
  data$temp_spec_clean <- gsub("^(midstream|peripheral|central|clean|voided|fresh)\\s+", "", data$temp_spec_clean, ignore.case = TRUE)

  # Remove content after hyphen or colon if it looks like a remark
  data$temp_spec_clean <- gsub("\\s*-\\s*(specimen|sample|report|test|site).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*:\\s*(specimen|sample|report|test|site).*$", "", data$temp_spec_clean, ignore.case = TRUE)

  # Remove common filler words and descriptors
  data$temp_spec_clean <- gsub("\\b(specimen|sample|site|swab from|aspirate from|tip|fluid from)\\b", "", data$temp_spec_clean, ignore.case = TRUE)

  # Clean up extra whitespace
  data$temp_spec_clean <- gsub("\\s+", " ", data$temp_spec_clean)
  data$temp_spec_clean <- trimws(data$temp_spec_clean)

  # Handle common abbreviations
  data$temp_spec_clean <- gsub("\\bet\\b", "endotracheal", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\beta\\b", "endotracheal aspirate", data$temp_spec_clean, ignore.case = TRUE)

  # Handle BAL specifically - remove generic "fluid" when BAL is present
  data$temp_spec_clean <- gsub("\\bbal\\s+fluid\\b", "bal", data$temp_spec_clean, ignore.case = TRUE)

  # For common specimen types, try to extract just the primary keyword
  # This helps with cases like "peripheral blood" -> "blood", "mid-stream urine" -> "urine"
  core_specimen_types <- c(
    "blood", "urine", "csf", "sputum", "stool", "pus", "bile",
    "faeces", "wound", "throat", "nasal", "vaginal", "cervical",
    "bal", "eta", "endotracheal"
  )

  for (spec_type in core_specimen_types) {
    # If the cleaned name contains this specimen type, prioritize it
    if (any(grepl(paste0("\\b", spec_type, "\\b"), data$temp_spec_clean))) {
      # Extract just this specimen type for matching
      idx <- grepl(paste0("\\b", spec_type, "\\b"), data$temp_spec_clean)
      if (any(idx)) {
        # Keep the specimen type keyword, remove other descriptors
        data$temp_spec_clean[idx] <- gsub(paste0(".*\\b(", spec_type, ")\\b.*"), "\\1", data$temp_spec_clean[idx])
      }
    }
  }

  data$temp_spec_clean <- trimws(data$temp_spec_clean)

  # Create reference lookup
  specimen_ref$temp_ref_lower <- tolower(trimws(specimen_ref$specimen_name))

  # Helper function: extract keywords from specimen name
  extract_spec_keywords <- function(name) {
    if (is.na(name) || name == "") {
      return(character(0))
    }

    # Remove special characters, split by spaces/slashes
    cleaned <- tolower(gsub("[^a-z0-9\\s/]", " ", name))
    words <- unlist(strsplit(cleaned, "[\\s/]+"))

    # Remove very short words and common fillers (but keep 2-3 letter abbreviations like "bal", "csf", "eta")
    filler <- c("a", "an", "the", "of", "from", "to", "and", "or", "is", "in", "on")
    words <- words[!words %in% filler & nchar(words) >= 2]

    return(unique(words))
  }

  # Create keyword index for reference specimens
  specimen_ref$keywords <- lapply(specimen_ref$temp_ref_lower, extract_spec_keywords)

  # Initialize mapping
  data$specimen_normalized <- NA_character_
  data$sample_category <- NA_character_
  data$sterile_classification <- NA_character_

  unique_inputs <- unique(data$temp_spec_clean[!is.na(data$temp_spec_clean) & data$temp_spec_clean != ""])
  specimen_map <- setNames(rep(NA_integer_, length(unique_inputs)), unique_inputs)

  # Match each unique input
  for (input_spec in unique_inputs) {
    # First try exact match
    exact_match <- which(specimen_ref$temp_ref_lower == input_spec)
    if (length(exact_match) > 0) {
      specimen_map[input_spec] <- exact_match[1]
      next
    }

    # Try keyword matching
    input_keywords <- extract_spec_keywords(input_spec)

    if (length(input_keywords) == 0) {
      next
    }

    # Score each reference specimen by keyword overlap
    scores <- sapply(1:nrow(specimen_ref), function(i) {
      ref_keywords <- specimen_ref$keywords[[i]]

      if (length(ref_keywords) == 0) {
        return(0)
      }

      # Count overlapping keywords
      overlap <- sum(input_keywords %in% ref_keywords)

      # Partial string match bonus
      if (grepl(specimen_ref$temp_ref_lower[i], input_spec, fixed = TRUE) ||
        grepl(input_spec, specimen_ref$temp_ref_lower[i], fixed = TRUE)) {
        overlap <- overlap + 2
      }

      # Jaccard similarity
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) {
        return(0)
      }

      return(overlap / total_unique)
    })

    # Pick best match if score > 0.3
    max_score <- max(scores)
    if (max_score > 0.3) {
      top_matches <- which(scores == max_score)
      specimen_map[input_spec] <- top_matches[1]
    } else {
      # Try string distance as fallback
      distances <- sapply(specimen_ref$temp_ref_lower, function(ref) {
        adist(input_spec, ref, ignore.case = TRUE)[1, 1]
      })

      min_dist_idx <- which.min(distances)
      min_dist <- distances[min_dist_idx]

      # Accept if distance is small (≤ 2 characters)
      if (min_dist <= 2) {
        specimen_map[input_spec] <- min_dist_idx
      }
    }
  }

  # Apply mapping to data using vectorized join (much faster than row-by-row loop)
  # Create lookup table from specimen_map
  valid_mappings <- specimen_map[!is.na(specimen_map)]
  if (length(valid_mappings) > 0) {
    lookup_df <- data.frame(
      temp_spec_clean = names(valid_mappings),
      ref_idx = as.integer(valid_mappings),
      stringsAsFactors = FALSE
    )
    lookup_df$specimen_normalized <- specimen_ref$specimen_name[lookup_df$ref_idx]
    lookup_df$sample_category <- specimen_ref$sample_category[lookup_df$ref_idx]
    lookup_df$sterile_classification <- specimen_ref$sterile_classification[lookup_df$ref_idx]
    lookup_df$ref_idx <- NULL

    # Remove initialized columns and join
    data$specimen_normalized <- NULL
    data$sample_category <- NULL
    data$sterile_classification <- NULL

    data <- dplyr::left_join(data, lookup_df, by = "temp_spec_clean")
  }

  # Clean up temporary columns
  data$temp_spec_input <- NULL
  data$temp_spec_clean <- NULL

  # Report
  n_matched <- sum(!is.na(data$specimen_normalized))
  message(sprintf(
    "Specimen normalization: %d/%d matched (%.1f%%)",
    n_matched, nrow(data), 100 * n_matched / nrow(data)
  ))

  if (add_categories) {
    n_sterile <- sum(data$sterile_classification == "Sterile site", na.rm = TRUE)
    n_nonsterile <- sum(data$sterile_classification == "Non-sterile site", na.rm = TRUE)
    message(sprintf(
      "Sterile classification: %d sterile, %d non-sterile",
      n_sterile, n_nonsterile
    ))
  }

  return(data)
}


#' Get Contaminant List from Reference File
#'
#' Loads contaminant organisms from common_commensals.csv and filters
#' by syndrome and/or specimen type.
#'
#' @param syndrome Character. Optional syndrome name to filter by
#'   (e.g., "Bloodstream infections", "Urinary tract infections").
#' @param specimen_type Character. Optional specimen type to filter by
#'   (e.g., "Blood culture", "Urine culture").
#' @param return_all Logical. If TRUE, returns all contaminants from all
#'   syndromes. Default FALSE.
#'
#' @return Character vector of contaminant organism names (lowercase)
#'
#' @export
#'
#' @examples
#' # Get all blood culture contaminants
#' get_contaminant_list(syndrome = "Bloodstream infections")
#'
#' # Get all contaminants for a specific specimen type
#' get_contaminant_list(specimen_type = "Blood culture")
#'
#' # Get all contaminants across all syndromes
#' get_contaminant_list(return_all = TRUE)
get_contaminant_list <- function(syndrome = NULL,
                                 specimen_type = NULL,
                                 return_all = FALSE) {
  # Load common commensals reference file
  commensals_path <- find_extdata_file("common_commensals.csv")

  if (commensals_path == "" || !file.exists(commensals_path)) {
    warning("common_commensals.csv not found in inst/extdata/. Returning empty list.")
    return(character(0))
  }

  # Read the CSV file
  commensals_data <- readr::read_csv(
    commensals_path,
    show_col_types = FALSE
  )

  # Filter by syndrome if provided
  if (!is.null(syndrome) && !return_all) {
    commensals_data <- commensals_data %>%
      dplyr::filter(Syndrome == syndrome)
  }

  # Filter by specimen type if provided
  if (!is.null(specimen_type) && !return_all) {
    commensals_data <- commensals_data %>%
      dplyr::filter(`Type of culture/specimen` == specimen_type)
  }

  # Extract and parse the contaminant lists
  contaminants <- commensals_data %>%
    dplyr::pull(`Common commensals`) %>%
    stringr::str_split(";\\s*") %>%
    unlist() %>%
    stringr::str_trim() %>%
    unique()

  # Remove any empty strings
  contaminants <- contaminants[contaminants != ""]

  # Generate flexible patterns for each contaminant
  # This allows matching variations like "Staph" → "Staphylococcus"
  contaminant_patterns <- lapply(contaminants, function(name) {
    name_lower <- tolower(name)

    # Split into genus and species
    parts <- strsplit(name_lower, "\\s+")[[1]]

    patterns <- c(name_lower) # Original full name

    if (length(parts) >= 1) {
      genus <- parts[1]

      # Only add genus as a standalone pattern for genus-level entries
      # (e.g. "Micrococcus species", "Aerococcus spp.") — not for specific species
      # like "Staphylococcus epidermidis", which would incorrectly match S. aureus
      is_genus_level <- length(parts) == 1 ||
        grepl("^sp(p\\.?|ecies)$", parts[2], ignore.case = TRUE)

      # Abbreviation for combining with species: "^s\." matches "s. epidermidis"
      genus_abbrev_base <- paste0("^", substr(genus, 1, 1), "\\.")
      # Standalone abbreviation requires dot + whitespace to avoid matching
      # full genus names: "^a\.\s" matches "A. spp." but NOT "Acinetobacter spp."
      genus_abbrev_standalone <- paste0("^", substr(genus, 1, 1), "\\.\\s")

      if (is_genus_level) {
        patterns <- c(patterns, genus, genus_abbrev_standalone)
      } else {
        patterns <- c(patterns, genus_abbrev_standalone)
      }

      # If there's a species, add genus+species combinations
      if (length(parts) >= 2) {
        species <- parts[2]

        # Only add species as standalone if it is a real binomial species name:
        # - Skip "spp." / "species" — these match any organism with that suffix
        # - Skip when genus contains a hyphen (descriptor format like "coagulase-negative")
        #   because the second word is another genus name (e.g., "Staphylococcus"),
        #   not a true species, and would match unrelated organisms
        is_real_species <- !grepl("^sp(p\\.?|ecies)$", species, ignore.case = TRUE) &&
          !grepl("-", genus)

        if (is_real_species) {
          patterns <- c(patterns, species)
        }
        patterns <- c(
          patterns,
          paste0(genus_abbrev_base, "\\s*", species), # e.g., "^s\.\\s*epidermidis"
          paste(genus, species) # Full name
        )
      }

      # Handle common abbreviations
      if (grepl("staphylococcus", genus)) {
        if (is_genus_level) {
          # Genus-level entry: "staph" alone is safe (e.g., "Staphylococcus spp.")
          patterns <- c(patterns, "staph", "coag.*neg", "coagulase.*negative")
        } else {
          # Species-specific entry: bind "staph" to the species to avoid matching S. aureus
          patterns <- c(patterns, paste0("staph.*", parts[2]), "coag.*neg", "coagulase.*negative")
        }
      }
      if (grepl("streptococcus", genus)) {
        patterns <- c(patterns, "strep")
      }
      if (grepl("escherichia", genus)) {
        patterns <- c(patterns, "e\\.?\\s*coli")
      }
      if (grepl("pseudomonas", genus)) {
        patterns <- c(patterns, "pseudo")
      }
      if (grepl("klebsiella", genus)) {
        patterns <- c(patterns, "kleb")
      }
      if (grepl("acinetobacter", genus)) {
        patterns <- c(patterns, "acin")
      }
      if (grepl("enterococcus", genus)) {
        patterns <- c(patterns, "entero")
      }
      if (grepl("corynebacterium", genus)) {
        patterns <- c(patterns, "coryno", "diphtheroids?")
      }
      if (grepl("bacillus", genus)) {
        patterns <- c(patterns, "bacil")
      }
      if (grepl("micrococcus", genus)) {
        patterns <- c(patterns, "micro")
      }
      if (grepl("cutibacterium|propionibacterium", genus)) {
        patterns <- c(patterns, "propioni", "cuti", "p\\.?\\s*acnes")
      }
      if (grepl("lactobacillus", genus)) {
        patterns <- c(patterns, "lacto")
      }
    }

    # Return unique patterns
    list(
      original = name,
      patterns = unique(patterns)
    )
  })

  # Return structure with both original names and matching patterns
  result <- list(
    names = contaminants,
    patterns = contaminant_patterns
  )

  return(result)
}


#' Check if Organism is a Contaminant
#'
#' Checks if a given organism name matches any known contaminant for a
#' specific syndrome or specimen type using flexible pattern matching.
#'
#' @param organism_name Character vector. Organism name(s) to check.
#' @param syndrome Character. Optional syndrome name to filter contaminants.
#' @param specimen_type Character. Optional specimen type to filter contaminants.
#'
#' @return Logical vector indicating if each organism is a contaminant
#'
#' @export
#'
#' @examples
#' # Check if organism is a blood culture contaminant
#' is_contaminant("Staph epidermidis", syndrome = "Bloodstream infections")
#' is_contaminant("E. coli", syndrome = "Bloodstream infections")
#'
#' # Check multiple organisms
#' is_contaminant(c("Staph", "Klebsiella"), specimen_type = "Blood culture")
is_contaminant <- function(organism_name,
                           syndrome = NULL,
                           specimen_type = NULL) {
  # Get contaminant list with patterns
  contaminant_data <- get_contaminant_list(
    syndrome = syndrome,
    specimen_type = specimen_type
  )

  # Check each organism name against all contaminant patterns
  sapply(organism_name, function(org) {
    if (is.na(org) || org == "") {
      return(FALSE)
    }

    org_lower <- tolower(trimws(org))

    # Check against each contaminant's patterns
    for (contam_info in contaminant_data$patterns) {
      for (pattern in contam_info$patterns) {
        if (grepl(pattern, org_lower, perl = TRUE)) {
          return(TRUE)
        }
      }
    }

    return(FALSE)
  }, USE.NAMES = FALSE)
}


#' Clean Antibiotic Susceptibility Values
#'
#' Extracts clean S/I/R values from messy antibiotic result columns.
#' Handles common patterns like "S (HIGH LEVEL)", "R   Escherichia coli", etc.
#'
#' @param data Data frame with antibiotic susceptibility data
#' @param value_col Character. Column name containing susceptibility values.
#'   Default "antibiotic_value".
#' @param strict Logical. If TRUE, only accept S/I/R values. If FALSE, attempt
#'   to parse from messy strings. Default FALSE.
#'
#' @return Data frame with cleaned antibiotic_value column
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   antibiotic_value = c("S", "R   E. coli", "S (HIGH LEVEL)", "I")
#' )
#' clean_antibiotic_values(data)
#' }
clean_antibiotic_values <- function(data,
                                    value_col = "antibiotic_value",
                                    strict = FALSE) {
  if (!value_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", value_col))
  }

  n_before <- nrow(data)
  n_unique_before <- dplyr::n_distinct(data[[value_col]], na.rm = TRUE)

  message(sprintf(
    "Cleaning antibiotic values: %d unique values found",
    n_unique_before
  ))

  # Create a cleaned column
  data$temp_value_clean <- trimws(as.character(data[[value_col]]))

  if (strict) {
    # Strict mode: only accept exact S/I/R
    data[[value_col]] <- ifelse(
      data$temp_value_clean %in% c("S", "I", "R"),
      data$temp_value_clean,
      NA_character_
    )
  } else {
    # Lenient mode: extract S/I/R from messy strings
    data[[value_col]] <- sapply(data$temp_value_clean, function(val) {
      if (is.na(val) || val == "") {
        return(NA_character_)
      }

      # Convert to uppercase for matching
      val_upper <- toupper(val)

      # Extract first character that is S, I, or R
      # This handles cases like:
      # - "R   Escherichia coli" → "R"
      # - "S (HIGH LEVEL)" → "S"
      # - "I (HIGH LEVEL SYNERGY)" → "I"

      # Check if first character is S/I/R
      first_char <- substr(val_upper, 1, 1)
      if (first_char %in% c("S", "I", "R")) {
        return(first_char)
      }

      # Check if string contains S, I, or R anywhere
      if (grepl("^R\\s", val_upper) || grepl("^R\\(", val_upper) || grepl("^R$", val_upper)) {
        return("R")
      }
      if (grepl("^S\\s", val_upper) || grepl("^S\\(", val_upper) || grepl("^S$", val_upper)) {
        return("S")
      }
      if (grepl("^I\\s", val_upper) || grepl("^I\\(", val_upper) || grepl("^I$", val_upper)) {
        return("I")
      }

      # Check for lowercase s, i, r
      if (val == "s") {
        return("S")
      }
      if (val == "i") {
        return("I")
      }
      if (val == "r") {
        return("R")
      }

      # Special cases
      if (grepl("resistant", val_upper)) {
        return("R")
      }
      if (grepl("susceptible", val_upper) || grepl("sensitive", val_upper)) {
        return("S")
      }
      if (grepl("intermediate", val_upper)) {
        return("I")
      }

      # Cannot parse - return NA
      return(NA_character_)
    }, USE.NAMES = FALSE)
  }

  # Clean up temp column
  data$temp_value_clean <- NULL

  n_unique_after <- dplyr::n_distinct(data[[value_col]], na.rm = TRUE)
  n_na <- sum(is.na(data[[value_col]]))

  message(sprintf(
    "Cleaned: %d unique values → %d (S/I/R)",
    n_unique_before, n_unique_after
  ))

  if (n_na > 0) {
    message(sprintf(
      "⚠ %d values could not be parsed (%.1f%%)",
      n_na, 100 * n_na / n_before
    ))
  }

  # Show distribution
  value_dist <- table(data[[value_col]], useNA = "ifany")
  message("\nValue distribution:")
  print(value_dist)

  return(data)
}


#' Recode Intermediate (I) Susceptibility Values
#'
#' Converts "I" (Intermediate) values to either "S" or "R" based on antibiotic type.
#' For Colistin: I -> S (following clinical guidelines)
#' For all other antibiotics: I -> R (conservative approach for surveillance)
#'
#' @param data Data frame with antibiotic susceptibility data
#' @param antibiotic_col Character. Column name with antibiotic names.
#'   Default "antibiotic_name".
#' @param value_col Character. Column name with S/I/R values.
#'   Default "antibiotic_value".
#' @param colistin_to_s Logical. Convert Colistin I to S. Default TRUE.
#' @param others_to_r Logical. Convert other antibiotics' I to R. Default TRUE.
#'
#' @return Data frame with recoded intermediate values
#' @export
#'
#' @examples
#' \dontrun{
#' # Recode all I values according to standard rules
#' data <- recode_intermediate(data)
#'
#' # Only recode Colistin
#' data <- recode_intermediate(data, others_to_r = FALSE)
#' }
recode_intermediate <- function(data,
                                antibiotic_col = "antibiotic_name",
                                value_col = "antibiotic_value",
                                colistin_to_s = TRUE,
                                others_to_r = TRUE) {
  # Input validation
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", antibiotic_col))
  }
  if (!value_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", value_col))
  }

  # Count I values before recoding
  n_intermediate_before <- sum(data[[value_col]] == "I", na.rm = TRUE)

  if (n_intermediate_before == 0) {
    message("No intermediate (I) values found. Returning data unchanged.")
    return(data)
  }

  # Create working copy
  data$temp_value <- data[[value_col]]
  data$temp_antibiotic <- tolower(trimws(data[[antibiotic_col]]))

  # Recode Colistin I -> S
  if (colistin_to_s) {
    colistin_i_idx <- which(
      data$temp_value == "I" &
        grepl("colistin", data$temp_antibiotic, ignore.case = TRUE)
    )

    if (length(colistin_i_idx) > 0) {
      data$temp_value[colistin_i_idx] <- "S"
      message(sprintf("Recoded %d Colistin I -> S", length(colistin_i_idx)))
    }
  }

  # Recode all other antibiotics I -> R
  if (others_to_r) {
    other_i_idx <- which(
      data$temp_value == "I" &
        !grepl("colistin", data$temp_antibiotic, ignore.case = TRUE)
    )

    if (length(other_i_idx) > 0) {
      data$temp_value[other_i_idx] <- "R"
      message(sprintf("Recoded %d non-Colistin I -> R", length(other_i_idx)))
    }
  }

  # Update original column
  data[[value_col]] <- data$temp_value

  # Clean up
  data$temp_value <- NULL
  data$temp_antibiotic <- NULL

  # Summary
  n_intermediate_after <- sum(data[[value_col]] == "I", na.rm = TRUE)
  message(sprintf(
    "Intermediate values: %d -> %d",
    n_intermediate_before, n_intermediate_after
  ))

  return(data)
}


#' Normalize Antibiotic Names
#'
#' Standardizes antibiotic names using fuzzy matching against WHO reference.
#' Similar approach to organism normalization.
#'
#' @param data Data frame with antibiotic data
#' @param antibiotic_col Character. Column name with antibiotic names.
#'   Default "antibiotic_name".
#' @param who_table Data frame. WHO AWaRe classification table. If NULL,
#'   loads from inst/extdata/WHO_aware_class.csv.
#' @param add_class Logical. Add antibiotic_class column. Default TRUE.
#' @param add_aware Logical. Add aware_category column. Default TRUE.
#'
#' @return Data frame with antibiotic_normalized, antibiotic_class, aware_category
#' @export
#'
#' @examples
#' \dontrun{
#' data <- normalize_antibiotic(data, antibiotic_col = "antibiotic_name")
#' }
normalize_antibiotic <- function(data,
                                 antibiotic_col = "antibiotic_name",
                                 who_table = NULL,
                                 add_class = TRUE,
                                 add_aware = TRUE) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", antibiotic_col))
  }

  # Load WHO reference if not provided
  if (is.null(who_table)) {
    who_path <- find_extdata_file("WHO_aware_class.csv")

    if (who_path == "" || !file.exists(who_path)) {
      warning("WHO_aware_class.csv not found in inst/extdata/. Antibiotic normalization skipped.")
      data$antibiotic_normalized <- data[[antibiotic_col]]
      if (add_class) data$antibiotic_class <- NA_character_
      if (add_aware) data$aware_category <- NA_character_
      return(data)
    }

    who_table <- readr::read_csv(who_path, show_col_types = FALSE)
  }

  # Standardize WHO table column names
  if ("Antibiotic" %in% names(who_table) && !"antibiotic_name" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(antibiotic_name = Antibiotic)
  }
  if ("Class" %in% names(who_table) && !"antibiotic_class" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(antibiotic_class = Class)
  }
  if ("Category " %in% names(who_table) && !"aware_category" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(aware_category = `Category `)
  } else if ("Category" %in% names(who_table) && !"aware_category" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(aware_category = Category)
  }

  # Ensure WHO table has required columns
  if (!all(c("antibiotic_name", "antibiotic_class") %in% names(who_table))) {
    stop("WHO table must have 'antibiotic_name' and 'antibiotic_class' columns")
  }

  n_unique_before <- dplyr::n_distinct(data[[antibiotic_col]], na.rm = TRUE)

  message(sprintf(
    "Normalizing %d unique antibiotic names against WHO reference (%d antibiotics)...",
    n_unique_before, nrow(who_table)
  ))

  # Helper function: extract keywords from antibiotic name
  extract_abx_keywords <- function(name) {
    if (is.na(name) || name == "") {
      return(character(0))
    }

    # Remove special characters, lowercase, split by spaces/slashes/hyphens
    cleaned <- tolower(gsub("[^a-z0-9\\s/-]", " ", name))
    words <- unlist(strsplit(cleaned, "[\\s/-]+"))

    # Remove common filler words
    filler <- c("iv", "oral", "injection", "tablet", "mg", "strip", "mic", "done", "by")
    words <- words[!words %in% filler & nchar(words) > 2]

    return(unique(words))
  }

  # Create lookup: extract keywords from WHO antibiotics
  who_table$keywords <- lapply(who_table$antibiotic_name, extract_abx_keywords)
  who_table$ref_lower <- tolower(trimws(who_table$antibiotic_name))

  # Process each unique input antibiotic
  data$temp_abx_input <- trimws(as.character(data[[antibiotic_col]]))
  data$antibiotic_normalized <- NA_character_

  unique_inputs <- unique(data$temp_abx_input[!is.na(data$temp_abx_input) & data$temp_abx_input != ""])

  # Create a mapping for each unique input
  antibiotic_map <- setNames(rep(NA_character_, length(unique_inputs)), unique_inputs)

  for (input_abx in unique_inputs) {
    input_keywords <- extract_abx_keywords(input_abx)

    if (length(input_keywords) == 0) {
      antibiotic_map[input_abx] <- tolower(input_abx)
      next
    }

    # Score each WHO antibiotic by keyword overlap
    scores <- sapply(1:nrow(who_table), function(i) {
      ref_keywords <- who_table$keywords[[i]]

      if (length(ref_keywords) == 0) {
        return(0)
      }

      # Count overlapping keywords
      overlap <- sum(input_keywords %in% ref_keywords)

      # Exact match bonus
      if (tolower(input_abx) == who_table$ref_lower[i]) {
        return(1000)
      }

      # Partial string match bonus
      if (grepl(who_table$ref_lower[i], tolower(input_abx), fixed = TRUE) ||
        grepl(tolower(input_abx), who_table$ref_lower[i], fixed = TRUE)) {
        overlap <- overlap + 2
      }

      # Jaccard similarity
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) {
        return(0)
      }

      return(overlap / total_unique)
    })

    # Pick best match if score > 0.3
    max_score <- max(scores)
    if (max_score > 0.3) {
      # Tie-breaker: use string distance
      top_matches <- which(scores == max_score)

      if (length(top_matches) == 1) {
        antibiotic_map[input_abx] <- who_table$ref_lower[top_matches[1]]
      } else {
        input_lower <- tolower(input_abx)
        tie_distances <- sapply(top_matches, function(idx) {
          adist(input_lower, who_table$ref_lower[idx], ignore.case = TRUE)[1, 1]
        })
        best_tie_idx <- top_matches[which.min(tie_distances)]
        antibiotic_map[input_abx] <- who_table$ref_lower[best_tie_idx]
      }
    } else {
      # Fallback: try string distance matching
      input_lower <- tolower(input_abx)

      distances <- sapply(who_table$ref_lower, function(ref) {
        adist(input_lower, ref, ignore.case = TRUE)[1, 1]
      })

      min_dist_idx <- which.min(distances)
      min_dist <- distances[min_dist_idx]

      # Accept if distance is small (≤ 3 characters)
      if (min_dist <= 3) {
        antibiotic_map[input_abx] <- who_table$ref_lower[min_dist_idx]
      } else {
        # No good match - keep original lowercased
        antibiotic_map[input_abx] <- tolower(input_abx)
      }
    }
  }

  # Apply mapping to data
  data$antibiotic_normalized <- antibiotic_map[data$temp_abx_input]
  data$antibiotic_normalized[is.na(data$temp_abx_input) | data$temp_abx_input == ""] <- NA_character_
  # Guard: empty-string normalized results must also be NA (prevents spurious join matches)
  data$antibiotic_normalized[!is.na(data$antibiotic_normalized) &
    trimws(data$antibiotic_normalized) == ""] <- NA_character_

  # Add antibiotic_class and aware_category
  if (add_class || add_aware) {
    join_cols <- c("antibiotic_normalized" = "ref_lower")
    select_cols <- "ref_lower"

    if (add_class) select_cols <- c(select_cols, "antibiotic_class")
    if (add_aware && "aware_category" %in% names(who_table)) {
      select_cols <- c(select_cols, "aware_category")
    }

    data <- data %>%
      dplyr::left_join(
        who_table %>% dplyr::select(dplyr::all_of(select_cols)),
        by = join_cols
      )
  }

  # Bug fix: empty/NA original antibiotic names must not receive a class or category
  orig_empty <- is.na(data$temp_abx_input) | data$temp_abx_input == ""
  if (add_class && "antibiotic_class" %in% names(data)) {
    data$antibiotic_class[orig_empty] <- NA_character_
  }
  if (add_aware && "aware_category" %in% names(data)) {
    data$aware_category[orig_empty] <- NA_character_
  }

  # Special case overrides: certain antibiotics retain their name as the class
  if (add_class && "antibiotic_class" %in% names(data)) {
    orig_lower <- tolower(data$temp_abx_input)

    # Colistin: no standard antibiotic class - keep as "Colistin" (not Polymyxins)
    data$antibiotic_class[!is.na(orig_lower) & orig_lower == "colistin"] <- "Colistin"


    # Mupirocin High level: not in WHO table - use drug name as class
    mupirocin_hl <- !is.na(orig_lower) &
      grepl("mupirocin", orig_lower) &
      grepl("high", orig_lower) &
      grepl("level", orig_lower)
    data$antibiotic_class[mupirocin_hl] <- "Mupirocin High level"

    # Nalidixic acid: not in WHO table - use drug name as class
    data$antibiotic_class[!is.na(orig_lower) & grepl("nalidixic", orig_lower)] <- "Nalidixic acid"

    # Polymyxin B: WHO table only has route-specific entries (_IV, _oral) - use class name
    data$antibiotic_class[!is.na(orig_lower) & grepl("polymyxin.*b", orig_lower)] <- "Polymyxins"
  }

  # Clean up
  data$temp_abx_input <- NULL

  # Report
  n_unique_after <- dplyr::n_distinct(data$antibiotic_normalized, na.rm = TRUE)
  n_matched <- sum(!is.na(data$antibiotic_normalized))

  message(sprintf(
    "Normalized: %d unique names → %d",
    n_unique_before, n_unique_after
  ))

  if (add_class) {
    n_with_class <- sum(!is.na(data$antibiotic_class))
    message(sprintf(
      "With antibiotic_class: %d (%.1f%%)",
      n_with_class, 100 * n_with_class / nrow(data)
    ))
  }

  if (add_aware && "aware_category" %in% names(data)) {
    n_with_aware <- sum(!is.na(data$aware_category))
    message(sprintf(
      "With AWaRe category: %d (%.1f%%)",
      n_with_aware, 100 * n_with_aware / nrow(data)
    ))
  }

  return(data)
}
