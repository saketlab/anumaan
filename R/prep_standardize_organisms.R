# prep_standardize_organisms.R
# Layer 4a: Organism name standardization
#
# Functions moved here from prep_clean_and_standardize.R:
#   - prep_standardize_organisms
#   - prep_extract_genus
#   - prep_extract_species
#   - prep_assign_organism_group
#
# Reference data: inst/extdata/organisms.csv
#
# New functions: prep_normalize_sp_variants, prep_flag_organism_unmatched


#' Normalize Organism Names
#'
#' Standardizes organism names using organisms.csv reference file.
#' Automatically handles abbreviations (E. coli), case variations, and typos.
#'
#' @param data Data frame with organism column.
#' @param organism_col Character. Organism column name. Default "organism_name".
#' @param add_organism_group Logical. Add organism_group column. Default TRUE.
#' @param add_resistance_flags Logical. Add MRSA/MRCONS flags. Default TRUE.
#'
#' @return Data frame with organism_normalized, organism_group, and optionally
#'   resistance flag columns.
#' @export
prep_standardize_organisms <- function(data,
                                       organism_col         = "organism_name",
                                       add_organism_group   = TRUE,
                                       add_resistance_flags = TRUE) {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", organism_col))
  }

  csv_path <- find_extdata_file("organisms.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    warning("organisms.csv not found in inst/extdata/. Organism normalization skipped.")
    data$organism_normalized <- data[[organism_col]]
    if (add_organism_group) data$organism_group <- NA_character_
    return(data)
  }

  org_ref <- readr::read_csv(csv_path, show_col_types = FALSE)

  extract_keywords <- function(name) {
    if (is.na(name) || name == "") return(character(0))
    cleaned <- tolower(gsub("[^a-z0-9\\s]", " ", name))
    words   <- unlist(strsplit(cleaned, "\\s+"))
    filler  <- c("spp", "sp", "species", "not", "specified", "other", "resistant",
                 "positive", "negative", "gen", "generation")
    words   <- words[!words %in% filler & nchar(words) > 1]
    unique(words)
  }

  org_ref$keywords  <- lapply(org_ref$organism_name, extract_keywords)
  org_ref$ref_lower <- tolower(trimws(org_ref$organism_name))

  data$temp_org_input <- trimws(as.character(data[[organism_col]]))

  # Resistance flags
  if (add_resistance_flags) {
    data$is_MRSA   <- 0L; data$is_MSSA  <- 0L
    data$is_MRCONS <- 0L; data$is_MSCONS <- 0L

    mrsa_pat   <- "\\bmrsa\\b|\\(mrsa\\)|methicillin[- ]*resist[a-z]*[- ]*s[a-z]*[- ]*aureus"
    mssa_pat   <- "\\bmssa\\b|\\(mssa\\)|methicillin[- ]*sens[a-z]*[- ]*s[a-z]*[- ]*aureus|methicillin[- ]*suscept[a-z]*[- ]*s[a-z]*[- ]*aureus"
    mrcons_pat <- "\\bmr[- ]?cons\\b|\\bmrcos\\b|\\(mr[- ]?cons\\)|methicillin[- ]*resist[a-z]*[- ]*coagulase[- ]*neg"
    mscons_pat <- "\\bms[- ]?cons\\b|\\(ms[- ]?cons\\)|methicillin[- ]*sens[a-z]*[- ]*coagulase[- ]*neg"
    cons_pat   <- "\\bcons\\b|\\bcns\\b|\\(cons\\)|\\(cns\\)|coagulase[- ]*neg|\\bco[a]?g[- ]*neg"

    data$is_MRSA[grepl(mrsa_pat, data$temp_org_input, ignore.case = TRUE)]     <- 1L
    data$is_MSSA[grepl(mssa_pat, data$temp_org_input, ignore.case = TRUE)]     <- 1L
    data$is_MRCONS[grepl(mrcons_pat, data$temp_org_input, ignore.case = TRUE)] <- 1L
    data$is_MSCONS[grepl(mscons_pat, data$temp_org_input, ignore.case = TRUE)] <- 1L

    data$temp_org_input[data$is_MRSA == 1 | data$is_MSSA == 1] <- "Staphylococcus aureus"
    data$temp_org_input[data$is_MRCONS == 1 | data$is_MSCONS == 1 |
                          grepl(cons_pat, data$temp_org_input, ignore.case = TRUE)] <-
      "Coagulase-negative staphylococci"
  } else {
    cons_pat <- "\\bcons\\b|\\bcns\\b|\\(cons\\)|\\(cns\\)|coagulase[- ]*neg|\\bco[a]?g[- ]*neg|\\bmr[- ]?cons\\b|\\bmrcos\\b|\\bms[- ]?cons\\b"
    data$temp_org_input[grepl(cons_pat, data$temp_org_input, ignore.case = TRUE)] <-
      "Coagulase-negative staphylococci"
    mrsa_mssa_pat <- "\\bmrsa\\b|\\(mrsa\\)|\\bmssa\\b|\\(mssa\\)"
    data$temp_org_input[grepl(mrsa_mssa_pat, data$temp_org_input, ignore.case = TRUE)] <-
      "Staphylococcus aureus"
  }

  # Normalize sp/spp
  data$temp_org_input <- prep_normalize_sp_variants(data$temp_org_input)
  # Collapse culture-negative statements to missing organism.
  no_growth_pat <- "^(no\\s*growth(\\s*in\\s*culture)?\\.?|culture\\s*negative\\.?|sterile\\s*culture\\.?|nil\\s*growth\\.?)$"
  data$temp_org_input[grepl(no_growth_pat, tolower(trimws(data$temp_org_input)), perl = TRUE)] <- NA_character_

  data$organism_normalized <- NA_character_

  unique_inputs <- unique(data$temp_org_input[!is.na(data$temp_org_input) & data$temp_org_input != ""])
  organism_map  <- setNames(rep(NA_character_, length(unique_inputs)), unique_inputs)

  for (input_org in unique_inputs) {
    input_keywords <- extract_keywords(input_org)

    if (length(input_keywords) == 0) {
      organism_map[input_org] <- tolower(input_org)
      next
    }

    scores <- sapply(seq_len(nrow(org_ref)), function(i) {
      ref_keywords <- org_ref$keywords[[i]]
      if (length(ref_keywords) == 0) return(0)
      overlap <- sum(input_keywords %in% ref_keywords)
      if (tolower(input_org) == org_ref$ref_lower[i]) return(1000)
      if (grepl(org_ref$ref_lower[i], tolower(input_org), fixed = TRUE) ||
            grepl(tolower(input_org), org_ref$ref_lower[i], fixed = TRUE))
        overlap <- overlap + 2
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) return(0)
      overlap / total_unique
    })

    max_score <- max(scores)
    if (max_score > 0.3) {
      top_matches <- which(scores == max_score)
      if (length(top_matches) == 1) {
        organism_map[input_org] <- org_ref$ref_lower[top_matches[1]]
      } else {
        input_lower   <- tolower(input_org)
        tie_distances <- sapply(top_matches, function(idx) {
          adist(input_lower, org_ref$ref_lower[idx], ignore.case = TRUE)[1, 1]
        })
        organism_map[input_org] <- org_ref$ref_lower[top_matches[which.min(tie_distances)]]
      }
    } else {
      input_lower  <- tolower(input_org)
      input_words  <- strsplit(input_lower, "\\s+")[[1]]
      input_genus  <- input_words[1]
      if (input_genus %in% c("non", "multi", "methicillin"))
        input_genus <- if (length(input_words) > 1) input_words[2] else input_words[1]

      genus_info      <- lapply(org_ref$ref_lower, function(ref) {
        ref_words  <- strsplit(ref, "\\s+")[[1]]
        ref_genus  <- ref_words[1]
        if (ref_genus %in% c("non", "multi", "methicillin"))
          ref_genus <- if (length(ref_words) > 1) ref_words[2] else ref_words[1]
        list(ref_genus = ref_genus,
             genus_dist = adist(input_genus, ref_genus, ignore.case = TRUE)[1, 1])
      })
      genus_distances <- sapply(genus_info, function(x) x$genus_dist)
      genus_threshold <- if (nchar(input_genus) >= 4) 2 else 1
      genus_matches   <- which(genus_distances <= genus_threshold)

      if (length(genus_matches) > 0) {
        best_idx  <- genus_matches[1]
        best_dist <- Inf
        for (idx in genus_matches) {
          full_dist <- adist(input_lower, org_ref$ref_lower[idx], ignore.case = TRUE)[1, 1]
          if (full_dist < best_dist) { best_dist <- full_dist; best_idx <- idx }
        }
        organism_map[input_org] <- org_ref$ref_lower[best_idx]
      } else {
        organism_map[input_org] <- tolower(input_org)
      }
    }
  }

  data$organism_normalized <- organism_map[data$temp_org_input]
  data$organism_normalized[is.na(data$temp_org_input) | data$temp_org_input == ""] <- NA_character_

  if (add_organism_group) {
    data <- data %>%
      dplyr::left_join(
        org_ref %>% dplyr::select(ref_lower, organism_group),
        by = c("organism_normalized" = "ref_lower")
      )

    norm_lower <- tolower(data$organism_normalized)
    manual_groups <- list(
      "Gram-negative bacilli" = c("chryseobacterium", "vibrio", "ralstonia", "delftia",
                                  "elizabethkingia", "ochrobacterium", "shewanella",
                                  "non fermenter", "non lactose fermenting"),
      "Enterobacterales"      = c("pantoea", "kluyvera", "leclercia")
    )
    for (grp in names(manual_groups)) {
      for (pat in manual_groups[[grp]])
        data$organism_group[!is.na(norm_lower) & grepl(pat, norm_lower)] <- grp
    }
  }

  data$temp_org_input <- NULL

  n_matched <- sum(!is.na(data$organism_normalized))
  message(sprintf("Normalized %d/%d organisms (%.1f%%)",
                  n_matched, nrow(data), 100 * n_matched / nrow(data)))
  message(sprintf("Result: %d unique organisms",
                  dplyr::n_distinct(data$organism_normalized, na.rm = TRUE)))

  if (add_resistance_flags) {
    for (flag in c("is_MRSA", "is_MSSA", "is_MRCONS", "is_MSCONS")) {
      if (flag %in% names(data))
        message(sprintf("%s: %d", flag, sum(data[[flag]] == 1, na.rm = TRUE)))
    }
  }

  return(data)
}


#' Extract Genus from Organism Name
#'
#' @param data Data frame.
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @return Data frame with org_genus column added.
#' @export
prep_extract_genus <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) stop(sprintf("Column '%s' not found", organism_col))
  data$org_genus <- stringr::str_extract(data[[organism_col]], "^[a-z]+")
  message(sprintf("Extracted genus: %d records, %d unique genera",
                  sum(!is.na(data$org_genus)),
                  length(unique(data$org_genus[!is.na(data$org_genus)]))))
  return(data)
}


#' Extract Species from Organism Name
#'
#' @param data Data frame.
#' @param organism_col Character. Organism column. Default "organism_normalized".
#' @return Data frame with org_species column added.
#' @export
prep_extract_species <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) stop(sprintf("Column '%s' not found", organism_col))
  data$org_species <- stringr::str_extract(data[[organism_col]], "(?<=\\s)[a-z]+")
  data$org_species <- ifelse(data$org_species == "spp", "species", data$org_species)
  message(sprintf("Extracted species: %d records, %d unique",
                  sum(!is.na(data$org_species)),
                  length(unique(data$org_species[!is.na(data$org_species)]))))
  return(data)
}


#' Assign Organism Group
#'
#' @param data Data frame.
#' @param organism_col Character. Normalized organism column. Default "organism_normalized".
#' @return Data frame with org_group column added.
#' @export
prep_assign_organism_group <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) stop(sprintf("Column '%s' not found", organism_col))
  taxonomy <- get_organism_taxonomy()
  data <- data %>%
    dplyr::left_join(taxonomy,
                     by = stats::setNames("organism_name", organism_col))
  data$org_group <- ifelse(is.na(data$org_group), "Other", data$org_group)
  group_counts <- table(data$org_group)
  message("Organism group distribution:")
  print(group_counts)
  return(data)
}


# ---------------------------------------------------------------------------
# New functions (Layer 4a)
# ---------------------------------------------------------------------------

#' Normalize sp/spp Variants
#'
#' Internal helper: normalizes "sp", "sp.", "spp", "spp." to "spp." and
#' handles common abbreviations like "E.coli" -> "Escherichia coli".
#'
#' @param x Character vector of organism names.
#' @return Character vector with normalized variants.
#' @keywords internal
#' @noRd
prep_normalize_sp_variants <- function(x) {
  x <- gsub("\\bsp\\.\\s*$",   "spp.", x, ignore.case = TRUE)
  x <- gsub("\\bsp\\s+",       "spp. ", x, ignore.case = TRUE)
  x <- gsub("\\bsp\\b",        "spp.", x, ignore.case = TRUE)
  x <- gsub("\\bspp\\b(?!\\.)", "spp.", x, ignore.case = TRUE, perl = TRUE)

  # Common abbreviation expansions
  x <- gsub("\\bE\\.?\\s*coli\\b",              "Escherichia coli",       x, ignore.case = TRUE)
  x <- gsub("\\bK\\.?\\s*pneumoniae\\b",         "Klebsiella pneumoniae",  x, ignore.case = TRUE)
  x <- gsub("\\bP\\.?\\s*aeruginosa\\b",         "Pseudomonas aeruginosa", x, ignore.case = TRUE)
  x <- gsub("\\bA\\.?\\s*baumannii\\b",          "Acinetobacter baumannii",x, ignore.case = TRUE)
  x <- gsub("\\bS\\.?\\s*aureus\\b(?!.*cons)",   "Staphylococcus aureus",  x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("non[- ]?ferm[ea]nt[ia]ng",          "Non-fermenting",         x, ignore.case = TRUE)

  x
}


#' Flag Organisms Unmatched in Reference
#'
#' Flags organisms that could not be matched to the reference organisms.csv.
#' Adds a logical column \code{is_organism_unmatched}.
#'
#' @param data Data frame with organism_normalized column.
#' @param organism_col Character. Standard organism column. Default "organism_normalized".
#'
#' @return Data frame with \code{is_organism_unmatched} column.
#' @export
prep_flag_organism_unmatched <- function(data, organism_col = "organism_normalized") {
  if (!organism_col %in% names(data)) {
    stop(sprintf("Column '%s' not found.", organism_col))
  }

  csv_path <- find_extdata_file("organisms.csv")
  if (csv_path == "" || !file.exists(csv_path)) {
    warning("organisms.csv not found. Cannot flag unmatched organisms.")
    data$is_organism_unmatched <- NA
    return(data)
  }

  org_ref      <- readr::read_csv(csv_path, show_col_types = FALSE)
  ref_names    <- tolower(trimws(org_ref$organism_name))
  norm_vals    <- tolower(trimws(as.character(data[[organism_col]])))
  has_value    <- !is.na(norm_vals) & nzchar(norm_vals)

  data$is_organism_unmatched <- has_value & !norm_vals %in% ref_names

  n_unmatched <- sum(data$is_organism_unmatched, na.rm = TRUE)
  if (n_unmatched > 0) {
    message(sprintf("[prep_flag_organism_unmatched] %d organism(s) not in reference.", n_unmatched))
    top_unmatched <- utils::head(sort(table(
      data[[organism_col]][data$is_organism_unmatched]), decreasing = TRUE), 5)
    message("Top unmatched:")
    print(top_unmatched)
  } else {
    message("[prep_flag_organism_unmatched] All organisms matched to reference.")
  }

  return(data)
}


# ---------------------------------------------------------------------------
# Moved function: get_organism_taxonomy (from prep_clean_and_standardize.R)
# ---------------------------------------------------------------------------

#' Get Organism Taxonomy Mapping
#'
#' Reads organism taxonomy from inst/extdata/organisms.csv and returns
#' a data frame mapping organism names to organism groups.
#'
#' @return Data frame with columns organism_name and org_group.
#' @keywords internal
get_organism_taxonomy <- function() {
  file_path <- find_extdata_file("organisms.csv")
  if (file_path == "") {
    warning("organisms.csv not found. Returning empty taxonomy.")
    return(data.frame(organism_name = character(), org_group = character(),
                      stringsAsFactors = FALSE))
  }
  taxonomy <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  if ("organism_group" %in% names(taxonomy) && !"org_group" %in% names(taxonomy))
    names(taxonomy)[names(taxonomy) == "organism_group"] <- "org_group"
  taxonomy[, c("organism_name", "org_group")]
}
