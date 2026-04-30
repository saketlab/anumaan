# prep_standardize_specimens.R
# Layer 4c: Specimen type standardization
#
# Functions moved here from prep_clean_and_standardize.R:
#   - prep_standardize_specimens
#
# Reference data: inst/extdata/sample_type_classification.csv


#' Normalize Specimen/Sample Type
#'
#' Normalizes specimen/sample type names and adds sample_category and
#' sterile_classification from the reference CSV file. Includes rule-based text
#' cleaning and fuzzy matching for common shorthand and minor misspellings.
#'
#' @param data Data frame with specimen column.
#' @param specimen_col Character. Specimen column name. Default "specimen_type".
#' @param add_categories Logical. Add sample_category and sterile_classification.
#'   Default TRUE.
#'
#' @return Data frame with specimen_normalized, sample_category, and
#'   sterile_classification columns.
#' @export
prep_standardize_specimens <- function(data,
                                       specimen_col    = "specimen_type",
                                       add_categories  = TRUE) {
  if (!specimen_col %in% names(data)) {
    warning(sprintf("Column '%s' not found. Skipping specimen normalization.", specimen_col))
    return(data)
  }

  csv_path <- find_extdata_file("sample_type_classification.csv")

  if (csv_path == "" || !file.exists(csv_path)) {
    warning("sample_type_classification.csv not found. Specimen normalization skipped.")
    data$specimen_normalized <- data[[specimen_col]]
    if (add_categories) {
      data$sample_category      <- NA_character_
      data$sterile_classification <- NA_character_
    }
    return(data)
  }

  specimen_ref <- readr::read_csv(csv_path, show_col_types = FALSE)

  data$temp_spec_input <- tolower(trimws(data[[specimen_col]]))
  data$temp_spec_clean <- data$temp_spec_input

  # String cleaning
  data$temp_spec_clean <- gsub("\\s*culture\\s*/\\s*sensitivity.*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*culture.*$",                       "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*c/s.*$",                           "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*sensitivity.*$",                   "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*\\(.*?\\)\\s*",                    " ", data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\s*\\[.*?\\]\\s*",                    " ", data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\s*from\\s+.*$",                      "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*(catheter|foley|clean catch|voided).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*-\\s*(catheter|foley|midstream|peripheral|central).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("^(midstream|peripheral|central|clean|voided|fresh)\\s+", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s*-\\s*(specimen|sample|report|test|site).*$", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\b(specimen|sample|site|swab from|aspirate from|tip|fluid from)\\b", "", data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\s+", " ", data$temp_spec_clean)
  data$temp_spec_clean <- trimws(data$temp_spec_clean)
  data$temp_spec_clean <- gsub("\\bet\\b",              "endotracheal",           data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\beta\\b",             "endotracheal aspirate",  data$temp_spec_clean, ignore.case = TRUE)
  data$temp_spec_clean <- gsub("\\bbal\\s+fluid\\b",    "bal",                    data$temp_spec_clean, ignore.case = TRUE)

  core_types <- c("blood", "urine", "csf", "sputum", "stool", "pus", "bile",
                  "faeces", "wound", "throat", "nasal", "vaginal", "cervical",
                  "bal", "eta", "endotracheal")
  for (sp in core_types) {
    idx <- grepl(paste0("\\b", sp, "\\b"), data$temp_spec_clean)
    if (any(idx))
      data$temp_spec_clean[idx] <- gsub(paste0(".*\\b(", sp, ")\\b.*"), "\\1", data$temp_spec_clean[idx])
  }

  data$temp_spec_clean <- trimws(data$temp_spec_clean)

  specimen_ref$temp_ref_lower <- tolower(trimws(specimen_ref$specimen_name))

  extract_spec_keywords <- function(name) {
    if (is.na(name) || name == "") return(character(0))
    cleaned <- tolower(gsub("[^a-z0-9\\s/]", " ", name))
    words   <- unlist(strsplit(cleaned, "[\\s/]+"))
    filler  <- c("a", "an", "the", "of", "from", "to", "and", "or", "is", "in", "on")
    words   <- words[!words %in% filler & nchar(words) >= 2]
    unique(words)
  }

  specimen_ref$keywords <- lapply(specimen_ref$temp_ref_lower, extract_spec_keywords)

  data$specimen_normalized    <- NA_character_
  data$sample_category        <- NA_character_
  data$sterile_classification <- NA_character_

  unique_inputs <- unique(data$temp_spec_clean[!is.na(data$temp_spec_clean) & data$temp_spec_clean != ""])
  specimen_map  <- setNames(rep(NA_integer_, length(unique_inputs)), unique_inputs)

  for (input_spec in unique_inputs) {
    exact_match <- which(specimen_ref$temp_ref_lower == input_spec)
    if (length(exact_match) > 0) { specimen_map[input_spec] <- exact_match[1]; next }

    input_keywords <- extract_spec_keywords(input_spec)
    if (length(input_keywords) == 0) next

    scores <- sapply(seq_len(nrow(specimen_ref)), function(i) {
      ref_keywords <- specimen_ref$keywords[[i]]
      if (length(ref_keywords) == 0) return(0)
      overlap <- sum(input_keywords %in% ref_keywords)
      if (grepl(specimen_ref$temp_ref_lower[i], input_spec, fixed = TRUE) ||
            grepl(input_spec, specimen_ref$temp_ref_lower[i], fixed = TRUE))
        overlap <- overlap + 2
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) return(0)
      overlap / total_unique
    })

    max_score <- max(scores)
    if (max_score > 0.3) {
      specimen_map[input_spec] <- which(scores == max_score)[1]
    } else {
      distances    <- sapply(specimen_ref$temp_ref_lower, function(ref)
        adist(input_spec, ref, ignore.case = TRUE)[1, 1])
      min_dist_idx <- which.min(distances)
      if (distances[min_dist_idx] <= 2) specimen_map[input_spec] <- min_dist_idx
    }
  }

  valid_mappings <- specimen_map[!is.na(specimen_map)]
  if (length(valid_mappings) > 0) {
    lookup_df <- data.frame(
      temp_spec_clean      = names(valid_mappings),
      ref_idx              = as.integer(valid_mappings),
      stringsAsFactors     = FALSE
    )
    lookup_df$specimen_normalized    <- specimen_ref$specimen_name[lookup_df$ref_idx]
    lookup_df$sample_category        <- specimen_ref$sample_category[lookup_df$ref_idx]
    lookup_df$sterile_classification <- specimen_ref$sterile_classification[lookup_df$ref_idx]
    lookup_df$ref_idx <- NULL

    data$specimen_normalized    <- NULL
    data$sample_category        <- NULL
    data$sterile_classification <- NULL

    data <- dplyr::left_join(data, lookup_df, by = "temp_spec_clean")
  }

  data$temp_spec_input <- NULL
  data$temp_spec_clean <- NULL

  n_matched <- sum(!is.na(data$specimen_normalized))
  message(sprintf("Specimen normalization: %d/%d matched (%.1f%%)",
                  n_matched, nrow(data), 100 * n_matched / nrow(data)))

  if (add_categories) {
    n_sterile    <- sum(data$sterile_classification == "Sterile site",     na.rm = TRUE)
    n_nonsterile <- sum(data$sterile_classification == "Non-sterile site", na.rm = TRUE)
    message(sprintf("Sterile classification: %d sterile, %d non-sterile", n_sterile, n_nonsterile))
  }

  return(data)
}
