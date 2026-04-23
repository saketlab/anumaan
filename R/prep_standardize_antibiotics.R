# prep_standardize_antibiotics.R
# Layer 4b: Antibiotic name standardization + AWaRe classification
#
# Functions moved here from prep_clean_and_standardize.R:
#   - prep_standardize_antibiotics
#   - prep_classify_antibiotic_class
#   - prep_classify_aware
#
# New functions: prep_decode_antibiotic_code


#' Normalize Antibiotic Names
#'
#' Standardizes antibiotic names using fuzzy matching against WHO AWaRe reference.
#'
#' @param data Data frame with antibiotic data.
#' @param antibiotic_col Character. Column name with antibiotic names. Default "antibiotic_name".
#' @param who_table Data frame. WHO AWaRe classification table. If NULL,
#'   loads from inst/extdata/WHO_aware_class.csv.
#' @param add_class Logical. Add antibiotic_class column. Default TRUE.
#' @param add_aware Logical. Add aware_category column. Default TRUE.
#'
#' @return Data frame with antibiotic_normalized, antibiotic_class, aware_category.
#' @export
prep_standardize_antibiotics <- function(data,
                                         antibiotic_col = "antibiotic_name",
                                         who_table      = NULL,
                                         add_class      = TRUE,
                                         add_aware      = TRUE) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", antibiotic_col))
  }

  if (is.null(who_table)) {
    who_path <- find_extdata_file("WHO_aware_class.csv")
    if (who_path == "" || !file.exists(who_path)) {
      warning("WHO_aware_class.csv not found. Antibiotic normalization skipped.")
      data$antibiotic_normalized <- data[[antibiotic_col]]
      if (add_class) data$antibiotic_class  <- NA_character_
      if (add_aware) data$aware_category <- NA_character_
      return(data)
    }
    who_table <- readr::read_csv(who_path, show_col_types = FALSE)
  }

  # Normalise column names in who_table
  if ("Antibiotic" %in% names(who_table) && !"antibiotic_name" %in% names(who_table))
    who_table <- who_table %>% dplyr::rename(antibiotic_name = Antibiotic)
  if ("Class" %in% names(who_table) && !"antibiotic_class" %in% names(who_table))
    who_table <- who_table %>% dplyr::rename(antibiotic_class = Class)
  if ("Category " %in% names(who_table) && !"aware_category" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(aware_category = `Category `)
  } else if ("Category" %in% names(who_table) && !"aware_category" %in% names(who_table)) {
    who_table <- who_table %>% dplyr::rename(aware_category = Category)
  }

  if (!all(c("antibiotic_name", "antibiotic_class") %in% names(who_table)))
    stop("WHO table must have 'antibiotic_name' and 'antibiotic_class' columns")

  n_unique_before <- dplyr::n_distinct(data[[antibiotic_col]], na.rm = TRUE)
  message(sprintf("Normalizing %d unique antibiotic names against WHO reference (%d antibiotics)...",
                  n_unique_before, nrow(who_table)))

  extract_abx_keywords <- function(name) {
    if (is.na(name) || name == "") return(character(0))
    cleaned <- tolower(gsub("[^a-z0-9\\s/-]", " ", name))
    words   <- unlist(strsplit(cleaned, "[\\s/-]+"))
    filler  <- c("iv", "oral", "injection", "tablet", "mg", "strip", "mic", "done", "by")
    words   <- words[!words %in% filler & nchar(words) > 2]
    unique(words)
  }

  who_table$keywords  <- lapply(who_table$antibiotic_name, extract_abx_keywords)
  who_table$ref_lower <- tolower(trimws(who_table$antibiotic_name))

  data$temp_abx_input    <- trimws(as.character(data[[antibiotic_col]]))
  data$antibiotic_normalized <- NA_character_

  unique_inputs  <- unique(data$temp_abx_input[!is.na(data$temp_abx_input) & data$temp_abx_input != ""])
  antibiotic_map <- setNames(rep(NA_character_, length(unique_inputs)), unique_inputs)

  for (input_abx in unique_inputs) {
    input_keywords <- extract_abx_keywords(input_abx)

    if (length(input_keywords) == 0) {
      antibiotic_map[input_abx] <- tolower(input_abx)
      next
    }

    scores <- sapply(seq_len(nrow(who_table)), function(i) {
      ref_keywords <- who_table$keywords[[i]]
      if (length(ref_keywords) == 0) return(0)
      overlap <- sum(input_keywords %in% ref_keywords)
      if (tolower(input_abx) == who_table$ref_lower[i]) return(1000)
      if (grepl(who_table$ref_lower[i], tolower(input_abx), fixed = TRUE) ||
            grepl(tolower(input_abx), who_table$ref_lower[i], fixed = TRUE))
        overlap <- overlap + 2
      total_unique <- length(union(input_keywords, ref_keywords))
      if (total_unique == 0) return(0)
      overlap / total_unique
    })

    max_score <- max(scores)
    if (max_score > 0.3) {
      top_matches <- which(scores == max_score)
      if (length(top_matches) == 1) {
        antibiotic_map[input_abx] <- who_table$ref_lower[top_matches[1]]
      } else {
        input_lower   <- tolower(input_abx)
        tie_distances <- sapply(top_matches, function(idx)
          adist(input_lower, who_table$ref_lower[idx], ignore.case = TRUE)[1, 1])
        antibiotic_map[input_abx] <- who_table$ref_lower[top_matches[which.min(tie_distances)]]
      }
    } else {
      input_lower <- tolower(input_abx)
      distances   <- sapply(who_table$ref_lower, function(ref)
        adist(input_lower, ref, ignore.case = TRUE)[1, 1])
      min_dist_idx <- which.min(distances)
      if (distances[min_dist_idx] <= 3) {
        antibiotic_map[input_abx] <- who_table$ref_lower[min_dist_idx]
      } else {
        antibiotic_map[input_abx] <- tolower(input_abx)
      }
    }
  }

  data$antibiotic_normalized <- antibiotic_map[data$temp_abx_input]
  data$antibiotic_normalized[is.na(data$temp_abx_input) | data$temp_abx_input == ""] <- NA_character_
  data$antibiotic_normalized[!is.na(data$antibiotic_normalized) &
                               trimws(data$antibiotic_normalized) == ""]              <- NA_character_

  if (add_class || add_aware) {
    select_cols <- "ref_lower"
    if (add_class) select_cols <- c(select_cols, "antibiotic_class")
    if (add_aware && "aware_category" %in% names(who_table))
      select_cols <- c(select_cols, "aware_category")

    data <- data %>%
      dplyr::left_join(who_table %>% dplyr::select(dplyr::all_of(select_cols)),
                       by = c("antibiotic_normalized" = "ref_lower"))
  }

  orig_empty <- is.na(data$temp_abx_input) | data$temp_abx_input == ""
  if (add_class && "antibiotic_class" %in% names(data)) {
    data$antibiotic_class[orig_empty] <- NA_character_
    # Special corrections
    orig_lower <- tolower(data$temp_abx_input)
    data$antibiotic_class[!is.na(orig_lower) & orig_lower == "colistin"]         <- "Colistin"
    data$antibiotic_class[!is.na(orig_lower) & grepl("nalidixic", orig_lower)]   <- "Nalidixic acid"
    data$antibiotic_class[!is.na(orig_lower) & grepl("polymyxin.*b", orig_lower)] <- "Polymyxins"
  }
  if (add_aware && "aware_category" %in% names(data))
    data$aware_category[orig_empty] <- NA_character_

  data$temp_abx_input <- NULL

  n_unique_after <- dplyr::n_distinct(data$antibiotic_normalized, na.rm = TRUE)
  message(sprintf("Normalized: %d unique names -> %d", n_unique_before, n_unique_after))
  if (add_class && "antibiotic_class" %in% names(data))
    message(sprintf("With antibiotic_class: %d (%.1f%%)",
                    sum(!is.na(data$antibiotic_class)),
                    100 * sum(!is.na(data$antibiotic_class)) / nrow(data)))
  if (add_aware && "aware_category" %in% names(data))
    message(sprintf("With AWaRe category: %d (%.1f%%)",
                    sum(!is.na(data$aware_category)),
                    100 * sum(!is.na(data$aware_category)) / nrow(data)))

  return(data)
}


#' Classify Antibiotic to WHO Class
#'
#' Maps antibiotic names to WHO antibiotic classes.
#'
#' @param data Data frame.
#' @param antibiotic_col Character. Normalized antibiotic column. Default "antibiotic_normalized".
#' @param who_table Data frame. WHO classification table. If NULL uses built-in.
#'
#' @return Data frame with antibiotic_class column added.
#' @export
prep_classify_antibiotic_class <- function(data,
                                           antibiotic_col = "antibiotic_normalized",
                                           who_table      = NULL) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Antibiotic column '%s' not found", antibiotic_col))
  }

  if (is.null(who_table)) {
    message("Using built-in WHO class mapping")
    data$antibiotic_class <- NA_character_
    data$class_source     <- "needs_who_table"
    warning("WHO table not provided.")
    return(data)
  }

  data <- data %>%
    dplyr::left_join(who_table %>% dplyr::select(Antibiotic, Class),
                     by = stats::setNames("Antibiotic", antibiotic_col)) %>%
    dplyr::rename(antibiotic_class = Class)

  n_classified <- sum(!is.na(data$antibiotic_class))
  message(sprintf("Classified antibiotics: %d/%d (%.1f%%)",
                  n_classified, nrow(data), 100 * n_classified / nrow(data)))
  return(data)
}


#' Classify AWaRe Category
#'
#' Assigns WHO AWaRe (Access, Watch, Reserve) categories to antibiotics.
#'
#' @param data Data frame.
#' @param antibiotic_col Character. Antibiotic column. Default "antibiotic_normalized".
#' @param who_table Data frame. WHO AWaRe table. If NULL, uses built-in.
#'
#' @return Data frame with aware_category column added.
#' @export
prep_classify_aware <- function(data,
                                antibiotic_col = "antibiotic_normalized",
                                who_table      = NULL) {
  if (!antibiotic_col %in% names(data)) {
    stop(sprintf("Antibiotic column '%s' not found", antibiotic_col))
  }

  if (is.null(who_table)) {
    message("Using built-in AWaRe mapping")
    data$aware_category <- NA_character_
    data$aware_source   <- "needs_who_table"
    warning("WHO AWaRe table not provided.")
    return(data)
  }

  data <- data %>%
    dplyr::left_join(who_table %>% dplyr::select(Antibiotic, Category),
                     by = stats::setNames("Antibiotic", antibiotic_col)) %>%
    dplyr::rename(aware_category = Category)

  message(sprintf("Classified AWaRe: %d records", sum(!is.na(data$aware_category))))
  print(table(data$aware_category, useNA = "ifany"))
  return(data)
}


# ---------------------------------------------------------------------------
# New function (Layer 4b)
# ---------------------------------------------------------------------------

#' Decode Antibiotic Short Codes to Full Names
#'
#' Converts antibiotic short codes (e.g. AMK, AMP, AMPSUL) to full antibiotic
#' names using a caller-supplied code map CSV.
#'
#' The code map must have columns \code{code} and \code{antibiotic_name_full}.
#' This function is dataset-agnostic: the code map and the column containing
#' the codes are supplied by the caller (typically from the analysis layer).
#'
#' @param data Data frame with a column of antibiotic codes.
#' @param code_col Character. Column containing antibiotic codes.
#'   Default "antibiotic_name_raw".
#' @param map_path Character. Path to the code map CSV file. Must have columns
#'   \code{code} and \code{antibiotic_name_full}. If NULL, falls back to
#'   \code{inst/extdata/antibiotic_code_map.csv} if it exists.
#' @param output_col Character. Name of the output column for decoded names.
#'   Default "antibiotic_name_std".
#'
#' @return Data frame with \code{output_col} column added.
#' @export
prep_decode_antibiotic_code <- function(data,
                                         code_col   = "antibiotic_name_raw",
                                         map_path   = NULL,
                                         output_col = "antibiotic_name_std") {
  if (!code_col %in% names(data)) {
    warning(sprintf("[prep_decode_antibiotic_code] Column '%s' not found. Skipping.", code_col))
    data[[output_col]] <- NA_character_
    return(data)
  }

  if (is.null(map_path) || map_path == "") {
    map_path <- find_extdata_file("antibiotic_code_map.csv")
  }

  if (is.null(map_path) || map_path == "" || !file.exists(map_path)) {
    warning("[prep_decode_antibiotic_code] Code map CSV not found. Passing raw values through.")
    data[[output_col]] <- data[[code_col]]
    return(data)
  }

  code_map <- readr::read_csv(map_path, show_col_types = FALSE)

  if (!all(c("code", "antibiotic_name_full") %in% names(code_map)))
    stop("[prep_decode_antibiotic_code] Code map must have columns 'code' and 'antibiotic_name_full'.")

  lookup <- stats::setNames(code_map$antibiotic_name_full, toupper(code_map$code))

  raw_upper            <- toupper(trimws(as.character(data[[code_col]])))
  decoded              <- lookup[raw_upper]
  decoded[is.na(decoded)] <- data[[code_col]][is.na(decoded)]
  data[[output_col]]   <- decoded

  n_mapped   <- sum(!is.na(lookup[raw_upper]), na.rm = TRUE)
  n_unmapped <- sum(is.na(lookup[raw_upper]), na.rm = TRUE)
  message(sprintf("[prep_decode_antibiotic_code] %d codes decoded, %d passed through unchanged.",
                  n_mapped, n_unmapped))

  return(data)
}
