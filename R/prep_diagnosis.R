# prep_diagnosis.R
# Layer: Optional diagnosis text -> ICD -> syndrome pipeline
#
# Functions:
#   prep_diagnosis_text()         - prepare and combine diagnosis text columns
#   prep_map_diagnosis_to_icd()   - map text to ICD candidates (exact/fuzzy/python_embedding)
#   prep_map_icd_to_syndrome()    - map ICD code to infectious syndrome
#   prep_assign_patient_syndrome()- assign one syndrome per patient/event via hierarchy
#   infer_patient_syndrome_wide() - wide-format entry point: pivot then assign syndrome
#   load_icd_reference()          - internal: load icd10_who.csv
#   load_syndrome_hierarchy()     - internal: load infectious_syndrome_hierarchy.csv
#
# The python_embedding method requires reticulate + Python alethia package.
# All other methods work without Python.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Normalize syndrome name for case-insensitive matching
#' @keywords internal
#' @noRd
normalize_syndrome_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\*+", "", x)
  x <- trimws(gsub("\\s+", " ", x))
  tolower(x)
}


# ---------------------------------------------------------------------------
# Internal reference loaders
# ---------------------------------------------------------------------------

#' Load ICD-10 Reference Table
#'
#' Loads \code{icd10_who.csv} from \code{inst/extdata} or a user-supplied path.
#'
#' @param path Character. Optional path to a custom CSV. If NULL, uses the
#'   shipped \code{inst/extdata/icd10_who.csv}.
#' @return Data frame with ICD-10 reference data.
#' @keywords internal
#' @noRd
load_icd_reference <- function(path = NULL) {
  if (!is.null(path)) {
    if (!file.exists(path))
      stop(sprintf("[load_icd_reference] File not found: %s", path))
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }

  pkg_path <- find_extdata_file("icd10_who.csv")
  if (pkg_path == "" || !file.exists(pkg_path))
    stop("[load_icd_reference] icd10_who.csv not found in inst/extdata.")

  utils::read.csv(pkg_path, stringsAsFactors = FALSE)
}


#' Load Infectious Syndrome Hierarchy
#'
#' Loads \code{infectious_syndrome_hierarchy.csv} and enriches it with
#' \code{syndrome_norm} (lowercased, whitespace-squeezed) and \code{aggregate_to}
#' (maps "Other unspecified respiratory site infections" -> "Lower respiratory
#' infections") columns used for case-insensitive matching and syndrome collapsing.
#'
#' @param path Character. Optional path to a custom CSV.
#' @return Data frame with syndrome hierarchy, sorted by rank.
#' @keywords internal
#' @noRd
load_syndrome_hierarchy <- function(path = NULL) {
  if (!is.null(path)) {
    if (!file.exists(path))
      stop(sprintf("[load_syndrome_hierarchy] File not found: %s", path))
    h <- utils::read.csv(path, stringsAsFactors = FALSE)
  } else {
    pkg_path <- find_extdata_file("infectious_syndrome_hierarchy.csv")
    if (pkg_path == "" || !file.exists(pkg_path))
      stop("[load_syndrome_hierarchy] infectious_syndrome_hierarchy.csv not found in inst/extdata.")
    h <- utils::read.csv(pkg_path, stringsAsFactors = FALSE)
  }

  required_cols <- c("rank", "infectious_syndrome", "contributes_to_amr_burden")
  missing_cols  <- setdiff(required_cols, names(h))
  if (length(missing_cols) > 0)
    stop(sprintf("[load_syndrome_hierarchy] Missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))

  h$syndrome_norm <- normalize_syndrome_name(h$infectious_syndrome)
  resp_norm       <- normalize_syndrome_name("Other unspecified respiratory site infections")
  h$aggregate_to  <- ifelse(h$syndrome_norm == resp_norm,
                            "Lower respiratory infections", NA_character_)
  h[order(h$rank), ]
}


# ---------------------------------------------------------------------------
# Layer 1: Prepare diagnosis text
# ---------------------------------------------------------------------------

#' Prepare Diagnosis Text
#'
#' Combines a primary diagnosis column with an optional fallback column,
#' trims whitespace, and collapses repeated separators. Optionally appends
#' the organism name to the text (as the notebook does for context-aware
#' embedding).
#'
#' @param data Data frame.
#' @param diagnosis_col Character. Primary diagnosis column.
#' @param fallback_col Character or NULL. Used when \code{diagnosis_col} is NA.
#'   Default NULL.
#' @param output_col Character. Name for the output text column.
#'   Default \code{"diagnosis_text"}.
#' @param include_organism Logical. Append organism name to the text.
#'   Default FALSE.
#' @param organism_col Character. Organism column. Used only when
#'   \code{include_organism = TRUE}. Default \code{"organism_normalized"}.
#' @param keep_original Logical. Keep the original diagnosis columns.
#'   Default TRUE.
#'
#' @return Data frame with \code{output_col} added.
#' @export
prep_diagnosis_text <- function(data,
                                diagnosis_col,
                                fallback_col      = NULL,
                                output_col        = "diagnosis_text",
                                include_organism  = FALSE,
                                organism_col      = "organism_normalized",
                                keep_original     = TRUE) {
  if (!diagnosis_col %in% names(data))
    stop(sprintf("[prep_diagnosis_text] Column '%s' not found.", diagnosis_col))

  if (!is.null(fallback_col) && !fallback_col %in% names(data))
    stop(sprintf("[prep_diagnosis_text] Fallback column '%s' not found.", fallback_col))

  if (include_organism && !organism_col %in% names(data))
    stop(sprintf("[prep_diagnosis_text] Organism column '%s' not found.", organism_col))

  # Combine primary + fallback
  text <- as.character(data[[diagnosis_col]])
  text <- trimws(text)
  text[text %in% c("", "NA")] <- NA_character_

  if (!is.null(fallback_col)) {
    fallback <- trimws(as.character(data[[fallback_col]]))
    fallback[fallback %in% c("", "NA")] <- NA_character_
    text <- dplyr::coalesce(text, fallback)
  }

  # Clean repeated separators (", ," / ";;" / "  " etc.)
  text <- gsub(",\\s*,", ",", text)
  text <- gsub(";\\s*;", ";", text)
  text <- gsub("\\s{2,}", " ", text)
  text <- trimws(text)

  # Optionally append organism
  if (include_organism) {
    org <- trimws(as.character(data[[organism_col]]))
    org[org %in% c("", "NA")] <- NA_character_
    text <- dplyr::case_when(
      !is.na(text) & !is.na(org) ~ paste0(text, " | organism: ", org),
      !is.na(text)               ~ text,
      TRUE                       ~ NA_character_
    )
  }

  data[[output_col]] <- text

  n_missing <- sum(is.na(text))
  message(sprintf(
    "[prep_diagnosis_text] Created '%s': %d rows, %d missing (%.1f%%)",
    output_col, nrow(data), n_missing, 100 * n_missing / nrow(data)
  ))

  data
}


# ---------------------------------------------------------------------------
# Layer 2: Map diagnosis text to ICD candidates
# ---------------------------------------------------------------------------

#' Map Diagnosis Text to ICD Candidates
#'
#' Maps free-text diagnosis strings to ICD-10 candidate descriptions using
#' one of three methods:
#'
#' \describe{
#'   \item{\code{"exact"}}{Case-insensitive exact string match against ICD
#'     descriptions. Fast, high precision, low recall.}
#'   \item{\code{"fuzzy"}}{String distance matching via \pkg{stringdist}.
#'     Handles typos and minor variations. Requires \pkg{stringdist}.}
#'   \item{\code{"python_embedding"}}{Semantic embedding similarity using the
#'     Python \code{alethia} package via \pkg{reticulate}. Highest recall.
#'     Requires Python, \code{alethia}, and a sentence-transformers model.}
#' }
#'
#' @param data Data frame.
#' @param text_col Character. Column containing prepared diagnosis text
#'   (output of \code{prep_diagnosis_text()}). Default \code{"diagnosis_text"}.
#' @param reference Data frame or NULL. ICD-10 reference table. If NULL,
#'   loads from \code{inst/extdata/icd10_who.csv}.
#' @param method Character. One of \code{"exact"}, \code{"fuzzy"},
#'   \code{"python_embedding"}. Default \code{"exact"}.
#' @param icd_desc_col Character. Column in \code{reference} to match against.
#'   Default \code{"description3"} (most concise ICD labels).
#' @param icd_code_col Character. Column in \code{reference} holding ICD codes.
#'   Default \code{"icd_code_who_eq"}.
#' @param top_k Integer. Maximum ICD candidates to return per input string.
#'   Default 5. Ignored for \code{"exact"} (returns all exact matches).
#' @param threshold Numeric. Minimum similarity score (0-1) to retain a
#'   candidate. Default 0.0 (keep all). For \code{"fuzzy"}: similarity is
#'   \code{1 - normalised_distance}.
#' @param model Character. Sentence-transformers model name. Used only for
#'   \code{"python_embedding"}. Default \code{"FremyCompany/BioLORD-2023"}.
#' @param id_col Character or NULL. Identifier column to carry through into
#'   the output. Default NULL (output contains only match columns).
#'
#' @return Long data frame with columns:
#'   \code{diagnosis_text}, \code{icd_prediction}, \code{icd_code},
#'   \code{icd_score}, \code{icd_rank}, \code{icd_method}.
#'   If \code{id_col} is supplied, it is included as the first column.
#'
#' @export
prep_map_diagnosis_to_icd <- function(data,
                                      text_col     = "diagnosis_text",
                                      reference    = NULL,
                                      method       = c("exact", "fuzzy", "python_embedding"),
                                      icd_desc_col = "description3",
                                      icd_code_col = "icd_code_who_eq",
                                      top_k        = 5L,
                                      threshold    = 0.0,
                                      model        = "FremyCompany/BioLORD-2023",
                                      id_col       = NULL) {
  method <- match.arg(method)

  if (!text_col %in% names(data))
    stop(sprintf("[prep_map_diagnosis_to_icd] Column '%s' not found.", text_col))

  if (!is.null(id_col) && !id_col %in% names(data))
    stop(sprintf("[prep_map_diagnosis_to_icd] id_col '%s' not found.", id_col))

  if (is.null(reference))
    reference <- load_icd_reference()

  if (!icd_desc_col %in% names(reference))
    stop(sprintf("[prep_map_diagnosis_to_icd] ICD description column '%s' not in reference.", icd_desc_col))

  if (!icd_code_col %in% names(reference))
    stop(sprintf("[prep_map_diagnosis_to_icd] ICD code column '%s' not in reference.", icd_code_col))

  # Work on unique non-NA diagnosis texts
  all_texts <- data[[text_col]]
  unique_texts <- unique(all_texts[!is.na(all_texts)])
  message(sprintf(
    "[prep_map_diagnosis_to_icd] Mapping %d unique texts using method '%s'...",
    length(unique_texts), method
  ))

  ref_desc  <- reference[[icd_desc_col]]
  ref_codes <- reference[[icd_code_col]]

  results <- switch(method,

    exact = {
      lapply(unique_texts, function(txt) {
        hits <- which(tolower(ref_desc) == tolower(txt))
        if (length(hits) == 0L) return(NULL)
        hits <- utils::head(hits, top_k)
        data.frame(
          diagnosis_text = txt,
          icd_prediction = ref_desc[hits],
          icd_code       = ref_codes[hits],
          icd_score      = 1.0,
          icd_rank       = seq_along(hits),
          icd_method     = "exact",
          stringsAsFactors = FALSE
        )
      })
    },

    fuzzy = {
      if (!requireNamespace("stringdist", quietly = TRUE))
        stop("[prep_map_diagnosis_to_icd] Package 'stringdist' required for method='fuzzy'. Install it with install.packages('stringdist').")

      unique_refs <- unique(ref_desc)
      lapply(unique_texts, function(txt) {
        dists <- stringdist::stringdist(tolower(txt), tolower(unique_refs), method = "jw")
        scores <- 1 - dists
        keep <- order(scores, decreasing = TRUE)[seq_len(min(top_k, length(scores)))]
        keep <- keep[scores[keep] >= threshold]
        if (length(keep) == 0L) return(NULL)

        matched_desc  <- unique_refs[keep]
        matched_codes <- ref_codes[match(matched_desc, ref_desc)]
        data.frame(
          diagnosis_text = txt,
          icd_prediction = matched_desc,
          icd_code       = matched_codes,
          icd_score      = scores[keep],
          icd_rank       = seq_along(keep),
          icd_method     = "fuzzy_jw",
          stringsAsFactors = FALSE
        )
      })
    },

    python_embedding = {
      if (!requireNamespace("reticulate", quietly = TRUE))
        stop("[prep_map_diagnosis_to_icd] Package 'reticulate' required for method='python_embedding'.")

      alethia_available <- tryCatch({
        reticulate::import("alethia")
        TRUE
      }, error = function(e) FALSE)

      if (!alethia_available)
        stop("[prep_map_diagnosis_to_icd] Python package 'alethia' not found. Install it in the active Python environment.")

      alethia  <- reticulate::import("alethia")
      ref_uniq <- unique(ref_desc)

      result_df <- alethia$alethia(
        dirty_entries     = as.list(unique_texts),
        reference_entries = as.list(ref_uniq),
        model             = model,
        threshold         = threshold
      )
      result_df <- reticulate::py_to_r(result_df)

      # alethia returns one row per match above threshold (may be multiple per input)
      # Columns: given_entity, alethia_prediction, alethia_score, alethia_method, alethia_backend
      result_df <- result_df %>%
        dplyr::rename(
          diagnosis_text = given_entity,
          icd_prediction = alethia_prediction,
          icd_score      = alethia_score
        ) %>%
        dplyr::group_by(diagnosis_text) %>%
        dplyr::arrange(dplyr::desc(icd_score), .by_group = TRUE) %>%
        dplyr::mutate(icd_rank = dplyr::row_number()) %>%
        dplyr::filter(icd_rank <= top_k) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          icd_code   = ref_codes[match(icd_prediction, ref_desc)],
          icd_method = paste0("embedding_", model)
        ) %>%
        dplyr::select(diagnosis_text, icd_prediction, icd_code,
                      icd_score, icd_rank, icd_method)

      list(result_df)
    }
  )

  matched <- dplyr::bind_rows(results)

  if (nrow(matched) == 0L) {
    warning("[prep_map_diagnosis_to_icd] No matches found.")
    return(matched)
  }

  # Join id_col back if requested (join on diagnosis_text)
  if (!is.null(id_col)) {
    id_map <- data %>%
      dplyr::select(dplyr::all_of(c(id_col, text_col))) %>%
      dplyr::distinct()
    matched <- dplyr::left_join(matched, id_map,
                                by = stats::setNames(text_col, "diagnosis_text"))
    matched <- dplyr::relocate(matched, dplyr::all_of(id_col))
  }

  message(sprintf(
    "[prep_map_diagnosis_to_icd] Returned %d candidate rows for %d unique texts.",
    nrow(matched), dplyr::n_distinct(matched$diagnosis_text)
  ))

  matched
}


# ---------------------------------------------------------------------------
# Layer 3: Map ICD code to syndrome
# ---------------------------------------------------------------------------

#' Map ICD Codes to Infectious Syndromes
#'
#' Joins an ICD-to-syndrome reference table onto the data. The reference must
#' have at minimum two columns: an ICD code/description column and a syndrome
#' column. Typically this is a manually curated mapping CSV.
#'
#' @param data Data frame (output of \code{prep_map_diagnosis_to_icd()} or any
#'   data frame with an ICD column).
#' @param icd_col Character. Column in \code{data} holding ICD codes or
#'   descriptions to match on. Default \code{"icd_code"}.
#' @param icd_to_syndrome_ref Data frame. Reference table with at minimum
#'   \code{icd_col} values and a syndrome column.
#' @param ref_icd_col Character. The ICD column name in \code{icd_to_syndrome_ref}.
#'   Default \code{"icd_code"}.
#' @param ref_syndrome_col Character. The syndrome column name in the reference.
#'   Default \code{"infectious_syndrome"}.
#' @param output_col Character. Name for the syndrome column added to \code{data}.
#'   Default \code{"syndrome"}.
#' @param unmatched_label Character. Label for rows with no syndrome match.
#'   Default \code{NA_character_}.
#'
#' @return Data frame with \code{output_col} added.
#' @export
prep_map_icd_to_syndrome <- function(data,
                                     icd_col           = "icd_code",
                                     icd_to_syndrome_ref,
                                     ref_icd_col       = "icd_code",
                                     ref_syndrome_col  = "infectious_syndrome",
                                     output_col        = "syndrome",
                                     unmatched_label   = NA_character_) {
  if (!icd_col %in% names(data))
    stop(sprintf("[prep_map_icd_to_syndrome] Column '%s' not found in data.", icd_col))

  if (!ref_icd_col %in% names(icd_to_syndrome_ref))
    stop(sprintf("[prep_map_icd_to_syndrome] Column '%s' not found in reference.", ref_icd_col))

  if (!ref_syndrome_col %in% names(icd_to_syndrome_ref))
    stop(sprintf("[prep_map_icd_to_syndrome] Column '%s' not found in reference.", ref_syndrome_col))

  ref_slim <- icd_to_syndrome_ref %>%
    dplyr::select(dplyr::all_of(c(ref_icd_col, ref_syndrome_col))) %>%
    dplyr::distinct() %>%
    dplyr::rename(!!output_col := dplyr::all_of(ref_syndrome_col))

  n_before <- nrow(data)
  data <- dplyr::left_join(
    data, ref_slim,
    by = stats::setNames(ref_icd_col, icd_col)
  )

  if (!is.na(unmatched_label)) {
    data[[output_col]][is.na(data[[output_col]])] <- unmatched_label
  }

  n_matched   <- sum(!is.na(data[[output_col]]))
  n_unmatched <- sum(is.na(data[[output_col]]))
  message(sprintf(
    "[prep_map_icd_to_syndrome] Syndrome assigned: %d matched, %d unmatched.",
    n_matched, n_unmatched
  ))

  data
}


# ---------------------------------------------------------------------------
# Layer 4: Assign one syndrome per patient / event
# ---------------------------------------------------------------------------

#' Assign One Syndrome Per Patient or Event
#'
#' Selects the highest-priority syndrome per patient using the infectious
#' syndrome hierarchy. Matching is case-insensitive (uses normalized names).
#'
#' Priority order per patient:
#' \enumerate{
#'   \item Hierarchy rank (lower = higher priority)
#'   \item ICD match score (higher = higher priority), when \code{score_col} present
#'   \item Alphabetical tie-break
#' }
#'
#' This function absorbs the logic of the former \code{infer_patient_syndrome_long()}:
#' use \code{collapse_unspecified_respiratory = TRUE} and
#' \code{keep_only_burden_syndromes = TRUE} to replicate that behaviour.
#'
#' @param data Data frame with one row per syndrome candidate.
#' @param patient_col Character. Patient or event ID column. Default \code{"patient_id"}.
#' @param syndrome_col Character. Syndrome column. Default \code{"syndrome"}.
#' @param score_col Character or NULL. ICD match score column for secondary priority.
#'   Default \code{"icd_score"}.
#' @param hierarchy_ref Data frame or NULL. Syndrome hierarchy (from
#'   \code{load_syndrome_hierarchy()}). If NULL, loads automatically.
#' @param hierarchy_syndrome_col Character. Syndrome column in hierarchy. Default \code{"infectious_syndrome"}.
#' @param hierarchy_rank_col Character. Rank column in hierarchy. Default \code{"rank"}.
#' @param keep_all_candidates Logical. TRUE -> return all rows with \code{syndrome_selected}
#'   flag; FALSE -> return one row per \code{patient_col} (selected only). Default TRUE.
#' @param collapse_unspecified_respiratory Logical. Map "Other unspecified respiratory
#'   site infections" -> "Lower respiratory infections" before selection. Default FALSE.
#' @param keep_only_burden_syndromes Logical. Retain only syndromes where
#'   \code{contributes_to_amr_burden == TRUE}. Default FALSE.
#'
#' @return Data frame. When \code{keep_all_candidates = FALSE}, includes
#'   \code{contributes_to_amr_burden} column from the hierarchy.
#' @export
prep_assign_patient_syndrome <- function(data,
                                         patient_col                      = "patient_id",
                                         syndrome_col                     = "syndrome",
                                         score_col                        = "icd_score",
                                         hierarchy_ref                    = NULL,
                                         hierarchy_syndrome_col           = "infectious_syndrome",
                                         hierarchy_rank_col               = "rank",
                                         keep_all_candidates              = TRUE,
                                         collapse_unspecified_respiratory = FALSE,
                                         keep_only_burden_syndromes       = FALSE) {
  if (!patient_col %in% names(data))
    stop(sprintf("[prep_assign_patient_syndrome] Column '%s' not found.", patient_col))
  if (!syndrome_col %in% names(data))
    stop(sprintf("[prep_assign_patient_syndrome] Column '%s' not found.", syndrome_col))

  use_score <- !is.null(score_col) && score_col %in% names(data)

  if (is.null(hierarchy_ref))
    hierarchy_ref <- load_syndrome_hierarchy()

  if (!hierarchy_syndrome_col %in% names(hierarchy_ref))
    stop(sprintf("[prep_assign_patient_syndrome] Hierarchy column '%s' not found.", hierarchy_syndrome_col))
  if (!hierarchy_rank_col %in% names(hierarchy_ref))
    stop(sprintf("[prep_assign_patient_syndrome] Rank column '%s' not found.", hierarchy_rank_col))

  # Ensure hierarchy has normalized names (load_syndrome_hierarchy adds these,
  # but a user-supplied ref might not)
  if (!"syndrome_norm" %in% names(hierarchy_ref))
    hierarchy_ref$syndrome_norm <- normalize_syndrome_name(hierarchy_ref[[hierarchy_syndrome_col]])
  if (!"aggregate_to" %in% names(hierarchy_ref))
    hierarchy_ref$aggregate_to <- NA_character_

  # Build named-vector lookups indexed by normalized syndrome name
  rank_map      <- stats::setNames(hierarchy_ref[[hierarchy_rank_col]], hierarchy_ref$syndrome_norm)
  aggregate_map <- stats::setNames(hierarchy_ref$aggregate_to,          hierarchy_ref$syndrome_norm)
  burden_map    <- if ("contributes_to_amr_burden" %in% names(hierarchy_ref))
    stats::setNames(hierarchy_ref$contributes_to_amr_burden, hierarchy_ref$syndrome_norm) else NULL

  max_rank <- max(rank_map, na.rm = TRUE) + 1L

  # Normalize input syndrome values for lookup
  data$.norm_syn <- normalize_syndrome_name(data[[syndrome_col]])

  # Optionally collapse "Other unspecified respiratory site infections"
  if (collapse_unspecified_respiratory) {
    agg_target <- aggregate_map[data$.norm_syn]
    collapse_rows <- !is.na(agg_target)
    if (any(collapse_rows, na.rm = TRUE)) {
      data[[syndrome_col]][collapse_rows] <- agg_target[collapse_rows]
      data$.norm_syn <- normalize_syndrome_name(data[[syndrome_col]])
      message(sprintf("[prep_assign_patient_syndrome] %d row(s) collapsed to aggregate syndrome.",
                      sum(collapse_rows)))
    }
  }

  # Assign rank and optional burden flag
  data$.syndrome_rank <- dplyr::coalesce(rank_map[data$.norm_syn], max_rank)
  if (!is.null(burden_map))
    data$contributes_to_amr_burden <- burden_map[data$.norm_syn]

  # Optionally restrict to burden syndromes
  if (keep_only_burden_syndromes && !is.null(burden_map)) {
    data <- data[!is.na(data$contributes_to_amr_burden) &
                   as.logical(data$contributes_to_amr_burden) == TRUE, ]
  }

  # Select best syndrome per patient
  if (use_score) {
    data <- data %>%
      dplyr::group_by(!!rlang::sym(patient_col)) %>%
      dplyr::arrange(.syndrome_rank, dplyr::desc(!!rlang::sym(score_col)),
                     !!rlang::sym(syndrome_col), .by_group = TRUE) %>%
      dplyr::mutate(syndrome_selected = dplyr::row_number() == 1L) %>%
      dplyr::ungroup()
  } else {
    data <- data %>%
      dplyr::group_by(!!rlang::sym(patient_col)) %>%
      dplyr::arrange(.syndrome_rank, !!rlang::sym(syndrome_col), .by_group = TRUE) %>%
      dplyr::mutate(syndrome_selected = dplyr::row_number() == 1L) %>%
      dplyr::ungroup()
  }

  data <- data %>% dplyr::select(-.syndrome_rank, -.norm_syn)

  n_assigned <- sum(data$syndrome_selected, na.rm = TRUE)
  message(sprintf("[prep_assign_patient_syndrome] Syndrome selected for %d %s(s).",
                  n_assigned, patient_col))

  if (!keep_all_candidates) {
    data <- data %>%
      dplyr::filter(syndrome_selected) %>%
      dplyr::select(-syndrome_selected)
  }

  data
}


# ---------------------------------------------------------------------------
# Layer 5: Wide-format entry point
# ---------------------------------------------------------------------------

#' Assign Syndrome from Wide-Format Syndrome Flags
#'
#' Converts a wide-format data frame (one column per syndrome, values indicating
#' presence) to long format and assigns the highest-priority syndrome per patient
#' using the infectious syndrome hierarchy.
#'
#' This is a convenience wrapper around \code{prep_assign_patient_syndrome()}.
#' All syndrome-selection logic and parameters (hierarchy, burden filter,
#' respiratory collapsing) are handled there.
#'
#' @param data Data frame in wide format where syndrome columns contain presence
#'   indicators (\code{1}/\code{TRUE}/\code{"Yes"} etc.).
#' @param patient_col Character. Patient ID column. Default \code{"patient_id"}.
#' @param syndrome_cols Character vector. Columns to treat as syndrome flags.
#'   If NULL, all columns except \code{patient_col} are used. Default NULL.
#' @param positive_values Character/logical/numeric vector. Values treated as
#'   "syndrome present". Default covers common representations of TRUE/1/Yes.
#' @param collapse_unspecified_respiratory Logical. Passed to
#'   \code{prep_assign_patient_syndrome()}. Default TRUE.
#' @param keep_only_burden_syndromes Logical. Passed to
#'   \code{prep_assign_patient_syndrome()}. Default FALSE.
#'
#' @return Data frame with one row per patient and the selected syndrome.
#' @export
infer_patient_syndrome_wide <- function(data,
                                        patient_col                      = "patient_id",
                                        syndrome_cols                    = NULL,
                                        positive_values                  = c(1, "1", TRUE, "TRUE",
                                                                             "True", "true",
                                                                             "Yes", "YES", "yes"),
                                        collapse_unspecified_respiratory = TRUE,
                                        keep_only_burden_syndromes       = FALSE) {
  if (!patient_col %in% names(data))
    stop(sprintf("[infer_patient_syndrome_wide] Column '%s' not found.", patient_col))

  if (is.null(syndrome_cols))
    syndrome_cols <- setdiff(names(data), patient_col)

  missing_cols <- setdiff(syndrome_cols, names(data))
  if (length(missing_cols) > 0)
    stop(sprintf("[infer_patient_syndrome_wide] Missing syndrome columns: %s",
                 paste(missing_cols, collapse = ", ")))

  long_df <- data %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(syndrome_cols),
      names_to  = "syndrome",
      values_to = ".present"
    ) %>%
    dplyr::filter(as.character(.data$.present) %in% as.character(positive_values)) %>%
    dplyr::select(-.present)

  message(sprintf("[infer_patient_syndrome_wide] %d positive syndrome rows after pivot.",
                  nrow(long_df)))

  prep_assign_patient_syndrome(
    data                             = long_df,
    patient_col                      = patient_col,
    syndrome_col                     = "syndrome",
    score_col                        = NULL,
    keep_all_candidates              = FALSE,
    collapse_unspecified_respiratory = collapse_unspecified_respiratory,
    keep_only_burden_syndromes       = keep_only_burden_syndromes
  )
}
