# zzz_utils_internal.R
# Internal utility helpers (not exported)

#' Parse Age Bin Labels to Breaks
#'
#' Parses age-bin labels into numeric breakpoints for interval binning.
#'
#' @param labels Character vector of age-bin labels, such as \code{"<1"},
#'   \code{"1-5"}, or \code{"85+"}.
#'
#' @return A list with \code{breaks} (numeric breakpoints) and \code{labels}
#'   (original labels).
#'
#' @keywords internal
parse_age_bin_labels <- function(labels) {
  clean_labels <- as.character(labels)
  breaks <- numeric(0)

  for (i in seq_along(clean_labels)) {
    label <- clean_labels[i]

    if (is.na(label) || !nzchar(trimws(label))) {
      stop(sprintf("Cannot parse age bin label: %s", label))
    }

    if (grepl("^<\\s*[0-9.]+\\s*$", label)) {
      # Handle "<1" format.
      upper <- as.numeric(gsub("^<\\s*", "", label))
      breaks <- c(breaks, -Inf, upper)
    } else if (grepl("^\\s*[0-9.]+\\s*-\\s*[0-9.]+\\s*$", label)) {
      # Handle "1-5" format.
      parts <- strsplit(gsub("\\s+", "", label), "-", fixed = TRUE)[[1]]
      lower <- as.numeric(parts[1])
      upper <- as.numeric(parts[2])
      breaks <- c(breaks, lower, upper)
    } else if (grepl("^\\s*[0-9.]+\\s*\\+\\s*$", label)) {
      # Handle "85+" format.
      lower <- as.numeric(gsub("\\+\\s*$", "", gsub("\\s+", "", label)))
      breaks <- c(breaks, lower, Inf)
    } else {
      stop(sprintf("Cannot parse age bin label: %s", label))
    }
  }

  list(breaks = unique(breaks), labels = clean_labels)
}
# utils.R
# Small utility functions used across modules

#' Largest-Remainder Rounding
#'
#' Rounds a numeric vector so that the individual rounded values sum exactly to
#' a target total. Uses the largest-remainder method (Hamilton method).
#' Input is coerced with \code{as.numeric()}; non-numeric values become
#' \code{NA}.
#'
#' @param x Numeric vector to round.
#' @param target Integer target sum. Default is
#'   \code{round(sum(x, na.rm = TRUE))}.
#'
#' @return Integer vector of the same length as \code{x} whose sum equals
#'   \code{target} across non-missing entries; missing values are preserved as
#'   \code{NA}.
#' @export
#'
#' @examples
#' round_to_sum(c(3.3, 3.3, 3.4), target = 10)
round_to_sum <- function(x, target = round(sum(x, na.rm = TRUE))) {
  x_num <- suppressWarnings(as.numeric(x))
  out <- rep(NA_integer_, length(x_num))
  idx <- which(!is.na(x_num))

  if (!length(idx)) {
    return(out)
  }

  floors <- floor(x_num[idx])
  fracs <- x_num[idx] - floors
  n <- length(floors)
  remainder <- as.integer(round(target - sum(floors)))

  if (remainder > 0L) {
    full <- remainder %/% n
    extra <- remainder %% n
    floors <- floors + full
    if (extra > 0L) {
      ranks <- order(fracs, decreasing = TRUE)
      floors[ranks[seq_len(extra)]] <- floors[ranks[seq_len(extra)]] + 1L
    }
  } else if (remainder < 0L) {
    take <- -remainder
    full <- take %/% n
    extra <- take %% n
    floors <- floors - full
    if (extra > 0L) {
      ranks <- order(fracs, decreasing = FALSE)
      floors[ranks[seq_len(extra)]] <- floors[ranks[seq_len(extra)]] - 1L
    }
  }

  out[idx] <- as.integer(floors)
  out
}


#' Shorten Antibiotic Class Names
#'
#' Maps long antibiotic class names to common abbreviations used in GBD-style
#' figures.
#'
#' @param x Character vector of antibiotic class names.
#'
#' @return Character vector of shortened names.
#' @export
#'
#' @examples
#' shorten_drug_class(c("Carbapenems", "Third-generation-cephalosporins"))
shorten_drug_class <- function(x) {
  x <- sub("_R$", "", x)
  dplyr::case_when(
    x == "Third-generation-cephalosporins" ~ "3GC",
    x == "Fourth-generation-cephalosporins" ~ "4GC",
    x == "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal" ~ "BL/BLI-AP",
    x == "Beta-lactam/beta-lactamase-inhibitor" ~ "BL/BLI",
    x == "Trimethoprim-sulfamethoxazole" ~ "TMP-SMX",
    x == "Anti-pseudomonal-penicillins" ~ "Anti-pseudo-PCN",
    TRUE ~ x
  )
}
normalize_join <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- gsub("\\s+", " ", x)
  x
}
