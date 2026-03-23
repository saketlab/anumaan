# utils.R
# Small utility functions used across modules

#' Largest-Remainder Rounding
#'
#' Rounds a numeric vector so that the individual rounded values sum exactly to
#' a target total.  Uses the largest-remainder method (Hamilton method).
#'
#' @param x Numeric vector to round.
#' @param target Integer target sum. Default is \code{round(sum(x))}.
#'
#' @return Integer vector of the same length as \code{x} whose sum equals
#'   \code{target}.
#' @export
#'
#' @examples
#' round_to_sum(c(3.3, 3.3, 3.4), target = 10)
round_to_sum <- function(x, target = round(sum(x))) {
  floors <- floor(x)
  remainder <- target - sum(floors)
  ranks <- order(x - floors, decreasing = TRUE)
  floors[ranks[seq_len(remainder)]] <- floors[ranks[seq_len(remainder)]] + 1L
  floors
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
