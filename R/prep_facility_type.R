#' Assign a facility type / sector column from a user-supplied mapping
#'
#' Adds a column that labels each row's facility with a category defined by the
#' caller (for example "Government" vs "Private", "Urban" vs "Rural", or a
#' region). The function is deliberately generic: it contains no hardcoded
#' facility names, so any dataset can be classified by passing an appropriate
#' \code{mapping}.
#'
#' @param data A data frame.
#' @param facility_col Name of the column holding facility names.
#'   Default \code{"center_name"}.
#' @param mapping A named character vector where names are facility values and
#'   values are the type/sector to assign, e.g.
#'   \code{c(AFMC = "Government", SGRH = "Private")}.
#' @param new_col Name of the column to create. Default \code{"facility_type"}.
#' @param default Value assigned to facilities not found in \code{mapping}.
#'   Default \code{NA_character_}.
#' @param case_insensitive Logical. If \code{TRUE} (default), matching ignores
#'   case and surrounding whitespace, so \code{"AFMC"}, \code{"afmc"} and
#'   \code{" AFMC "} all match.
#'
#' @return \code{data} with \code{new_col} added.
#'
#' @examples
#' df <- data.frame(center_name = c("AFMC", "SGRH", "Unknown"))
#' m  <- c(AFMC = "Government", SGRH = "Private")
#' prep_assign_facility_type(df, mapping = m)
#'
#' @export
prep_assign_facility_type <- function(data,
                                      facility_col     = "center_name",
                                      mapping,
                                      new_col          = "facility_type",
                                      default          = NA_character_,
                                      case_insensitive = TRUE) {
  if (!is.data.frame(data))
    stop("`data` must be a data frame.")
  if (!facility_col %in% names(data))
    stop(sprintf("Column '%s' not found in data.", facility_col))
  if (missing(mapping) || is.null(names(mapping)) || length(mapping) == 0L)
    stop("`mapping` must be a non-empty named character vector (names = facilities).")

  key <- as.character(data[[facility_col]])

  if (isTRUE(case_insensitive)) {
    lk  <- stats::setNames(unname(as.character(mapping)),
                           tolower(trimws(names(mapping))))
    out <- unname(lk[tolower(trimws(key))])
  } else {
    lk  <- stats::setNames(unname(as.character(mapping)), names(mapping))
    out <- unname(lk[key])
  }

  out[is.na(out)] <- default
  data[[new_col]] <- out
  data
}
