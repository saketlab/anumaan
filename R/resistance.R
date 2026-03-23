# resistance.R
# Functions for computing resistance profiles per pathogen.
#
# Workflow (three steps):
#
#   Step 1 -- compute_marginal_resistance()
#     Collapses drug-level data to class level (R if ANY drug in class R),
#     computes marginal resistance per pathogen x class, and flags classes
#     whose marginal resistance is at or below zero_threshold. The flagged
#     list is informational -- downstream profiling can choose to exclude them.
#
#   Step 2 -- compute_pairwise_coresistance()
#     Uses the collapsed class-level data from Step 1 to compute pairwise
#     co-resistance matrices for ALL tested classes per pathogen:
#
#       T[k, c_i, c_j] = # isolates tested for both class c_i and c_j
#       R[k, c_i, c_j] = # isolates resistant to both
#       Prev[k, c_i, c_j] = R / T
#
#   Step 3 -- compute_resistance_profiles()
#     Enumerates all 2^n binary resistance profiles delta in {0,1}^n and
#     estimates profile probabilities via a simplex-constrained weighted
#     least-squares QP (GBD eq. 7.5.1.3), using the marginal and pairwise
#     co-resistance estimates as constraints.
#
# Unit of analysis: isolate (unique value of isolate_col per pathogen).
#
# Class-level resistance rule (GBD):
#   R_{e,k,c} = 1  if the isolate is resistant to ANY drug in class c.
#
# References:
#   Antimicrobial Resistance Collaborators. Lancet. 2022.


# -- Internal helper -----------------------------------------------------------

#' Validate that required columns exist in a data frame
#' @keywords internal
.check_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf(
      "Column(s) not found in data: %s",
      paste(missing, collapse = ", ")
    ))
  }
}


# -- Step 1 --------------------------------------------------------------------

#' Compute Marginal Resistance per Pathogen and Antibiotic Class
#'
#' Collapses drug-level susceptibility data to antibiotic-class level per
#' isolate (\eqn{R_{e,k,c} = 1} if resistant to \strong{any} drug in class
#' \eqn{c}), then computes marginal resistance for every pathogen x class
#' combination found in the data.
#'
#' Classes whose marginal resistance is at or below \code{zero_threshold} are
#' listed in \code{$near_zero} as a flag for downstream use -- they are
#' \strong{not} removed here.
#'
#' @param data Data frame. Pre-processed AMR data at isolate x antibiotic
#'   level (one row per isolate-antibiotic combination). Results must already
#'   be binary (\code{"S"} / \code{"R"}); no reclassification is applied.
#' @param pathogen_col Character. Column with pathogen names.
#'   Default \code{"organism_name"}.
#' @param org_group_col Character. Column with organism group labels.
#'   Default \code{"org_group"}.
#' @param isolate_col Character. Column uniquely identifying each isolate.
#'   Default \code{"isolate_id"}.
#' @param antibiotic_class_col Character. Column with the antibiotic class
#'   for each drug. Default \code{"antibiotic_class"}.
#' @param antibiotic_value_col Character. Column with susceptibility result
#'   (\code{"S"} or \code{"R"}). Default \code{"antibiotic_value"}.
#' @param zero_threshold Numeric. Classes with
#'   \code{marginal_resistance <= zero_threshold} are listed in
#'   \code{$near_zero}. Default \code{0}.
#' @param min_n_tested Integer or \code{NULL}. Minimum number of isolates that
#'   must have been tested for a pathogen-class combination to be retained.
#'   Combinations with \code{n_tested < min_n_tested} are dropped from
#'   \code{$marginal} \strong{and} from \code{$class_long}, so the exclusion
#'   propagates automatically into \code{compute_pairwise_coresistance()} and
#'   \code{compute_resistance_profiles()}. Set to \code{NULL} or \code{0} to
#'   disable the filter. Default \code{30}.
#' @param facility_col Character or \code{NULL}. Name of the column identifying
#'   the facility/site. When provided together with \code{facility_name}, data
#'   are filtered to the specified facility \strong{before} any computation.
#'   Both \code{facility_col} and \code{facility_name} must be supplied together
#'   or both left \code{NULL}. When provided, \code{facility_col} is also
#'   retained in \code{$class_long} and \code{$marginal} so that downstream
#'   steps can apply the same filter. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. The facility value to retain
#'   (matched via \code{==} against \code{facility_col}). Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Name of the column containing
#'   patient outcomes (e.g. \code{"final_outcome"}). When provided together with
#'   \code{outcome_value}, data are filtered to isolates with the specified
#'   outcome \strong{before} any computation. Both \code{outcome_col} and
#'   \code{outcome_value} must be supplied together or both left \code{NULL}.
#'   When provided, \code{outcome_col} is retained in \code{$class_long} and
#'   \code{$marginal} so that downstream steps can apply the same filter.
#'   Default \code{NULL}.
#' @param outcome_value Character or \code{NULL}. The outcome value to retain
#'   (e.g. \code{"discharged"}, \code{"dead"}; matched via \code{==} against
#'   \code{outcome_col}). Default \code{NULL}.
#'
#' @return Named list:
#' \describe{
#'   \item{\code{marginal}}{Data frame with columns: \code{pathogen_col},
#'     \code{org_group_col}, \code{antibiotic_class_col}, \code{n_tested},
#'     \code{n_resistant}, \code{marginal_resistance}. Sorted descending by
#'     \code{marginal_resistance} within each pathogen.}
#'   \item{\code{near_zero}}{Subset of \code{marginal} where
#'     \code{marginal_resistance <= zero_threshold}. These classes are
#'     candidates for exclusion in downstream profiling.}
#'   \item{\code{class_long}}{Collapsed isolate x pathogen x class data frame
#'     (columns: \code{isolate_col}, \code{pathogen_col}, \code{org_group_col},
#'     \code{antibiotic_class_col}, \code{class_result}). Pass this directly
#'     to \code{compute_pairwise_coresistance()}.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(
#'   data                 = amr_clean,
#'   pathogen_col         = "organism_name",
#'   org_group_col        = "org_group",
#'   isolate_col          = "isolate_id",
#'   antibiotic_class_col = "antibiotic_class",
#'   antibiotic_value_col = "antibiotic_value"
#' )
#'
#' marg$marginal # full marginal resistance table
#' marg$near_zero # classes flagged as near-zero
#' }
compute_marginal_resistance <- function(
  data,
  pathogen_col = "organism_name",
  org_group_col = "org_group",
  isolate_col = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  antibiotic_value_col = "antibiotic_value",
  zero_threshold = 0,
  min_n_tested = 30,
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL
) {
  .check_cols(data, c(
    pathogen_col, org_group_col, isolate_col,
    antibiotic_class_col, antibiotic_value_col
  ))

  # -- Optional facility filter ----------------------------------------------

  if (!is.null(facility_col) || !is.null(facility_name)) {
    if (is.null(facility_col) || is.null(facility_name)) {
      stop("Both facility_col and facility_name must be provided together, or both NULL.")
    }
    .check_cols(data, facility_col)
    n_before <- nrow(data)
    data <- data[data[[facility_col]] == facility_name, , drop = FALSE]
    if (nrow(data) == 0L) {
      stop(sprintf("No rows found for %s = '%s'.", facility_col, facility_name))
    }
    message(sprintf(
      "Facility filter applied: %s = '%s' (%d -> %d rows).",
      facility_col, facility_name, n_before, nrow(data)
    ))
  }

  # -- Optional outcome filter -----------------------------------------------

  if (!is.null(outcome_col) || !is.null(outcome_value)) {
    if (is.null(outcome_col) || is.null(outcome_value)) {
      stop("Both outcome_col and outcome_value must be provided together, or both NULL.")
    }
    .check_cols(data, outcome_col)
    n_before <- nrow(data)
    data <- data[data[[outcome_col]] == outcome_value, , drop = FALSE]
    if (nrow(data) == 0L) {
      stop(sprintf("No rows found for %s = '%s'.", outcome_col, outcome_value))
    }
    message(sprintf(
      "Outcome filter applied: %s = '%s' (%d -> %d rows).",
      outcome_col, outcome_value, n_before, nrow(data)
    ))
  }

  # -- Collapse to isolate x pathogen x class -------------------------------
  #
  # R_{e,k,c} = 1 if resistant to ANY drug in class c.

  iso_grp <- c(isolate_col, pathogen_col, org_group_col, antibiotic_class_col)
  if (!is.null(facility_col)) iso_grp <- c(facility_col, iso_grp)
  if (!is.null(outcome_col)) iso_grp <- c(outcome_col, iso_grp)

  class_long <- data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(iso_grp))) %>%
    dplyr::summarise(
      class_result = dplyr::if_else(
        any(!!rlang::sym(antibiotic_value_col) == "R"), "R", "S"
      ),
      .groups = "drop"
    )

  message(sprintf(
    "Collapsed to class level: %d isolates | %d pathogens | %d classes.",
    dplyr::n_distinct(class_long[[isolate_col]]),
    dplyr::n_distinct(class_long[[pathogen_col]]),
    dplyr::n_distinct(class_long[[antibiotic_class_col]])
  ))

  # -- Marginal resistance per pathogen x class -----------------------------

  marg_grp <- c(pathogen_col, org_group_col, antibiotic_class_col)
  if (!is.null(facility_col)) marg_grp <- c(facility_col, marg_grp)
  if (!is.null(outcome_col)) marg_grp <- c(outcome_col, marg_grp)

  marginal <- class_long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(marg_grp))) %>%
    dplyr::summarise(
      n_tested    = dplyr::n(),
      n_resistant = sum(class_result == "R", na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      marginal_resistance = n_resistant / n_tested
    ) %>%
    dplyr::arrange(
      !!rlang::sym(pathogen_col),
      dplyr::desc(marginal_resistance)
    )

  # -- Flag near-zero classes ------------------------------------------------

  near_zero <- marginal %>%
    dplyr::filter(marginal_resistance <= zero_threshold)

  if (nrow(near_zero) > 0) {
    message(sprintf(
      "%d pathogen-class combination(s) have marginal_resistance <= %g (listed in $near_zero):",
      nrow(near_zero), zero_threshold
    ))
    print(
      near_zero[, c(
        pathogen_col, antibiotic_class_col,
        "n_tested", "n_resistant", "marginal_resistance"
      )],
      row.names = FALSE
    )
  } else {
    message(sprintf(
      "No classes with marginal_resistance <= %g.", zero_threshold
    ))
  }

  # -- Minimum-tests threshold (optional) ------------------------------------
  # Drop pathogen-class combinations with too few isolates tested.
  # class_long is filtered to match so the exclusion carries through to
  # compute_pairwise_coresistance() and compute_resistance_profiles().

  if (!is.null(min_n_tested) && min_n_tested > 0) {
    n_before <- nrow(marginal)
    excluded <- marginal[marginal$n_tested < min_n_tested, , drop = FALSE]
    marginal <- marginal[marginal$n_tested >= min_n_tested, , drop = FALSE]

    if (nrow(excluded) > 0L) {
      message(sprintf(
        "min_n_tested = %d: %d of %d pathogen-class combination(s) excluded (n_tested < %d):",
        min_n_tested, nrow(excluded), n_before, min_n_tested
      ))
      print(
        excluded[, c(pathogen_col, antibiotic_class_col, "n_tested"),
          drop = FALSE
        ],
        row.names = FALSE
      )
    } else {
      message(sprintf(
        "min_n_tested = %d: all %d pathogen-class combination(s) passed the threshold.",
        min_n_tested, n_before
      ))
    }

    # Re-derive near_zero from the already-filtered marginal
    near_zero <- marginal[marginal$marginal_resistance <= zero_threshold, ,
      drop = FALSE
    ]

    # Align class_long: remove isolate rows for excluded pathogen-class combos.
    # Join keys: pathogen x class (facility/outcome are already constant if filtered).
    join_keys <- c(pathogen_col, antibiotic_class_col)
    if (!is.null(facility_col)) join_keys <- c(facility_col, join_keys)
    if (!is.null(outcome_col)) join_keys <- c(outcome_col, join_keys)
    class_long <- dplyr::semi_join(
      class_long,
      marginal[, join_keys, drop = FALSE],
      by = join_keys
    )

    message(sprintf(
      "class_long after min_n_tested filter: %d isolate-class row(s) retained.",
      nrow(class_long)
    ))
  }

  return(list(
    marginal   = marginal,
    near_zero  = near_zero,
    class_long = class_long
  ))
}


# -- Step 2 --------------------------------------------------------------------

#' Compute Pairwise Co-resistance Matrices per Pathogen
#'
#' For every pathogen and every pair of antibiotic classes \eqn{(c_i, c_j)}
#' that were both tested:
#' \deqn{
#'   T_{k,i,j} = \sum_e \mathbf{1}(c_i \text{ tested} \land c_j \text{ tested})
#' }
#' \deqn{
#'   R_{k,i,j} = \sum_e \mathbf{1}(R_{e,k,c_i}=1 \land R_{e,k,c_j}=1)
#' }
#' \deqn{
#'   \text{Prev}_{k,i,j} = R_{k,i,j} \;/\; T_{k,i,j}
#' }
#'
#' Matrices are computed for \strong{all} tested classes -- no filtering by
#' marginal resistance or GBD core list is applied here.
#'
#' @param marginal_output The list returned by
#'   \code{compute_marginal_resistance()}. The \code{$class_long} element is
#'   used as input.
#' @param pathogen_col Character. Must match the column name used in Step 1.
#'   Default \code{"organism_name"}.
#' @param org_group_col Character. Default \code{"org_group"}.
#' @param isolate_col Character. Default \code{"isolate_id"}.
#' @param antibiotic_class_col Character. Default \code{"antibiotic_class"}.
#' @param min_co_tested Integer. Pairwise cells with fewer than
#'   \code{min_co_tested} co-tested isolates are set to \code{NA} in the
#'   prevalence matrix. Default \code{10}.
#' @param facility_col Character or \code{NULL}. Column identifying the
#'   facility/site. When provided together with \code{facility_name}, filters
#'   \code{class_long} to the specified facility before building matrices.
#'   For this to work, \code{compute_marginal_resistance()} must have been
#'   called with the same \code{facility_col} argument so that the column is
#'   present in \code{$class_long}. Both must be supplied together or both
#'   left \code{NULL}. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. Facility value to retain.
#'   Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Column containing patient
#'   outcomes. When provided together with \code{outcome_value}, filters
#'   \code{class_long} to isolates with the specified outcome before building
#'   matrices. \code{compute_marginal_resistance()} must have been called with
#'   the same \code{outcome_col} so that the column is present in
#'   \code{$class_long}. Both must be supplied together or both left
#'   \code{NULL}. Default \code{NULL}.
#' @param outcome_value Character or \code{NULL}. Outcome value to retain
#'   (e.g. \code{"discharged"}, \code{"dead"}). Default \code{NULL}.
#'
#' @return Named list. One entry per pathogen (keyed by pathogen name), each
#'   containing:
#' \describe{
#'   \item{\code{prevalence}}{Square symmetric matrix of pairwise co-resistance
#'     rates. Diagonal and cells with \code{n < min_co_tested} are \code{NA}.}
#'   \item{\code{T_matrix}}{Integer matrix of co-tested isolate counts.}
#'   \item{\code{R_matrix}}{Integer matrix of co-resistant isolate counts.}
#'   \item{\code{classes}}{Character vector of class names (row/column order).}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(amr_clean, ...)
#' co_res <- compute_pairwise_coresistance(marg)
#'
#' # Prevalence matrix for K. pneumoniae
#' co_res[["Klebsiella pneumoniae"]]$prevalence
#'
#' # Co-tested counts
#' co_res[["Klebsiella pneumoniae"]]$T_matrix
#' }
compute_pairwise_coresistance <- function(
  marginal_output,
  pathogen_col = "organism_name",
  org_group_col = "org_group",
  isolate_col = "isolate_id",
  antibiotic_class_col = "antibiotic_class",
  min_co_tested = 10,
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL
) {
  if (!is.list(marginal_output) || !"class_long" %in% names(marginal_output)) {
    stop("marginal_output must be the list returned by compute_marginal_resistance().")
  }

  class_long <- marginal_output$class_long

  .check_cols(class_long, c(
    isolate_col, pathogen_col,
    org_group_col, antibiotic_class_col
  ))

  # -- Optional facility filter ----------------------------------------------

  if (!is.null(facility_col) || !is.null(facility_name)) {
    if (is.null(facility_col) || is.null(facility_name)) {
      stop("Both facility_col and facility_name must be provided together, or both NULL.")
    }
    if (!facility_col %in% names(class_long)) {
      stop(sprintf(
        paste0(
          "facility_col '%s' not found in class_long. ",
          "Re-run compute_marginal_resistance() with the same ",
          "facility_col and facility_name to filter upstream."
        ),
        facility_col
      ))
    }
    n_before <- nrow(class_long)
    class_long <- class_long[class_long[[facility_col]] == facility_name, , drop = FALSE]
    if (nrow(class_long) == 0L) {
      stop(sprintf("No rows in class_long for %s = '%s'.", facility_col, facility_name))
    }
    message(sprintf(
      "Facility filter applied: %s = '%s' (%d -> %d rows in class_long).",
      facility_col, facility_name, n_before, nrow(class_long)
    ))
  }

  # -- Optional outcome filter -----------------------------------------------

  if (!is.null(outcome_col) || !is.null(outcome_value)) {
    if (is.null(outcome_col) || is.null(outcome_value)) {
      stop("Both outcome_col and outcome_value must be provided together, or both NULL.")
    }
    if (!outcome_col %in% names(class_long)) {
      stop(sprintf(
        paste0(
          "outcome_col '%s' not found in class_long. ",
          "Re-run compute_marginal_resistance() with the same ",
          "outcome_col and outcome_value to filter upstream."
        ),
        outcome_col
      ))
    }
    n_before <- nrow(class_long)
    class_long <- class_long[class_long[[outcome_col]] == outcome_value, , drop = FALSE]
    if (nrow(class_long) == 0L) {
      stop(sprintf("No rows in class_long for %s = '%s'.", outcome_col, outcome_value))
    }
    message(sprintf(
      "Outcome filter applied: %s = '%s' (%d -> %d rows in class_long).",
      outcome_col, outcome_value, n_before, nrow(class_long)
    ))
  }

  pathogens <- sort(unique(class_long[[pathogen_col]]))
  out <- list()

  for (path in pathogens) {
    org_data <- class_long[class_long[[pathogen_col]] == path, ]
    classes <- sort(unique(org_data[[antibiotic_class_col]]))
    n_c <- length(classes)

    if (n_c < 2) {
      message(sprintf("'%s': fewer than 2 classes tested, skipping.", path))
      next
    }

    # Wide: rows = isolates, cols = classes, values = "R" / "S" / NA (not tested)
    iso_wide <- org_data %>%
      tidyr::pivot_wider(
        id_cols     = !!rlang::sym(isolate_col),
        names_from  = !!rlang::sym(antibiotic_class_col),
        values_from = class_result
        # isolates not tested for a class -> NA
      )

    # Binary matrix: 1 = R, 0 = S, NA = not tested
    bin_mat <- iso_wide[, classes, drop = FALSE] %>%
      dplyr::mutate(dplyr::across(
        dplyr::everything(),
        ~ dplyr::case_when(.x == "R" ~ 1L, .x == "S" ~ 0L, TRUE ~ NA_integer_)
      )) %>%
      as.matrix()

    # Pairwise T and R
    pairwise_T <- matrix(0L, n_c, n_c, dimnames = list(classes, classes))
    pairwise_R <- matrix(0L, n_c, n_c, dimnames = list(classes, classes))

    for (i in seq_len(n_c)) {
      for (j in i:n_c) {
        both_tested <- !is.na(bin_mat[, i]) & !is.na(bin_mat[, j])
        pairwise_T[i, j] <- pairwise_T[j, i] <- sum(both_tested)
        pairwise_R[i, j] <- pairwise_R[j, i] <- sum(
          both_tested & bin_mat[, i] == 1L & bin_mat[, j] == 1L
        )
      }
    }

    # Prevalence matrix: NA where co-tested < min_co_tested or on diagonal
    prev_mat <- ifelse(
      pairwise_T < min_co_tested,
      NA_real_,
      pairwise_R / pairwise_T
    )
    dimnames(prev_mat) <- list(classes, classes)
    diag(prev_mat) <- NA_real_

    out[[path]] <- list(
      prevalence = prev_mat,
      T_matrix   = pairwise_T,
      R_matrix   = pairwise_R,
      classes    = classes
    )

    message(sprintf(
      "'%s': %d classes, %d isolates -- co-resistance matrix built.",
      path, n_c, nrow(bin_mat)
    ))
  }

  return(out)
}


# -- Step 3 --------------------------------------------------------------------

#' Compute Resistance Profile Probabilities per Pathogen
#'
#' For each pathogen \eqn{k} with \eqn{n_k} antibiotic classes, enumerates all
#' \eqn{2^{n_k}} binary resistance profiles \eqn{\delta \in \{0,1\}^{n_k}} and
#' estimates their probabilities by solving a simplex-constrained weighted
#' least-squares Quadratic Programme (GBD equation 7.5.1.3):
#'
#' \deqn{
#'   \hat{p} = \arg\min_{p \in \Delta_{2^n}}
#'             \sum_{i=1}^{m} \frac{(m_i^\top p - v_i)^2}{\sigma_i^2}
#' }
#'
#' where \eqn{m = n(n+1)/2} data-derived linear constraints encode
#' \strong{n marginal} resistance rates and \strong{n(n-1)/2 pairwise}
#' co-resistance rates, and \eqn{\Delta} is the standard probability simplex.
#'
#' \subsection{Constraint rows in M}{
#'   \describe{
#'     \item{Marginal (rows 1 to n)}{
#'       Row \eqn{d}: \eqn{M_{d,\delta} = 1} iff \eqn{\delta_d = 1}
#'       (i.e., class \eqn{d} is resistant in profile \eqn{\delta}).
#'       Constraint: \eqn{\sum_\delta M_{d,\delta}\,p_\delta = \hat{r}_{kd}}.
#'     }
#'     \item{Pairwise (rows n+1 to m)}{
#'       Row \eqn{(d_1,d_2)}: \eqn{M_{d_1 d_2,\delta} = 1} iff
#'       \eqn{\delta_{d_1} = 1 \land \delta_{d_2} = 1}.
#'       Constraint: \eqn{\sum_\delta M_{d_1 d_2,\delta}\,p_\delta =
#'       \hat{r}_{k,d_1 d_2}}.
#'       When a pairwise estimate is unavailable (too few co-tested isolates),
#'       the product of marginals (independence assumption) is used as fallback.
#'     }
#'   }
#' }
#'
#' The QP is solved via \code{quadprog::solve.QP}. A small ridge term
#' (\code{ridge}) is added to the Hessian to guarantee strict
#' positive-definiteness. On solver failure the pathogen gets a uniform
#' distribution over all profiles.
#'
#' \subsection{Performance Notes}{
#'   This function has been optimized for speed with vectorized profile generation
#'   and label creation (10-100x faster than previous versions). However,
#'   computational complexity is still exponential in the number of classes:
#'   \itemize{
#'     \item \strong{n <= 14}: Fast (seconds to minutes)
#'     \item \strong{n = 15}: Moderate (minutes)
#'     \item \strong{n >= 16}: Slow and memory-intensive (use \code{top_n_classes})
#'   }
#'   For large datasets with many pathogens, consider using \code{top_n_classes}
#'   to limit each pathogen to its most-tested classes (e.g., \code{top_n_classes = 12}).
#' }
#'
#' @param marginal_output    List returned by \code{compute_marginal_resistance()}.
#' @param coresistance_output List returned by
#'   \code{compute_pairwise_coresistance()}.
#' @param pathogens          Character vector. Pathogen(s) to process.
#'   \code{NULL} (default) processes every pathogen in \code{marginal_output}.
#' @param top_n_pathogens    Integer or \code{NULL} (default). When set, only
#'   the top \code{top_n_pathogens} pathogens ranked by total isolates tested
#'   (sum of \code{n_tested} across all antibiotic classes, descending) are
#'   processed. Applied after the \code{pathogens} argument filter. Useful for
#'   focusing on the most data-rich pathogens -- e.g. \code{top_n_pathogens = 5}
#'   runs profiles for only the 5 most-tested pathogens.
#' @param exclude_near_zero  Logical. If \code{TRUE} (default), antibiotic
#'   classes that appear in \code{marginal_output$near_zero} for a given
#'   pathogen are excluded from profile enumeration.
#' @param top_n_classes      Integer or \code{NULL} (default). When set, only
#'   the top \code{top_n_classes} antibiotic classes ranked by \code{n_tested}
#'   (descending) are kept per pathogen before profile enumeration. Useful for
#'   capping the combinatorial explosion (2^n profiles) for pathogens tested
#'   against many drug classes -- e.g. \code{top_n_classes = 5} gives at most
#'   32 profiles. Applied after \code{exclude_near_zero}.
#' @param sigma_sq           Positive numeric. Assumed variance for each
#'   constraint (uniform). Default \code{1}.
#' @param ridge              Positive numeric. Ridge term added to the QP
#'   Hessian for numerical stability. Default \code{1e-8}.
#' @param pathogen_col       Character. Column name for pathogens. Must match
#'   the column used in Steps 1-2. Default \code{"organism_name"}.
#' @param antibiotic_class_col Character. Column name for antibiotic classes.
#'   Default \code{"antibiotic_class"}.
#' @param facility_col Character or \code{NULL}. Column identifying the
#'   facility/site. When provided together with \code{facility_name}, filters
#'   \code{marginal_output$marginal} and \code{marginal_output$near_zero} to
#'   the specified facility before profile enumeration. For this to work,
#'   \code{compute_marginal_resistance()} must have been called with the same
#'   \code{facility_col} argument. Both must be supplied together or both left
#'   \code{NULL}. Default \code{NULL}.
#' @param facility_name Character or \code{NULL}. Facility value to retain.
#'   Default \code{NULL}.
#' @param outcome_col Character or \code{NULL}. Column containing patient
#'   outcomes. When provided together with \code{outcome_value}, filters
#'   \code{marginal_output$marginal} and \code{marginal_output$near_zero} to
#'   the specified outcome before profile enumeration.
#'   \code{compute_marginal_resistance()} must have been called with the same
#'   \code{outcome_col} argument. Both must be supplied together or both left
#'   \code{NULL}. Default \code{NULL}.
#' @param outcome_value Character or \code{NULL}. Outcome value to retain
#'   (e.g. \code{"discharged"}, \code{"dead"}). Default \code{NULL}.
#' @param n_cores Integer. Number of CPU cores for parallel computation.
#'   Default \code{1L} (sequential).
#'
#' @return Named list, one entry per pathogen, each a list with:
#' \describe{
#'   \item{\code{profiles}}{Data frame: \code{profile} (character label,
#'     e.g.\ \code{"RSR"}), \code{probability} (\eqn{\hat{p}_\delta}), and
#'     one binary (0/1) integer column per antibiotic class indicating whether
#'     that class is resistant (\code{1}) or susceptible (\code{0}) in each
#'     profile.}
#'   \item{\code{classes}}{Character vector of antibiotic class names used
#'     (alphabetical; bit 0 = classes[1]).}
#'   \item{\code{n_classes}}{Integer. Number of classes used.}
#'   \item{\code{constraint_residuals}}{Named numeric vector of
#'     \eqn{m_i^\top \hat{p} - v_i} for each constraint. Small absolute
#'     values indicate good constraint satisfaction.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' marg <- compute_marginal_resistance(amr_clean)
#' co_res <- compute_pairwise_coresistance(marg)
#' rp <- compute_resistance_profiles(marg, co_res)
#'
#' # All profiles and probabilities for K. pneumoniae
#' rp[["Klebsiella pneumoniae"]]$profiles
#'
#' # Constraint residuals (quality check)
#' rp[["Klebsiella pneumoniae"]]$constraint_residuals
#'
#' # Single pathogen
#' rp_kp <- compute_resistance_profiles(
#'   marg, co_res,
#'   pathogens = "Klebsiella pneumoniae"
#' )
#' }
compute_resistance_profiles <- function(
  marginal_output,
  coresistance_output,
  pathogens = NULL,
  top_n_pathogens = NULL,
  exclude_near_zero = TRUE,
  top_n_classes = NULL,
  sigma_sq = 1,
  ridge = 1e-8,
  pathogen_col = "organism_name",
  antibiotic_class_col = "antibiotic_class",
  facility_col = NULL,
  facility_name = NULL,
  outcome_col = NULL,
  outcome_value = NULL,
  n_cores = 1L
) {
  # -- Input validation -------------------------------------------------------
  if (!is.list(marginal_output) ||
    !all(c("marginal", "near_zero", "class_long") %in% names(marginal_output))) {
    stop("marginal_output must be the list returned by compute_marginal_resistance().")
  }
  if (!is.list(coresistance_output)) {
    stop("coresistance_output must be the list returned by compute_pairwise_coresistance().")
  }
  has_osqp <- requireNamespace("osqp", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)
  has_quadprog <- requireNamespace("quadprog", quietly = TRUE)
  if (!has_osqp && !has_quadprog) {
    stop("Either 'osqp' + 'Matrix' (recommended) or 'quadprog' must be installed.")
  }

  if (!is.null(n_cores)) {
    n_cores <- as.integer(n_cores)
    if (is.na(n_cores) || n_cores < 1L) stop("n_cores must be a positive integer.")
  } else {
    n_cores <- 1L
  }

  if (!is.null(top_n_pathogens)) {
    if (!is.numeric(top_n_pathogens) || length(top_n_pathogens) != 1 ||
      top_n_pathogens < 1 || top_n_pathogens != round(top_n_pathogens)) {
      stop("top_n_pathogens must be a single positive integer.")
    }
    top_n_pathogens <- as.integer(top_n_pathogens)
  }

  if (!is.null(top_n_classes)) {
    if (!is.numeric(top_n_classes) || length(top_n_classes) != 1 ||
      top_n_classes < 1 || top_n_classes != round(top_n_classes)) {
      stop("top_n_classes must be a single positive integer.")
    }
    top_n_classes <- as.integer(top_n_classes)
  }

  # -- Optional facility filter ----------------------------------------------

  if (!is.null(facility_col) || !is.null(facility_name)) {
    if (is.null(facility_col) || is.null(facility_name)) {
      stop("Both facility_col and facility_name must be provided together, or both NULL.")
    }
    if (!facility_col %in% names(marginal_output$marginal)) {
      stop(sprintf(
        paste0(
          "facility_col '%s' not found in marginal_output$marginal. ",
          "Re-run compute_marginal_resistance() with the same ",
          "facility_col and facility_name to filter upstream."
        ),
        facility_col
      ))
    }
    marginal_output$marginal <- marginal_output$marginal[
      marginal_output$marginal[[facility_col]] == facility_name, ,
      drop = FALSE
    ]
    marginal_output$near_zero <- marginal_output$near_zero[
      marginal_output$near_zero[[facility_col]] == facility_name, ,
      drop = FALSE
    ]
    if (nrow(marginal_output$marginal) == 0L) {
      stop(sprintf("No rows in marginal for %s = '%s'.", facility_col, facility_name))
    }
    message(sprintf(
      "Facility filter applied: %s = '%s'.", facility_col, facility_name
    ))
  }

  # -- Optional outcome filter -----------------------------------------------

  if (!is.null(outcome_col) || !is.null(outcome_value)) {
    if (is.null(outcome_col) || is.null(outcome_value)) {
      stop("Both outcome_col and outcome_value must be provided together, or both NULL.")
    }
    if (!outcome_col %in% names(marginal_output$marginal)) {
      stop(sprintf(
        paste0(
          "outcome_col '%s' not found in marginal_output$marginal. ",
          "Re-run compute_marginal_resistance() with the same ",
          "outcome_col and outcome_value to filter upstream."
        ),
        outcome_col
      ))
    }
    marginal_output$marginal <- marginal_output$marginal[
      marginal_output$marginal[[outcome_col]] == outcome_value, ,
      drop = FALSE
    ]
    marginal_output$near_zero <- marginal_output$near_zero[
      marginal_output$near_zero[[outcome_col]] == outcome_value, ,
      drop = FALSE
    ]
    if (nrow(marginal_output$marginal) == 0L) {
      stop(sprintf("No rows in marginal for %s = '%s'.", outcome_col, outcome_value))
    }
    message(sprintf(
      "Outcome filter applied: %s = '%s'.", outcome_col, outcome_value
    ))
  }

  all_pathogens <- sort(unique(marginal_output$marginal[[pathogen_col]]))

  if (!is.null(pathogens)) {
    missing_p <- setdiff(pathogens, all_pathogens)
    if (length(missing_p) > 0) {
      stop(sprintf(
        "Pathogen(s) not found in marginal_output: %s",
        paste(missing_p, collapse = ", ")
      ))
    }
    all_pathogens <- pathogens
  }

  # -- Restrict to top N most-tested pathogens ---------------------------------
  if (!is.null(top_n_pathogens)) {
    path_tested <- marginal_output$marginal %>%
      dplyr::filter(!!rlang::sym(pathogen_col) %in% all_pathogens) %>%
      dplyr::group_by(!!rlang::sym(pathogen_col)) %>%
      dplyr::summarise(total_tested = sum(n_tested), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(total_tested))

    top_paths <- path_tested[[pathogen_col]][
      seq_len(min(top_n_pathogens, nrow(path_tested)))
    ]
    excluded <- setdiff(all_pathogens, top_paths)

    message(sprintf(
      "top_n_pathogens = %d: keeping %d pathogen(s), excluding %d.",
      top_n_pathogens, length(top_paths), length(excluded)
    ))
    if (length(excluded) > 0) {
      message("  Excluded: ", paste(excluded, collapse = ", "))
    }
    message("  Selected (ranked by total isolates tested):")
    print(path_tested[path_tested[[pathogen_col]] %in% top_paths, ],
      row.names = FALSE
    )

    all_pathogens <- top_paths # preserves descending rank order
  }

  out <- list()
  skipped_single <- list()

  n_pathogens_total <- length(all_pathogens)
  if (n_pathogens_total > 1) {
    message(sprintf(
      "\nProcessing %d pathogen(s) with %d core(s)...",
      n_pathogens_total, n_cores
    ))
  }

  # -- Per-pathogen worker (closure over outer env) --------------------------
  .process_one <- function(path_idx) {
    path <- all_pathogens[path_idx]

    if (n_pathogens_total > 1) {
      message(sprintf(
        "\n[%d/%d] Processing '%s'...",
        path_idx, n_pathogens_total, path
      ))
    }

    # -- Determine antibiotic classes ---------------------------------------
    marg_k <- marginal_output$marginal[
      marginal_output$marginal[[pathogen_col]] == path,
    ]

    if (exclude_near_zero) {
      nz_classes <- marginal_output$near_zero[
        marginal_output$near_zero[[pathogen_col]] == path,
        antibiotic_class_col,
        drop = TRUE
      ]
      marg_k <- marg_k[!marg_k[[antibiotic_class_col]] %in% nz_classes, ]
    }

    if (!is.null(top_n_classes) && nrow(marg_k) > top_n_classes) {
      marg_k <- marg_k[order(marg_k$n_tested, decreasing = TRUE), ][
        seq_len(top_n_classes),
      ]
      message(sprintf(
        "'%s': restricting to top %d classes by n_tested.", path, top_n_classes
      ))
    }

    classes <- sort(marg_k[[antibiotic_class_col]])
    n <- length(classes)

    if (n == 0) {
      message(sprintf("'%s': no classes remain after filtering, skipping.", path))
      return(list(type = "empty", path = path))
    }

    if (n < 2) {
      return(list(
        type = "skipped",
        path = path,
        data = data.frame(
          pathogen            = path,
          antibiotic_class    = marg_k[[antibiotic_class_col]],
          n_tested            = marg_k$n_tested,
          marginal_resistance = marg_k$marginal_resistance,
          stringsAsFactors    = FALSE
        )
      ))
    }

    n_profiles <- 2L^n

    if (n > 18L) {
      warning(sprintf(
        "'%s': %d classes -> 2^%d = %d profiles. This may be very slow and memory-intensive.",
        path, n, n, n_profiles
      ))
    }

    # -- Marginal resistance vector r_kd ------------------------------------
    r_marg <- setNames(
      marg_k$marginal_resistance,
      marg_k[[antibiotic_class_col]]
    )[classes]

    # -- Pairwise co-resistance matrix (may be NULL) ------------------------
    co_mat <- if (path %in% names(coresistance_output)) {
      coresistance_output[[path]]$prevalence
    } else {
      NULL
    }

    # -- Enumerate 2^n resistance profiles ---------------------------------
    profiles_mat <- matrix(
      as.integer(
        outer(0L:(n_profiles - 1L), 2L^(0L:(n - 1L)), bitwAnd) > 0L
      ),
      nrow = n_profiles, ncol = n,
      dimnames = list(NULL, classes)
    )

    char_mat <- matrix("S", nrow = n_profiles, ncol = n)
    char_mat[profiles_mat == 1L] <- "R"
    profile_labels <- do.call(paste0, as.data.frame(char_mat))

    # -- Constraint matrix M (m x 2^n) and target vector v -----------------
    M_marg <- t(profiles_mat)
    v_marg <- r_marg

    pairs_mat <- utils::combn(n, 2L)
    n_pair <- ncol(pairs_mat)
    d1_idx <- pairs_mat[1L, ]
    d2_idx <- pairs_mat[2L, ]
    c1_names <- classes[d1_idx]
    c2_names <- classes[d2_idx]
    pair_names <- paste0(c1_names, "_", c2_names)

    M_pair <- t(
      profiles_mat[, d1_idx, drop = FALSE] *
        profiles_mat[, d2_idx, drop = FALSE]
    )

    r1 <- r_marg[c1_names]
    r2 <- r_marg[c2_names]

    if (!is.null(co_mat) &&
      !is.null(rownames(co_mat)) && !is.null(colnames(co_mat))) {
      in_rows <- c1_names %in% rownames(co_mat)
      in_cols <- c2_names %in% colnames(co_mat)
      can_look <- in_rows & in_cols
      co_vals <- rep(NA_real_, n_pair)
      if (any(can_look)) {
        co_vals[can_look] <- co_mat[cbind(
          c1_names[can_look],
          c2_names[can_look]
        )]
      }
    } else {
      co_vals <- rep(NA_real_, n_pair)
    }

    cap_vals <- pmin(r1, r2)
    has_co <- !is.na(co_vals)
    capped_co <- pmin(co_vals, cap_vals)
    was_capped <- has_co & !is.na(capped_co) & (capped_co < co_vals)
    v_pair <- ifelse(has_co, capped_co, r1 * r2)

    if (any(was_capped)) {
      message(sprintf(
        "'%s': %d pairwise value(s) capped to min(marginal): %s",
        path, sum(was_capped),
        paste(
          sprintf(
            "(%s,%s) %.4f->%.4f",
            c1_names[was_capped], c2_names[was_capped],
            co_vals[was_capped], capped_co[was_capped]
          ),
          collapse = "; "
        )
      ))
    }
    if (any(!has_co)) {
      message(sprintf(
        "'%s': %d pair(s) used independence fallback: %s",
        path, sum(!has_co),
        paste(sprintf("(%s,%s)", c1_names[!has_co], c2_names[!has_co]),
          collapse = "; "
        )
      ))
    }

    M <- rbind(M_marg, M_pair)
    v <- c(v_marg, v_pair)
    storage.mode(M) <- "double"

    # -- QP: simplex-constrained weighted least-squares ---------------------
    #
    # Minimise  (1/2) p^T H p  -  d^T p
    #   H = (2/sigma_sq) * M^T M  +  ridge * I   (guaranteed pos-def)
    #   d = (2/sigma_sq) * M^T v
    # Subject to: sum(p) = 1, p >= 0
    #
    # OSQP path: constraint matrix A = rbind([1...1], I_{2^n}) stored sparse
    #   -> O(2^n) memory vs O(4^n) for the dense diag() used by quadprog.
    # quadprog path: retained as fallback when osqp/Matrix are unavailable.
    coef <- 2.0 / sigma_sq
    H_mat <- coef * crossprod(M)
    diag(H_mat) <- diag(H_mat) + ridge
    d_qp <- coef * drop(crossprod(M, v))

    p_hat <- tryCatch(
      {
        if (has_osqp) {
          A_sp <- rbind(
            Matrix::Matrix(1.0, nrow = 1L, ncol = n_profiles, sparse = TRUE),
            Matrix::Diagonal(n_profiles)
          )
          prob <- osqp::osqp(
            P = Matrix::forceSymmetric(Matrix::Matrix(H_mat)),
            q = -d_qp,
            A = A_sp,
            l = c(1.0, rep(0.0, n_profiles)),
            u = c(1.0, rep(Inf, n_profiles)),
            pars = osqp::osqpSettings(
              verbose  = FALSE,
              eps_abs  = 1e-8,
              eps_rel  = 1e-8,
              max_iter = 10000L,
              polish   = TRUE
            )
          )
          res <- prob$solve()
          if (!(res$info$status %in% c("solved", "solved_inaccurate"))) {
            stop(paste("OSQP status:", res$info$status))
          }
          pmax(res$x, 0.0)
        } else {
          # quadprog fallback -- dense Amat, slow for n > 12
          Amat <- cbind(rep(1.0, n_profiles), diag(n_profiles))
          bvec <- c(1.0, rep(0.0, n_profiles))
          sol <- quadprog::solve.QP(
            Dmat = H_mat, dvec = d_qp,
            Amat = Amat,  bvec = bvec, meq = 1L
          )
          pmax(sol$solution, 0.0)
        }
      },
      error = function(e) {
        warning(sprintf(
          "'%s': QP solver failed (%s). Returning uniform distribution.",
          path, conditionMessage(e)
        ))
        rep(1.0 / n_profiles, n_profiles)
      }
    )

    p_hat <- p_hat / sum(p_hat)

    # -- Constraint residuals -----------------------------------------------
    residuals <- drop(M %*% p_hat) - v
    names(residuals) <- c(
      paste0("marg_", classes),
      paste0("pair_", pair_names)
    )

    profiles_df <- data.frame(
      profile = profile_labels,
      probability = p_hat,
      stringsAsFactors = FALSE
    )
    profiles_df <- cbind(profiles_df, as.data.frame(profiles_mat))

    message(sprintf(
      "'%s': n=%d classes -> %d profiles. Max |residual| = %.5f.",
      path, n, n_profiles, max(abs(residuals))
    ))

    list(
      type = "success",
      path = path,
      result = list(
        profiles             = profiles_df,
        classes              = classes,
        n_classes            = n,
        constraint_residuals = residuals
      )
    )
  }

  # -- Execute: parallel or sequential ---------------------------------------
  if (n_cores > 1L) {
    all_results <- parallel::mclapply(
      seq_along(all_pathogens), .process_one,
      mc.cores = n_cores,
      mc.preschedule = FALSE # better load balancing across heterogeneous pathogens
    )
  } else {
    all_results <- lapply(seq_along(all_pathogens), .process_one)
  }

  # -- Collect results --------------------------------------------------------
  for (res in all_results) {
    if (is.null(res) || inherits(res, "try-error")) next
    if (res$type == "success") {
      out[[res$path]] <- res$result
    } else if (res$type == "skipped") {
      skipped_single[[res$path]] <- res$data
    }
  }

  # -- Summary: pathogens skipped because only 1 class remained --------------
  if (length(skipped_single) > 0) {
    skipped_tbl <- do.call(rbind, skipped_single)
    message(sprintf(
      "\n%d pathogen(s) skipped -- only 1 antibiotic class remained after filtering.",
      length(skipped_single)
    ))
    message(paste0(
      "For these pathogens the resistance profile IS the marginal resistance.\n",
      "  Tip: use top_n_classes, reduce zero_threshold, or set exclude_near_zero = FALSE",
      " to include them.\n"
    ))
    print(skipped_tbl, row.names = FALSE)
  }

  return(out)
}
