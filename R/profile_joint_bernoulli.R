# profile_joint_bernoulli.R
# Skeleton for the joint multivariate-Bernoulli resistance-profile model
# (generalized from the J2-L production model in
#  anumaan-analysis/Scripts/stewardship/resistance_profiles/).
#
# The equation being fit, for an event with covariates x and D chosen
# antibiotic classes, is a log-linear model over the 2^D joint resistance
# states s (each class R/S):
#
#   score(s) = sum_d B[s,d] * (x' beta_d)              <- D main effects, covariate-dependent
#            + sum_{k=2}^{D} sum_{c in combn(D,k)} V_k[s,c] * eta_{k,c}   <- globally shared interactions
#
#   P(state = s | tested classes) = exp(score(s)) / sum_{s' compatible} exp(score(s'))
#
# Untested classes are handled by exact marginalization: the denominator and
# numerator both sum only over states compatible with whatever was actually
# observed (this is what makes the model correct for partial AST panels,
# not an imputation).
#
# Number of parameters is driven entirely by D, the number of antibiotic
# classes chosen for the panel:
#   - D main-effect vectors beta_d (each length K = ncol(covariate design))
#   - choose(D,k) global interaction scalars for every order k = 2..D
#     (6 pairwise + 4 three-way + 1 four-way for D=4, as in J2-L)
# Covariates (the design matrix X) are entirely up to the caller -- this
# skeleton does not assume age/sex/hospital/ICU, it only needs *a* design
# matrix.
#
# ONE joint model is always fit across every configured stratum at once --
# never a silent separate fit per stratum. The interaction terms eta above
# have no stratum index at all, so they are automatically shared across
# whatever fixed-effect stratifier (hospital, state, district, country, or
# none) a caller's X encodes -- that choice is a config decision made in
# anumaan-analysis, never hardcoded to "hospital" here. See
# profile_build_fixed_effect_design() for the correctly-specified way to
# build that X (guards against the unconstrained-dummies-plus-intercept
# aliasing mistake).
#
# Functions:
#   - profile_enumerate_states    (2^D x D sign/state enumeration)
#   - profile_interaction_design  (state-space design matrix for one interaction order)
#   - profile_bernoulli_terms     (full parameter manifest: which columns are which order)
#   - profile_compatible_states   (marginalization bookkeeping for one partially-observed event)
#   - profile_build_fixed_effect_design (config-chosen fixed-effect stratifier -> X)
#   - profile_build_design        (the skeleton: Y + X -> Stan-ready data list)


#' Enumerate the 2^D Joint Resistance States
#'
#' Builds every combination of resistant (+1) / susceptible (-1) across
#' \code{D} antibiotic classes. Row order is the canonical state ordering
#' used by \code{\link{profile_interaction_design}} and
#' \code{\link{profile_compatible_states}} -- do not reorder after
#' generating it for a given model, the Stan data has to agree with this
#' ordering.
#'
#' @param D Integer. Number of antibiotic classes in the panel.
#' @return Integer matrix, \code{2^D} rows x \code{D} columns, entries in
#'   \code{c(-1, 1)} (susceptible / resistant).
#' @export
profile_enumerate_states <- function(D) {
  stopifnot(is.numeric(D), D >= 1, D == round(D))
  grid <- expand.grid(rep(list(c(-1L, 1L)), D))
  as.matrix(grid[, rev(seq_len(D)), drop = FALSE])
}


#' Interaction-Order Design Matrix Over the State Space
#'
#' For interaction order \code{k}, returns a \code{2^D x choose(D, k)} matrix
#' whose column for class-combination \code{c} is the elementwise product of
#' those \code{k} classes' sign columns in \code{states} -- the standard
#' saturated log-linear / Ising-model design. \code{k = 1} just returns
#' \code{states} itself (the main-effect design, \code{B} in the Stan model).
#'
#' @param states Output of \code{\link{profile_enumerate_states}}.
#' @param k Integer. Interaction order, \code{1 <= k <= ncol(states)}.
#' @return Numeric matrix, \code{nrow(states) x choose(ncol(states), k)},
#'   with \code{attr(, "combn")} recording which classes (by column index of
#'   \code{states}) each output column corresponds to.
#' @export
profile_interaction_design <- function(states, k) {
  D <- ncol(states)
  stopifnot(k >= 1, k <= D)
  if (k == 1) {
    out <- states
    attr(out, "combn") <- matrix(seq_len(D), nrow = 1)
    return(out)
  }
  combos <- utils::combn(D, k)
  out <- apply(combos, 2, function(cols) apply(states[, cols, drop = FALSE], 1, prod))
  out <- matrix(out, nrow = nrow(states))
  attr(out, "combn") <- combos
  out
}


#' Full Parameter Manifest for a D-Class Joint Bernoulli Model
#'
#' Enumerates every term the model will fit: the \code{D} covariate-dependent
#' main effects, plus one global scalar per class-combination for every
#' interaction order from 2 up to \code{max_order} (default: fully
#' saturated, i.e. up to \code{D}). This is the thing that makes the model
#' "automatic" in the number of antibiotic classes chosen -- adding or
#' removing a class changes \code{D}, which changes every count here without
#' the caller hand-editing interaction terms.
#'
#' @param D Integer. Number of antibiotic classes.
#' @param class_names Character vector of length \code{D}, optional. Used to
#'   label interaction terms (e.g. \code{"Aminoglycosides:Carbapenems"}).
#' @param max_order Integer, \code{2 <= max_order <= D}. Cap the highest
#'   interaction order fit (e.g. \code{max_order = 2} for pairwise-only --
#'   useful when \code{D} is large enough that the full \code{2^D} state
#'   space or \code{choose(D, D)} terms become impractical). Default: fully
#'   saturated (\code{= D}).
#' @return A data frame, one row per interaction term, columns
#'   \code{order} (1 = main effect), \code{term} (label), and
#'   \code{classes} (list-column of class-name vectors). Row count for
#'   orders 2..max_order is \code{sum_{k=2}^{max_order} choose(D, k)}.
#' @export
profile_bernoulli_terms <- function(D, class_names = NULL, max_order = D) {
  stopifnot(max_order >= 1, max_order <= D)
  if (is.null(class_names)) class_names <- paste0("class", seq_len(D))
  stopifnot(length(class_names) == D)

  rows <- list()
  for (k in seq_len(max_order)) {
    combos <- utils::combn(D, k)
    for (j in seq_len(ncol(combos))) {
      cls <- class_names[combos[, j]]
      rows[[length(rows) + 1]] <- data.frame(
        order = k,
        term = paste(cls, collapse = ":"),
        classes = I(list(cls))
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Which Joint States Are Compatible With a Partially-Observed Event
#'
#' For one event with possibly-missing classes, returns the indices (into
#' \code{states}, i.e. into \code{1:2^D}) of every joint state consistent
#' with what was actually tested. This is the exact-marginalization
#' machinery: untested classes are summed out by including every state that
#' agrees with the tested classes' signs, regardless of what the untested
#' classes would have been.
#'
#' @param states Output of \code{\link{profile_enumerate_states}}.
#' @param observed_sign Integer vector length \code{D}, entries in
#'   \code{c(-1, 0, 1)}: \code{-1}/\code{1} for susceptible/resistant,
#'   \code{0} (or \code{NA}) for untested.
#' @return Integer vector of row indices into \code{states}.
#' @export
profile_compatible_states <- function(states, observed_sign) {
  observed_sign[is.na(observed_sign)] <- 0L
  D <- ncol(states)
  stopifnot(length(observed_sign) == D)
  tested <- which(observed_sign != 0L)
  if (length(tested) == 0L) return(seq_len(nrow(states)))
  match_rows <- apply(states[, tested, drop = FALSE], 1, function(row) {
    all(row == observed_sign[tested])
  })
  which(match_rows)
}


#' Build a Covariate Design Matrix With a Config-Chosen Fixed-Effect Stratifier
#'
#' Builds the \code{X} matrix \code{\link{profile_build_design}} expects,
#' for the ONE joint model this family fits across every configured
#' stratum at once -- never a silent separate fit per stratum. The
#' interaction coefficients \code{\link{profile_build_design}} builds
#' (\code{eta}) have no stratum index at all (they are defined purely over
#' the state space), so they are automatically SHARED across every level
#' of whatever stratifier this function encodes here; nothing downstream
#' needs to know what that stratifier was.
#'
#' The stratifier itself is never assumed to be "hospital": it is whatever
#' column the caller's analysis config names via \code{group_col} --
#' hospital, state, district, country, or \code{NULL} for no fixed effects
#' at all (a fully pooled single-intercept model). That choice belongs to
#' whoever prepares the config in \code{anumaan-analysis}, exactly like
#' \code{group_col} is already generic everywhere else in this family
#' (\code{\link{profile_audit_support}}, \code{\link{profile_prepare_model_data}},
#' ...) -- this function does not special-case hospitals, and does not
#' decide the granularity for the caller.
#'
#' Guards against the failure mode an unconstrained design invites: a full
#' set of stratum dummy columns PLUS an intercept is rank-deficient
#' (aliased) -- \code{contrast = "reference"} drops one declared/default
#' level instead, and \code{contrast = "sum_to_zero"} uses deviation-from-
#' grand-mean coding with no dropped reference level at all. Either way,
#' the returned design's rank is checked and a rank deficiency (from this
#' or from collinear \code{covariates}) stops with a specific diagnosis
#' rather than silently fitting an unidentified model.
#'
#' @param ids Data frame, one row per event, matching
#'   \code{\link{profile_prepare_model_data}}'s \code{$ids} row order (its
#'   own \code{$ids} is the intended input, but any data frame with a
#'   stratum column works).
#' @param group_col Character or \code{NULL}. Column in \code{ids} naming
#'   the fixed-effect stratifier -- hospital, state, district, country,
#'   whatever the config specifies. Default \code{"group"}, matching
#'   \code{profile_prepare_model_data()$ids}'s own column name.
#'   \code{NULL} means NO fixed effects (fully pooled) -- an explicit,
#'   valid choice, not something to avoid by default.
#' @param covariates Data frame or \code{NULL}, same row order as
#'   \code{ids} (\code{profile_prepare_model_data()$covariates} is the
#'   intended input). Additional covariates (age, sex, ...) included
#'   alongside the stratifier. Default \code{NULL} (stratifier + intercept
#'   only).
#' @param contrast Character, \code{"reference"} (default) or
#'   \code{"sum_to_zero"}. \code{"reference"} drops one level
#'   (\code{reference_level}) and keeps an intercept -- standard treatment
#'   coding. \code{"sum_to_zero"} uses deviation-from-grand-mean coding
#'   instead, with no dropped reference level.
#' @param reference_level Character or \code{NULL}. Which level of
#'   \code{group_col} is dropped, when \code{contrast = "reference"}.
#'   \code{NULL} (default) uses the first level alphabetically -- R's own
#'   default, made explicit and recorded here rather than left implicit.
#'
#' @return A list: \code{X} (the design matrix), \code{term_names}
#'   (\code{colnames(X)}, in order), \code{group_col} (echoed, \code{NULL}
#'   if no stratifier), \code{reference_level} (the dropped level,
#'   \code{NA} for \code{sum_to_zero} or no stratifier), \code{contrast}
#'   (echoed), \code{levels} (every level of \code{group_col} present, in
#'   the order used), \code{rank} (\code{qr(X)$rank} -- always
#'   \code{== ncol(X)}, or the function has already stopped).
#' @export
profile_build_fixed_effect_design <- function(ids,
                                               group_col = "group",
                                               covariates = NULL,
                                               contrast = c("reference", "sum_to_zero"),
                                               reference_level = NULL) {
  contrast <- match.arg(contrast)
  N <- nrow(ids)
  if (!is.null(covariates) && nrow(covariates) != N) {
    stop("covariates must have the same number of rows as ids (one row per event).")
  }

  has_stratifier <- !is.null(group_col)
  if (has_stratifier && !group_col %in% names(ids)) {
    stop(sprintf("group_col '%s' not found in ids", group_col))
  }

  levels_used <- NULL
  reference_used <- NA_character_
  df <- data.frame(row.names = seq_len(N))

  if (has_stratifier) {
    if (anyNA(ids[[group_col]])) {
      stop(sprintf(
        "group_col '%s' has %d NA value(s) -- every event needs a known stratum; resolve or exclude those events explicitly before calling this function.",
        group_col, sum(is.na(ids[[group_col]]))
      ), call. = FALSE)
    }
    g <- factor(as.character(ids[[group_col]]))
    levels_used <- levels(g)
    if (length(levels_used) < 2L) {
      stop(sprintf(
        "group_col '%s' has only %d distinct level(s) in this data -- fixed effects need >= 2 levels to be identifiable; pass group_col = NULL for a fully pooled model instead.",
        group_col, length(levels_used)
      ))
    }
    if (contrast == "reference") {
      reference_used <- reference_level %||% levels_used[1]
      if (!reference_used %in% levels_used) {
        stop(sprintf(
          "reference_level '%s' is not a level of group_col '%s' present in this data: %s",
          reference_used, group_col, paste(levels_used, collapse = ", ")
        ))
      }
      g <- stats::relevel(g, ref = reference_used)
      cm <- stats::contr.treatment(length(levels(g)))
      rownames(cm) <- levels(g)
      colnames(cm) <- levels(g)[-1]
      stats::contrasts(g) <- cm
    } else {
      cm <- stats::contr.sum(length(levels_used))
      rownames(cm) <- levels_used
      colnames(cm) <- paste0("vs_mean_", levels_used[-length(levels_used)])
      stats::contrasts(g) <- cm
    }
    df$.fixed_effect_stratum <- g
  }

  if (!is.null(covariates)) {
    na_counts <- colSums(is.na(covariates))
    if (any(na_counts > 0L)) {
      stop(sprintf(
        "covariates contain NA value(s): %s. model.matrix()'s default na.action silently DROPS those rows, which would misalign X against Y (Y is not filtered the same way) -- impute or complete-case-filter the covariates (and Y/ids together) explicitly before calling this function, rather than letting a row-count mismatch surface later as a cryptic profile_build_design() error.",
        paste(sprintf("%s (%d NA)", names(na_counts)[na_counts > 0L], na_counts[na_counts > 0L]), collapse = ", ")
      ), call. = FALSE)
    }
    df <- cbind(df, covariates)
  }

  rhs <- if (ncol(df) == 0L) stats::as.formula("~1") else stats::as.formula("~.")
  X <- stats::model.matrix(rhs, data = df)
  colnames(X) <- sub("^\\.fixed_effect_stratum", "fixed_effect:", colnames(X))

  if (nrow(X) != N) {
    stop(sprintf(
      "Internal: model.matrix() returned %d rows for %d input events -- some row was silently dropped (most likely NA in a covariate or in group_col that slipped past the checks above). Refusing to return a misaligned X.",
      nrow(X), N
    ), call. = FALSE)
  }

  rank <- qr(X)$rank
  if (rank < ncol(X)) {
    stop(sprintf(
      "Design matrix is rank-deficient: %d columns, rank %d. Check for an unconstrained full set of stratum dummies plus an intercept (contrast = \"sum_to_zero\" avoids this), or redundant/collinear covariates.",
      ncol(X), rank
    ))
  }

  list(
    X = X, term_names = colnames(X),
    group_col = if (has_stratifier) group_col else NULL,
    reference_level = if (contrast == "reference") reference_used else NA_character_,
    contrast = contrast, levels = levels_used, rank = rank
  )
}


#' Build Stan-Ready Design for a D-Class Joint Bernoulli Resistance Model
#'
#' The general skeleton: takes a wide S/R outcome matrix and a caller-supplied
#' covariate design, and returns everything a Stan program needs to fit the
#' saturated (or order-capped) log-linear joint Bernoulli model with exact
#' marginalization of untested classes. \code{D} -- and therefore every
#' interaction count -- comes entirely from \code{ncol(Y)}; this function
#' does not hardcode a panel size or a covariate set. The interaction
#' coefficients (\code{eta}) have NO stratum index, so they are shared
#' across every level of whatever fixed-effect stratifier is encoded into
#' \code{X} (hospital, state, district, country, or none) -- this is the
#' SAME fit regardless of how many strata \code{X} represents, never a
#' separate fit per stratum (see \code{\link{profile_build_fixed_effect_design}},
#' the recommended way to build \code{X} with that stratifier correctly
#' specified).
#'
#' @param Y Numeric/integer matrix or data frame, \code{N} events x \code{D}
#'   classes, entries in \code{c(0, 1, NA)} (susceptible / resistant /
#'   untested). Column names become the class labels unless
#'   \code{class_names} is supplied.
#' @param X Numeric matrix, \code{N} x \code{K}: the covariate design for
#'   the main effects. Entirely up to the caller -- build it with
#'   \code{model.matrix()} from whatever formula is appropriate (age, sex,
#'   hospital, ICU, whatever). Must have \code{nrow(X) == nrow(Y)}.
#' @param class_names Character vector length \code{D}, optional (default:
#'   \code{colnames(Y)}).
#' @param max_order Integer, passed to \code{\link{profile_bernoulli_terms}}.
#'   Default: fully saturated (\code{= D}).
#' @param beta_sd Numeric scalar or length-\code{K} vector: prior SD for the
#'   main-effect coefficients. Default \code{2.5} for an intercept-like
#'   first column, \code{1.5} for the rest, matching J2-L's convention, but
#'   fully overridable.
#' @param interaction_sd Numeric scalar, or a function \code{order -> sd},
#'   giving the prior SD for each interaction order's global scalars.
#'   Default halves the SD per additional order
#'   (\code{function(k) 1.5 / (k - 1)}), matching J2-L's 1.5 / 0.75(*) / 0.5 / 0.35
#'   decreasing schedule in spirit ((*) J2-L's order-2 term is actually a
#'   separate ICU/location effect, not part of this skeleton's covariate
#'   block -- fold it into \code{X} instead).
#' @param grainsize Integer. Stan \code{reduce_sum} grainsize. Default 64L.
#'
#' @return A named list, directly usable as \code{cmdstanr}/\code{rstan}
#'   Stan data, plus a \code{manifest} data frame (one row per parameter,
#'   from \code{\link{profile_bernoulli_terms}}) and \code{design} (the
#'   state enumeration + per-order design matrices, for inspection/plotting).
#'   Stan list elements: \code{N_events, D, K, S} (\code{= 2^D}), \code{P}
#'   (total interaction-term count), \code{X}, \code{B} (\code{S x D} main-effect
#'   design), \code{V} (\code{S x P} concatenated higher-order design),
#'   \code{n_compatible}, \code{compatible_state} (\code{N_events x S}, only
#'   the first \code{n_compatible[e]} entries of row \code{e} are used),
#'   \code{beta_sd} (length \code{K}), \code{eta_sd} (length \code{P}).
#'
#' @seealso The reference Stan model this pairs with:
#'   \code{system.file("stan", "joint_bernoulli_general.stan", package = "anumaan")}
#' @export
profile_build_design <- function(Y, X,
                                 class_names = NULL,
                                 max_order = NULL,
                                 beta_sd = NULL,
                                 interaction_sd = function(k) 1.5 / (k - 1),
                                 grainsize = 64L) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  N <- nrow(Y)
  D <- ncol(Y)
  K <- ncol(X)
  if (nrow(X) != N) stop("nrow(X) must equal nrow(Y) (one covariate row per event).")
  if (is.null(class_names)) {
    class_names <- colnames(Y)
    if (is.null(class_names)) class_names <- paste0("class", seq_len(D))
  }
  if (is.null(max_order)) max_order <- D
  if (D < 2) stop("Need at least 2 antibiotic classes for a joint model (D = ", D, ").")

  # recode Y to {-1, 0 (untested), 1} sign coding used by the state design
  Ysign <- matrix(0L, N, D)
  Ysign[Y == 1] <- 1L
  Ysign[Y == 0] <- -1L
  Ysign[is.na(Y)] <- 0L

  states <- profile_enumerate_states(D)
  S <- nrow(states) # = 2^D

  manifest <- profile_bernoulli_terms(D, class_names, max_order = max_order)
  interaction_terms <- manifest[manifest$order >= 2, , drop = FALSE]

  B <- states # S x D, main-effect design (order 1), order matches class_names via colnames(Y)/arg order

  # concatenate every interaction order's design matrix into one S x P block,
  # and build the matching per-column prior-SD vector
  V_list <- list()
  eta_sd <- numeric(0)
  if (max_order >= 2) {
    for (k in 2:max_order) {
      Vk <- profile_interaction_design(states, k)
      V_list[[length(V_list) + 1]] <- Vk
      sd_k <- if (is.function(interaction_sd)) interaction_sd(k) else interaction_sd
      eta_sd <- c(eta_sd, rep(sd_k, ncol(Vk)))
    }
  }
  V <- if (length(V_list)) do.call(cbind, V_list) else matrix(numeric(0), S, 0)
  P <- ncol(V)
  stopifnot(nrow(interaction_terms) == P)

  if (is.null(beta_sd)) {
    beta_sd <- rep(1.5, K)
    if (K >= 1) beta_sd[1] <- 2.5 # intercept-like first column, matching J2-L's convention
  } else if (length(beta_sd) == 1L) {
    beta_sd <- rep(beta_sd, K)
  }
  stopifnot(length(beta_sd) == K)

  # exact-marginalization bookkeeping, one row per event
  n_compatible <- integer(N)
  # pad unused slots with 1 (always a valid state index) rather than 0 --
  # Stan's declared bound is <lower=1,upper=S> on the *whole* array, even
  # though the model only ever reads the first n_compatible[e] entries of
  # each row; a 0 padding value fails data validation at read time.
  compatible_state <- matrix(1L, N, S)
  for (e in seq_len(N)) {
    idx <- profile_compatible_states(states, Ysign[e, ])
    n_compatible[e] <- length(idx)
    compatible_state[e, seq_along(idx)] <- idx
  }

  list(
    stan_data = list(
      N_events = N, D = D, K = K, S = S, P = P,
      X = X, B = B, V = V,
      n_compatible = n_compatible,
      compatible_state = compatible_state,
      grainsize = grainsize,
      beta_sd = beta_sd,
      eta_sd = eta_sd
    ),
    manifest = manifest,
    design = list(states = states, class_names = class_names, max_order = max_order)
  )
}
