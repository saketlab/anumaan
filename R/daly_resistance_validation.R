# daly_resistance_validation.R
#
# Observed-versus-model and masked-AST validation for a fitted Bayesian
# multivariate probit model (fit_bayesian_multivariate_probit()). These
# functions answer a different question from compute_event_profile_
# probabilities()/aggregate_profiles_for_daly(): "does the fitted model
# reproduce the resistance patterns actually observed in the data?" rather
# than "what is the observed-plus-imputed resistance profile for DALY use?".
# None of the functions here feed into the DALY profile pipeline.
#
# All four checks (marginal, pairwise, complete-profile, masked-AST) use the
# model's own fitted probability Phi(mu_ed) for TESTED cells -- i.e. they
# deliberately do NOT apply the observed-cell-preserving logic from
# compute_event_profile_probabilities(). That preservation is correct for the
# DALY estimand (never overwrite what was actually measured) and wrong here
# (a validation check that can't disagree with the data is not a check).

# ---------------------------------------------------------------------------
# Internal: shared posterior-draw setup for validation functions
# ---------------------------------------------------------------------------
#' @keywords internal
.probit_validation_draws_setup <- function(fitted_model, n_posterior_draws, seed) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)

  set.seed(as.integer(seed))

  draws        <- fitted_model$draws
  class_cols   <- fitted_model$class_cols
  idx_maps     <- fitted_model$index_maps
  event_meta   <- fitted_model$event_metadata
  upper_re_col <- fitted_model$upper_re_col
  pathogen_col <- fitted_model$pathogen_col
  re_prep      <- fitted_model$random_effects_prep

  D <- length(idx_maps$class_levels)
  X_event_sim  <- if (!is.null(fitted_model$X_event)) fitted_model$X_event else fitted_model$X_design
  event_re_idx <- fitted_model$event_re_idx    # N_events x R flattened level index
  K <- ncol(X_event_sim)

  draws_mat <- posterior::as_draws_matrix(draws)
  n_total   <- nrow(draws_mat)
  S <- min(as.integer(n_posterior_draws), n_total)
  draw_idx <- if (S < n_total) sort(sample.int(n_total, S)) else seq_len(n_total)
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  .arr <- function(prefix, d1, d2) {
    cols <- as.vector(outer(seq_len(d1), seq_len(d2),
                            function(a, b) sprintf("%s[%d,%d]", prefix, a, b)))
    array(draws_mat[, cols, drop = FALSE], dim = c(S, d1, d2))
  }
  beta_arr   <- .arr("beta", K, D)
  re_eff_arr <- .arr("re_effect", D, re_prep$total_re_levels)
  residual_structure <- .null_default(fitted_model$residual_structure, "identity")
  L_omega_arr  <- if (identical(residual_structure, "correlated")) .arr("L_Omega", D, D) else NULL
  # Omega is a single global correlation matrix (LKJ prior on L_Omega), the
  # same for every event/hospital at a given posterior draw -- only mu varies
  # by RE block via re_effect. Regularised the same way as
  # compute_event_profile_probabilities() before use in Gibbs/pbivnorm.
  Omega_for_draw <- if (identical(residual_structure, "correlated")) {
    function(s) {
      Om <- tcrossprod(matrix(L_omega_arr[s, , ], nrow = D, ncol = D))
      (Om + t(Om)) / 2 + diag(1e-9, D)
    }
  } else NULL

  stopifnot(".event_idx" %in% names(event_meta))
  has_obs <- event_meta$.event_idx %in% fitted_model$data_long$ev_idx
  event_meta_obs <- event_meta[has_obs, , drop = FALSE]
  N_ev <- nrow(event_meta_obs)
  ev_row_idx <- event_meta_obs$.event_idx
  X_event  <- X_event_sim[ev_row_idx, , drop = FALSE]
  flat_re_idx_obs <- event_re_idx[ev_row_idx, , drop = FALSE]   # N_ev x R

  obs_ast_mat <- as.matrix(event_meta_obs[, class_cols, drop = FALSE])
  storage.mode(obs_ast_mat) <- "double"

  # Generic random-effect contribution: sums over an arbitrary number of
  # declared blocks via the shared re_contribution() helper.
  mu_all_for_draw <- function(s) {
    beta_s   <- matrix(beta_arr[s, , ],   nrow = K, ncol = D)
    re_eff_s <- matrix(re_eff_arr[s, , ], nrow = D, ncol = re_prep$total_re_levels)
    (X_event %*% beta_s) + re_contribution(re_eff_s, flat_re_idx_obs)
  }

  list(
    S = S, N_ev = N_ev, D = D, class_cols = class_cols,
    event_meta_obs = event_meta_obs, obs_ast_mat = obs_ast_mat,
    upper_re_col = upper_re_col, pathogen_col = pathogen_col,
    residual_structure = residual_structure,
    mu_all_for_draw = mu_all_for_draw,
    Omega_for_draw = Omega_for_draw
  )
}

# ---------------------------------------------------------------------------
# validate_marginal_calibration()
# ---------------------------------------------------------------------------

#' Observed-versus-Model Marginal Resistance Validation
#'
#' For every hospital x pathogen x antibiotic-class cell with adequate
#' observed testing support, compares the observed marginal resistance rate
#' against the model's fitted probability \eqn{\Phi(\mu_{ed})}, averaged over
#' the same tested-event cohort and over posterior draws. This is a
#' model-validation calculation: it uses the fitted probability for TESTED
#' cells (unlike \code{compute_event_profile_probabilities()}, which never
#' overwrites an observed value). A calibration residual near 0 and interval
#' coverage near the nominal CI level indicate the model reproduces the
#' observed marginal resistance rates; it does not by itself validate joint
#' (pairwise/profile) structure.
#'
#' @param fitted_model List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param n_posterior_draws_for_validation Integer. Posterior draws used.
#'   Default \code{2000L}.
#' @param seed Integer. Random seed (draw subsampling only). Default \code{123L}.
#' @param ci_level Numeric. Credible interval coverage. Default \code{0.95}.
#' @param min_tested,min_resistant,min_susceptible Integer or \code{NULL}.
#'   Eligibility thresholds. When all \code{NULL} (default), eligibility is
#'   taken from \code{fitted_model$eligibility_report$marginal} (the same
#'   approved rules used for profile panels) rather than re-deriving new
#'   thresholds. Supplying any of these overrides that report and applies the
#'   supplied thresholds instead.
#'
#' @return Tibble, one row per eligible hospital x pathogen x class cell:
#'   \code{n_events}, \code{n_tested}, \code{n_resistant}, \code{n_susceptible},
#'   \code{observed_resistance}, \code{model_resistance_mean/lower/upper},
#'   \code{absolute_error}, \code{calibration_residual},
#'   \code{interval_contains_observed}.
#' @export
validate_marginal_calibration <- function(
    fitted_model,
    n_posterior_draws_for_validation = 2000L,
    seed = 123L,
    ci_level = 0.95,
    min_tested = NULL,
    min_resistant = NULL,
    min_susceptible = NULL
) {
  setup <- .probit_validation_draws_setup(fitted_model, n_posterior_draws_for_validation, seed)
  upper_re_col <- setup$upper_re_col
  pathogen_col <- setup$pathogen_col
  class_cols   <- setup$class_cols
  event_meta_obs <- setup$event_meta_obs
  obs_ast_mat  <- setup$obs_ast_mat
  S <- setup$S

  elig <- fitted_model$eligibility_report$marginal
  use_custom <- !is.null(min_tested) || !is.null(min_resistant) || !is.null(min_susceptible)
  mt <- .null_default(min_tested, 30L)
  mr <- .null_default(min_resistant, 5L)
  ms <- .null_default(min_susceptible, 5L)

  hp_pairs <- unique(event_meta_obs[, c(upper_re_col, pathogen_col), drop = FALSE])

  cell_meta <- list()
  for (r in seq_len(nrow(hp_pairs))) {
    h_nm <- hp_pairs[[upper_re_col]][r]; k_nm <- hp_pairs[[pathogen_col]][r]
    sub_idx <- which(event_meta_obs[[upper_re_col]] == h_nm & event_meta_obs[[pathogen_col]] == k_nm)
    for (d in seq_along(class_cols)) {
      cc <- class_cols[d]
      vals <- obs_ast_mat[sub_idx, d]
      tested_local <- which(!is.na(vals))
      n_tested <- length(tested_local)
      if (n_tested == 0L) next
      n_resistant   <- sum(vals[tested_local] == 1)
      n_susceptible <- sum(vals[tested_local] == 0)
      eligible <- if (use_custom) {
        n_tested >= mt && n_resistant >= mr && n_susceptible >= ms
      } else {
        e <- elig[elig[[upper_re_col]] == h_nm & elig[[pathogen_col]] == k_nm &
                   elig$antibiotic_class == cc, ]
        nrow(e) > 0L && isTRUE(e$eligible[1])
      }
      if (!eligible) next
      key <- paste(h_nm, k_nm, cc, sep = "||")
      cell_meta[[key]] <- list(
        h = h_nm, k = k_nm, cc = cc, d = d,
        ev_idx = sub_idx[tested_local], n_events = length(sub_idx),
        n_tested = n_tested, n_resistant = n_resistant, n_susceptible = n_susceptible
      )
    }
  }

  if (length(cell_meta) == 0L) {
    warning("[validate_marginal_calibration] No hospital x pathogen x class cells meet eligibility.",
            call. = FALSE)
    return(tibble::tibble())
  }

  keys <- names(cell_meta)
  draw_rows <- vector("list", S)
  for (s in seq_len(S)) {
    p_all <- stats::pnorm(setup$mu_all_for_draw(s))
    draw_rows[[s]] <- tibble::tibble(
      key = keys, draw_s = s,
      model_resistance_s = vapply(cell_meta, function(cm) mean(p_all[cm$ev_idx, cm$d]), numeric(1L))
    )
  }
  draws_tbl <- dplyr::bind_rows(draw_rows)

  lo_q <- (1 - ci_level) / 2; hi_q <- 1 - lo_q
  summary_tbl <- draws_tbl %>%
    dplyr::group_by(.data$key) %>%
    dplyr::summarise(
      model_resistance_mean  = mean(.data$model_resistance_s),
      model_resistance_lower = stats::quantile(.data$model_resistance_s, lo_q),
      model_resistance_upper = stats::quantile(.data$model_resistance_s, hi_q),
      .groups = "drop"
    )

  meta_tbl <- dplyr::bind_rows(lapply(keys, function(key) {
    cm <- cell_meta[[key]]
    tibble::tibble(
      key = key,
      !!upper_re_col := cm$h, !!pathogen_col := cm$k, antibiotic_class = cm$cc,
      n_events = cm$n_events, n_tested = cm$n_tested,
      n_resistant = cm$n_resistant, n_susceptible = cm$n_susceptible,
      observed_resistance = cm$n_resistant / cm$n_tested
    )
  }))

  dplyr::left_join(meta_tbl, summary_tbl, by = "key") %>%
    dplyr::mutate(
      absolute_error = abs(.data$observed_resistance - .data$model_resistance_mean),
      calibration_residual = .data$observed_resistance - .data$model_resistance_mean,
      interval_contains_observed =
        .data$observed_resistance >= .data$model_resistance_lower &
        .data$observed_resistance <= .data$model_resistance_upper
    ) %>%
    dplyr::select(-"key") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 6L)))
}

# ---------------------------------------------------------------------------
# validate_pairwise_calibration()
# ---------------------------------------------------------------------------

#' Observed-versus-Model Pairwise Co-resistance Validation
#'
#' For every adequately co-tested hospital x pathogen x class-pair
#' combination, compares the observed pairwise resistance rate
#' (\eqn{n_{RR}/n_{\text{cotested}}}) against the model's implied joint
#' probability, computed differently by residual structure:
#' \describe{
#'   \item{Identity}{\eqn{P(Y_{d_1}=1, Y_{d_2}=1 \mid \theta) =
#'     \Phi(\mu_{ed_1})\,\Phi(\mu_{ed_2})} -- the two classes are combined
#'     only through the shared fixed and random effects in \eqn{\mu}, not
#'     through any residual correlation (the identity model has none).
#'     \strong{This checks whether the identity-residual model reproduces
#'     observed pairwise resistance through shared covariates and random
#'     effects; it is not evidence that a residual-independence assumption
#'     is correct.} Persistent pairwise miscalibration despite good marginal
#'     calibration suggests unmodelled class-to-class dependence that a
#'     correlated-residual fit might capture.}
#'   \item{Correlated (Commit 6)}{\eqn{P(Y_{d_1}=1, Y_{d_2}=1 \mid \theta) =
#'     \Phi_2(\mu_{ed_1}, \mu_{ed_2}; \rho_{d_1 d_2})} -- the standard
#'     bivariate-probit identity \eqn{P(Z_1>0,Z_2>0) = \Phi_2(\mu_1,\mu_2;\rho)}
#'     for \eqn{(Z_1,Z_2)\sim N((\mu_1,\mu_2), [[1,\rho],[\rho,1]])}, with
#'     \eqn{\rho_{d_1 d_2} = \Omega_{d_1 d_2}} from the fitted residual
#'     correlation matrix (computed via \pkg{pbivnorm}). This actually uses
#'     the fitted correlation structure, unlike the identity-model product
#'     formula -- checking whether \eqn{\Omega} reproduces observed
#'     co-resistance.}
#' }
#'
#' @inheritParams validate_marginal_calibration
#' @param min_cotested Integer or \code{NULL}. Eligibility threshold. When
#'   \code{NULL} (default), eligibility is taken from
#'   \code{fitted_model$eligibility_report$pairwise}.
#'
#' @return Tibble, one row per eligible hospital x pathogen x class-pair:
#'   \code{n_cotested}, \code{observed_RR/RS/SR/SS},
#'   \code{observed_pairwise_resistance}, \code{model_pairwise_mean/lower/upper},
#'   \code{absolute_error}, \code{interval_contains_observed}.
#' @export
validate_pairwise_calibration <- function(
    fitted_model,
    n_posterior_draws_for_validation = 2000L,
    seed = 123L,
    ci_level = 0.95,
    min_cotested = NULL
) {
  setup <- .probit_validation_draws_setup(fitted_model, n_posterior_draws_for_validation, seed)
  upper_re_col <- setup$upper_re_col
  pathogen_col <- setup$pathogen_col
  class_cols   <- setup$class_cols
  event_meta_obs <- setup$event_meta_obs
  obs_ast_mat  <- setup$obs_ast_mat
  S <- setup$S

  if (length(class_cols) < 2L) {
    warning("[validate_pairwise_calibration] Fewer than 2 class_cols; nothing to validate.", call. = FALSE)
    return(tibble::tibble())
  }

  pw <- fitted_model$eligibility_report$pairwise
  use_custom <- !is.null(min_cotested)
  mc <- .null_default(min_cotested, 30L)

  hp_pairs <- unique(event_meta_obs[, c(upper_re_col, pathogen_col), drop = FALSE])
  pairs <- utils::combn(class_cols, 2L, simplify = FALSE)

  cell_meta <- list()
  for (r in seq_len(nrow(hp_pairs))) {
    h_nm <- hp_pairs[[upper_re_col]][r]; k_nm <- hp_pairs[[pathogen_col]][r]
    sub_idx <- which(event_meta_obs[[upper_re_col]] == h_nm & event_meta_obs[[pathogen_col]] == k_nm)
    for (pr in pairs) {
      c1 <- pr[1L]; c2 <- pr[2L]
      d1 <- match(c1, class_cols); d2 <- match(c2, class_cols)
      v1 <- obs_ast_mat[sub_idx, d1]; v2 <- obs_ast_mat[sub_idx, d2]
      cotested_local <- which(!is.na(v1) & !is.na(v2))
      n_cotested <- length(cotested_local)
      if (n_cotested == 0L) next
      eligible <- if (use_custom) {
        n_cotested >= mc
      } else {
        e <- pw[pw[[upper_re_col]] == h_nm & pw$class_1 == c1 & pw$class_2 == c2, ]
        nrow(e) > 0L && isTRUE(e$sufficient[1])
      }
      if (!eligible) next
      ev_idx_ct <- sub_idx[cotested_local]
      v1c <- v1[cotested_local]; v2c <- v2[cotested_local]
      key <- paste(h_nm, k_nm, c1, c2, sep = "||")
      cell_meta[[key]] <- list(
        h = h_nm, k = k_nm, c1 = c1, c2 = c2, d1 = d1, d2 = d2,
        ev_idx = ev_idx_ct, n_cotested = n_cotested,
        n_RR = sum(v1c == 1 & v2c == 1), n_RS = sum(v1c == 1 & v2c == 0),
        n_SR = sum(v1c == 0 & v2c == 1), n_SS = sum(v1c == 0 & v2c == 0)
      )
    }
  }

  if (length(cell_meta) == 0L) {
    warning("[validate_pairwise_calibration] No hospital x class-pair combinations meet eligibility.",
            call. = FALSE)
    return(tibble::tibble())
  }

  is_correlated <- identical(setup$residual_structure, "correlated")
  if (is_correlated && !requireNamespace("pbivnorm", quietly = TRUE))
    stop("[validate_pairwise_calibration] Package 'pbivnorm' is required to compute the ",
         "correlated-residual joint probability Phi_2(mu_d1, mu_d2; rho).", call. = FALSE)

  keys <- names(cell_meta)
  draw_rows <- vector("list", S)
  for (s in seq_len(S)) {
    mu_all <- setup$mu_all_for_draw(s)
    if (is_correlated) {
      # Correlated residual: the model-implied joint P(Y_d1=1, Y_d2=1) is
      # Phi_2(mu_d1, mu_d2; rho_d1d2) -- the standard bivariate-probit
      # identity P(Z1>0,Z2>0) = Phi_2(mu1,mu2;rho) for (Z1,Z2) ~ N((mu1,mu2),
      # [[1,rho],[rho,1]]) -- NOT the independence-assumption product used
      # for identity residuals. rho_d1d2 = Omega[d1,d2], the SAME global
      # correlation matrix for every hospital/event at this draw.
      Omega_s <- setup$Omega_for_draw(s)
      draw_rows[[s]] <- tibble::tibble(
        key = keys, draw_s = s,
        model_pairwise_s = vapply(cell_meta, function(cm) {
          rho <- Omega_s[cm$d1, cm$d2]
          mean(pbivnorm::pbivnorm(mu_all[cm$ev_idx, cm$d1], mu_all[cm$ev_idx, cm$d2], rho = rho))
        }, numeric(1L))
      )
    } else {
      p_all <- stats::pnorm(mu_all)
      draw_rows[[s]] <- tibble::tibble(
        key = keys, draw_s = s,
        model_pairwise_s = vapply(cell_meta, function(cm)
          mean(p_all[cm$ev_idx, cm$d1] * p_all[cm$ev_idx, cm$d2]), numeric(1L))
      )
    }
  }
  draws_tbl <- dplyr::bind_rows(draw_rows)

  lo_q <- (1 - ci_level) / 2; hi_q <- 1 - lo_q
  summary_tbl <- draws_tbl %>%
    dplyr::group_by(.data$key) %>%
    dplyr::summarise(
      model_pairwise_mean  = mean(.data$model_pairwise_s),
      model_pairwise_lower = stats::quantile(.data$model_pairwise_s, lo_q),
      model_pairwise_upper = stats::quantile(.data$model_pairwise_s, hi_q),
      .groups = "drop"
    )

  meta_tbl <- dplyr::bind_rows(lapply(keys, function(key) {
    cm <- cell_meta[[key]]
    tibble::tibble(
      key = key,
      !!upper_re_col := cm$h, !!pathogen_col := cm$k,
      class_1 = cm$c1, class_2 = cm$c2, n_cotested = cm$n_cotested,
      observed_RR = cm$n_RR, observed_RS = cm$n_RS, observed_SR = cm$n_SR, observed_SS = cm$n_SS,
      observed_pairwise_resistance = cm$n_RR / cm$n_cotested
    )
  }))

  dplyr::left_join(meta_tbl, summary_tbl, by = "key") %>%
    dplyr::mutate(
      absolute_error = abs(.data$observed_pairwise_resistance - .data$model_pairwise_mean),
      interval_contains_observed =
        .data$observed_pairwise_resistance >= .data$model_pairwise_lower &
        .data$observed_pairwise_resistance <= .data$model_pairwise_upper
    ) %>%
    dplyr::select(-"key") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 6L)))
}

# ---------------------------------------------------------------------------
# validate_complete_profile_calibration()
# ---------------------------------------------------------------------------

#' Observed-versus-Model Complete-Profile Validation
#'
#' For hospital x pathogen panels (from \code{.resolve_profile_class_panel()},
#' the same approved-eligibility panel used for DALY profiles) with at least
#' \code{min_complete_events} events that have \strong{every} panel class
#' actually observed (no imputation involved), compares the empirical
#' complete-profile frequency distribution against the model-implied profile
#' probability distribution for that same event cohort (product of per-class
#' \eqn{\Phi(\mu_{ed})}, as if the classes had not been observed -- i.e. the
#' model's unconditional prediction, deliberately NOT using
#' \code{compute_event_profile_probabilities()}'s observed-cell-preserving
#' logic, since the point here is to check whether the model's predictions
#' agree with what was actually measured). Panels without enough complete
#' events are skipped and the reason is recorded in the \code{status} column
#' rather than silently omitted.
#'
#' @inheritParams validate_marginal_calibration
#' @param min_complete_events Integer. Minimum number of fully-observed-panel
#'   events required to evaluate a hospital-pathogen panel. Default \code{30L}.
#'
#' @return Tibble with one row per hospital x pathogen x profile for evaluated
#'   panels (\code{status = "evaluated"}: \code{empirical_count/frequency},
#'   \code{model_frequency_mean/lower/upper}, \code{absolute_error},
#'   \code{interval_contains_observed}) and one row per skipped panel
#'   (\code{status} starts with \code{"skipped_"}, numeric columns \code{NA}).
#' @export
validate_complete_profile_calibration <- function(
    fitted_model,
    n_posterior_draws_for_validation = 2000L,
    seed = 123L,
    ci_level = 0.95,
    min_complete_events = 30L
) {
  setup <- .probit_validation_draws_setup(fitted_model, n_posterior_draws_for_validation, seed)
  upper_re_col <- setup$upper_re_col
  pathogen_col <- setup$pathogen_col
  class_cols   <- setup$class_cols
  event_meta_obs <- setup$event_meta_obs
  obs_ast_mat  <- setup$obs_ast_mat
  S <- setup$S
  eligibility_report <- fitted_model$eligibility_report
  residual_structure  <- setup$residual_structure

  hp_pairs <- unique(event_meta_obs[, c(upper_re_col, pathogen_col), drop = FALSE])

  skipped_rows <- list()
  eligible_info <- list()
  for (r in seq_len(nrow(hp_pairs))) {
    h_nm <- hp_pairs[[upper_re_col]][r]; k_nm <- hp_pairs[[pathogen_col]][r]
    panel <- .resolve_profile_class_panel(
      class_cols = class_cols, hospital = h_nm, pathogen = k_nm,
      eligibility_report = eligibility_report,
      upper_re_col = upper_re_col, pathogen_col = pathogen_col,
      residual_structure = residual_structure
    )
    key <- paste(h_nm, k_nm, sep = "||")
    if (length(panel$classes) < 1L) {
      skipped_rows[[key]] <- tibble::tibble(
        !!upper_re_col := h_nm, !!pathogen_col := k_nm,
        profile_class_set = NA_character_, n_profile_classes = 0L,
        n_complete_events = 0L, status = "skipped_no_eligible_classes"
      )
      next
    }
    d_idx <- match(panel$classes, class_cols)
    sub_idx <- which(event_meta_obs[[upper_re_col]] == h_nm & event_meta_obs[[pathogen_col]] == k_nm)
    complete_local <- sub_idx[rowSums(is.na(obs_ast_mat[sub_idx, d_idx, drop = FALSE])) == 0L]
    n_complete <- length(complete_local)
    class_set <- paste(panel$classes, collapse = "|")
    if (n_complete < min_complete_events) {
      skipped_rows[[key]] <- tibble::tibble(
        !!upper_re_col := h_nm, !!pathogen_col := k_nm,
        profile_class_set = class_set, n_profile_classes = length(panel$classes),
        n_complete_events = n_complete,
        status = sprintf("skipped_insufficient_complete_events(need>=%d)", min_complete_events)
      )
      next
    }
    eligible_info[[key]] <- list(h = h_nm, k = k_nm, classes = panel$classes, d_idx = d_idx,
                                 ev_idx = complete_local, class_set = class_set, n_complete = n_complete)
  }

  skipped_tbl <- if (length(skipped_rows) > 0L) dplyr::bind_rows(skipped_rows) else tibble::tibble()

  if (length(eligible_info) == 0L) {
    warning("[validate_complete_profile_calibration] No hospital x pathogen panel has enough complete-profile events.",
            call. = FALSE)
    return(skipped_tbl)
  }

  emp_rows <- list()
  for (key in names(eligible_info)) {
    ci <- eligible_info[[key]]
    obs_sub <- obs_ast_mat[ci$ev_idx, ci$d_idx, drop = FALSE]
    labels  <- apply(obs_sub, 1L, function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))
    tab <- table(labels)
    enum_df <- enumerate_binary_profiles(ci$classes)
    emp_rows[[key]] <- tibble::tibble(
      key = key, profile_delta = enum_df$profile_delta,
      empirical_count = vapply(enum_df$profile_delta,
        function(lbl) if (lbl %in% names(tab)) as.integer(tab[[lbl]]) else 0L, integer(1L))
    ) %>% dplyr::mutate(empirical_frequency = .data$empirical_count / ci$n_complete)
  }
  emp_tbl <- dplyr::bind_rows(emp_rows)

  draw_rows <- vector("list", S)
  for (s in seq_len(S)) {
    p_all <- stats::pnorm(setup$mu_all_for_draw(s))
    per_key <- list()
    for (key in names(eligible_info)) {
      ci <- eligible_info[[key]]
      p_sub <- p_all[ci$ev_idx, ci$d_idx, drop = FALSE]
      enum_df <- enumerate_binary_profiles(ci$classes)
      profile_bin <- as.matrix(enum_df[, ci$classes, drop = FALSE])
      n_profiles <- nrow(profile_bin)
      prob_mat <- matrix(1, nrow(p_sub), n_profiles)
      for (d in seq_len(ncol(p_sub))) {
        col_is1 <- profile_bin[, d] == 1L
        f_d <- matrix(NA_real_, nrow(p_sub), n_profiles)
        if (any(col_is1))  f_d[, col_is1]  <- p_sub[, d]
        if (any(!col_is1)) f_d[, !col_is1] <- 1 - p_sub[, d]
        prob_mat <- prob_mat * f_d
      }
      per_key[[key]] <- tibble::tibble(key = key, profile_delta = enum_df$profile_delta,
                                       model_frequency_s = colMeans(prob_mat))
    }
    draw_rows[[s]] <- dplyr::bind_rows(per_key)
  }
  draws_tbl <- dplyr::bind_rows(draw_rows)

  lo_q <- (1 - ci_level) / 2; hi_q <- 1 - lo_q
  model_summary <- draws_tbl %>%
    dplyr::group_by(.data$key, .data$profile_delta) %>%
    dplyr::summarise(
      model_frequency_mean  = mean(.data$model_frequency_s),
      model_frequency_lower = stats::quantile(.data$model_frequency_s, lo_q),
      model_frequency_upper = stats::quantile(.data$model_frequency_s, hi_q),
      .groups = "drop"
    )

  meta_tbl <- dplyr::bind_rows(lapply(names(eligible_info), function(key) {
    ci <- eligible_info[[key]]
    tibble::tibble(
      key = key, !!upper_re_col := ci$h, !!pathogen_col := ci$k,
      profile_class_set = ci$class_set, n_profile_classes = length(ci$classes),
      n_complete_events = ci$n_complete, status = "evaluated"
    )
  }))

  evaluated_tbl <- meta_tbl %>%
    dplyr::left_join(emp_tbl, by = "key") %>%
    dplyr::left_join(model_summary, by = c("key", "profile_delta")) %>%
    dplyr::mutate(
      absolute_error = abs(.data$empirical_frequency - .data$model_frequency_mean),
      interval_contains_observed =
        .data$empirical_frequency >= .data$model_frequency_lower &
        .data$empirical_frequency <= .data$model_frequency_upper
    ) %>%
    dplyr::select(-"key") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 6L)))

  dplyr::bind_rows(skipped_tbl, evaluated_tbl)
}

# ---------------------------------------------------------------------------
# mask_and_validate_ast()
# ---------------------------------------------------------------------------
#' @keywords internal
.auroc <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2L) return(NA_real_)
  r  <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Masked-AST Holdout Validation
#'
#' Masks a reproducible, stratified subset of genuinely \strong{observed} AST
#' cells, predicts those cells' resistance probability from a fitted model,
#' and scores the predictions against their true values.
#'
#' \strong{Two validation modes:}
#' \describe{
#'   \item{\code{refit = FALSE} (default)}{Predicts masked cells using the
#'     \emph{already-fitted} \code{fitted_model}. Because \eqn{\mu_{ed}} never
#'     depends on class \eqn{d}'s own observed value (only on event-level
#'     covariates and random-effect memberships), this is honestly an
#'     \strong{in-sample conditional-prediction diagnostic}
#'     (\code{validation_mode = "in_sample_no_refit"}), not a genuine
#'     held-out test -- the masked cells still contributed to estimating
#'     \code{beta}/\code{hospital_effect} etc. It is useful for spotting
#'     miscalibrated subgroups quickly, but is optimistic relative to a real
#'     holdout.}
#'   \item{\code{refit = TRUE} (preferred, more expensive)}{Sets the selected
#'     cells to \code{NA} in a copy of \code{fitted_model$event_metadata},
#'     re-runs \code{fit_bayesian_multivariate_probit()} on the masked data,
#'     and predicts the held-out cells from that refit
#'     (\code{validation_mode = "refit_masked"}). This is a genuine holdout
#'     and is the recommended final procedure once computationally
#'     affordable.}
#' }
#'
#' \strong{Masking procedure:} reproducible via \code{seed}; stratified by
#' hospital x pathogen x class x observed R/S state; for every hospital x
#' pathogen x class cell, at most enough cells are masked to keep
#' \code{n_tested/n_resistant/n_susceptible} at or above
#' \code{min_tested_after_mask}/\code{min_resistant_after_mask}/
#' \code{min_susceptible_after_mask} (panel eligibility is not destroyed);
#' only observed (non-\code{NA}) cells are ever candidates; and no event is
#' ever left with zero observed classes across \emph{all} \code{class_cols}
#' (which would otherwise silently drop that event from a refit).
#'
#' \strong{Prediction} (Commit 6): for \code{fitted_model$residual_structure
#' == "identity"}, \eqn{\Phi(\mu_{ed})} is the exact model-implied
#' probability for a masked cell regardless of what else is known about that
#' event (classes are conditionally independent given \eqn{\mu}). For
#' \code{"correlated"}, using \eqn{\Phi(\mu_{ed})} alone would ignore every
#' OTHER still-observed class for that event -- wrong whenever classes are
#' correlated. Instead each masked event's full \code{class_cols} panel is
#' rebuilt with the masked cell(s) set to \code{NA} and every other class at
#' its true observed value, and \code{.gibbs_conditional_profile_probs()} is
#' used to get \eqn{P(\text{masked class}=R \mid \text{that event's other
#' observed classes}, \theta)} -- the same conditional-imputation machinery
#' used for DALY profile generation, not a marginal prediction.
#'
#' @param fitted_model List returned by \code{fit_bayesian_multivariate_probit()}.
#' @param fixed_effects Character vector. Required when \code{refit = TRUE}
#'   (not recoverable from \code{fitted_model} due to dummy-coding of
#'   categorical covariates); must match the original fitting call.
#' @param event_id_col Character. Column in \code{fitted_model$event_metadata}
#'   holding the original event id; only used when \code{refit = TRUE}.
#'   Default \code{"event_id"}.
#' @param outcome_col Character or \code{NULL}. Only used when
#'   \code{refit = TRUE}.
#' @param fraction_to_mask Numeric in (0, 1]. Target fraction of observed
#'   cells to mask per hospital x pathogen x class cell, before the
#'   eligibility floor and \code{max_mask_per_cell} caps are applied. Default
#'   \code{0.05}.
#' @param max_mask_per_cell Integer or \code{NULL}. Optional additional cap on
#'   cells masked per hospital x pathogen x class cell.
#' @param seed Integer. Random seed for mask selection (and posterior draw
#'   subsampling). Default \code{123L}.
#' @param refit Logical. See modes above. Default \code{FALSE}.
#' @param n_posterior_draws_for_validation Integer. Default \code{2000L}.
#' @param min_tested_after_mask,min_resistant_after_mask,min_susceptible_after_mask
#'   Integer. Panel-eligibility floors enforced during masking. Default
#'   \code{30L}, \code{5L}, \code{5L}.
#' @param panel_eligibility Named list. Forwarded to
#'   \code{fit_bayesian_multivariate_probit()} when \code{refit = TRUE}.
#' @param prior_config,sampler_config Named list or \code{NULL}. Forwarded
#'   when \code{refit = TRUE}; default to the original fit's
#'   \code{prior_config_used}/\code{sampler_config_used}.
#' @param n_gibbs_burnin,n_gibbs_kept Integer. Gibbs burn-in/kept-iteration
#'   counts forwarded to \code{.gibbs_conditional_profile_probs()} for
#'   correlated-residual fits' masked-cell prediction (Commit 6). Ignored for
#'   identity-residual fits. Default \code{10L}/\code{20L}.
#'
#' @return Named list: \code{predictions} (one row per masked cell: hospital,
#'   pathogen, event id, class, true value, predicted probability, log score,
#'   Brier score, predicted class at 0.5, calibration group, validation
#'   mode), \code{summary} (one-row aggregate: counts, mean log/Brier score,
#'   accuracy, AUROC \[only when both classes are present in the masked set\],
#'   masking configuration), and \code{refit_model} (the refit
#'   \code{fitted_model}, or \code{NULL} when \code{refit = FALSE}).
#' @export
mask_and_validate_ast <- function(
    fitted_model,
    fixed_effects  = NULL,
    event_id_col   = "event_id",
    outcome_col    = NULL,
    fraction_to_mask = 0.05,
    max_mask_per_cell = NULL,
    seed = 123L,
    refit = FALSE,
    n_posterior_draws_for_validation = 2000L,
    min_tested_after_mask      = 30L,
    min_resistant_after_mask   = 5L,
    min_susceptible_after_mask = 5L,
    panel_eligibility = list(),
    prior_config   = NULL,
    sampler_config = NULL,
    n_gibbs_burnin = 10L,
    n_gibbs_kept   = 20L
) {
  class_cols   <- fitted_model$class_cols
  upper_re_col <- fitted_model$upper_re_col
  pathogen_col <- fitted_model$pathogen_col
  event_meta   <- fitted_model$event_metadata

  if (!".event_idx" %in% names(event_meta))
    stop("[mask_and_validate_ast] fitted_model$event_metadata is missing `.event_idx`.", call. = FALSE)

  set.seed(as.integer(seed))

  n_observed_per_event <- rowSums(!is.na(event_meta[, class_cols, drop = FALSE]))

  pool <- lapply(seq_along(class_cols), function(d) {
    cc <- class_cols[d]
    vals <- event_meta[[cc]]
    # Only observed cells, and only where masking would NOT leave the event
    # with zero observed classes (which would silently drop it from a refit).
    obs_idx <- which(!is.na(vals) & n_observed_per_event > 1L)
    if (length(obs_idx) == 0L) return(NULL)
    tibble::tibble(
      row = obs_idx,
      event_idx = event_meta$.event_idx[obs_idx],
      hospital  = event_meta[[upper_re_col]][obs_idx],
      pathogen  = event_meta[[pathogen_col]][obs_idx],
      class_col = cc,
      value     = vals[obs_idx]
    )
  })
  pool_df <- dplyr::bind_rows(pool)
  if (nrow(pool_df) == 0L)
    stop("[mask_and_validate_ast] No maskable cells (every event has only one observed class).", call. = FALSE)

  cell_groups <- split(pool_df, list(pool_df$hospital, pool_df$pathogen, pool_df$class_col), drop = TRUE)
  masked_list <- lapply(cell_groups, function(cell_df) {
    n_total <- nrow(cell_df)
    n_res <- sum(cell_df$value == 1); n_sus <- sum(cell_df$value == 0)
    max_mask_total <- max(0L, n_total - min_tested_after_mask)
    max_mask_res   <- max(0L, n_res   - min_resistant_after_mask)
    max_mask_sus   <- max(0L, n_sus   - min_susceptible_after_mask)

    target_total <- floor(fraction_to_mask * n_total)
    if (!is.null(max_mask_per_cell)) target_total <- min(target_total, max_mask_per_cell)
    target_total <- min(target_total, max_mask_total)
    if (target_total <= 0L) return(NULL)

    res_df <- cell_df[cell_df$value == 1, ]; sus_df <- cell_df[cell_df$value == 0, ]
    target_res <- min(round(target_total * n_res / n_total), max_mask_res, nrow(res_df))
    target_sus <- min(target_total - target_res, max_mask_sus, nrow(sus_df))
    if (target_res + target_sus < target_total) {
      extra <- min(target_total - target_res - target_sus, max_mask_res - target_res, nrow(res_df) - target_res)
      target_res <- target_res + max(0L, extra)
    }
    pick_res <- if (target_res > 0L) res_df[sample.int(nrow(res_df), target_res), ] else res_df[0, ]
    pick_sus <- if (target_sus > 0L) sus_df[sample.int(nrow(sus_df), target_sus), ] else sus_df[0, ]
    dplyr::bind_rows(pick_res, pick_sus)
  })
  masked_df <- dplyr::bind_rows(masked_list)

  if (nrow(masked_df) == 0L)
    stop("[mask_and_validate_ast] No cells could be masked without breaching the eligibility floors -- ",
         "lower fraction_to_mask/min_*_after_mask or check data support.", call. = FALSE)

  message(sprintf("[mask_and_validate_ast] Masking %d observed cell(s) across %d hospital x pathogen x class cell(s).",
                  nrow(masked_df), length(cell_groups)))

  refit_model <- NULL
  if (isTRUE(refit)) {
    if (is.null(fixed_effects))
      stop("[mask_and_validate_ast] `fixed_effects` must be supplied when refit = TRUE ",
           "(not recoverable from fitted_model; dummy-coding of categoricals is ambiguous).", call. = FALSE)
    # Reuse the exact declared block spec (names/group_cols/terms) from the
    # original fit, not just its group-column names, so a refit during masked
    # validation is specified identically to the original model.
    random_effects <- fitted_model$random_effects_prep$blocks

    masked_data <- event_meta
    for (i in seq_len(nrow(masked_df))) {
      masked_data[[masked_df$class_col[i]]][masked_df$row[i]] <- NA
    }

    message(sprintf("[mask_and_validate_ast] Refitting on masked data (%d cell(s) withheld)...", nrow(masked_df)))
    refit_model <- fit_bayesian_multivariate_probit(
      event_class_data   = masked_data,
      class_cols         = class_cols,
      fixed_effects      = fixed_effects,
      random_effects     = random_effects,
      pathogen           = fitted_model$pathogen_fitted,
      pathogen_col       = pathogen_col,
      event_id_col       = event_id_col,
      outcome_col        = outcome_col,
      panel_eligibility  = panel_eligibility,
      residual_structure = .null_default(fitted_model$residual_structure, "identity"),
      estimand           = .null_default(fitted_model$estimand, "observed_stewardship_event_mix"),
      prior_config       = .null_default(prior_config, .null_default(fitted_model$prior_config_used, list())),
      sampler_config     = .null_default(sampler_config, .null_default(fitted_model$sampler_config_used, list())),
      show_messages      = FALSE
    )
    prediction_setup <- .probit_validation_draws_setup(refit_model, n_posterior_draws_for_validation, seed)
    validation_mode <- "refit_masked"
  } else {
    prediction_setup <- .probit_validation_draws_setup(fitted_model, n_posterior_draws_for_validation, seed)
    validation_mode <- "in_sample_no_refit"
  }

  pred_event_row <- match(masked_df$event_idx, prediction_setup$event_meta_obs$.event_idx)
  pred_class_idx <- match(masked_df$class_col, prediction_setup$class_cols)
  if (any(is.na(pred_event_row)) || any(is.na(pred_class_idx)))
    stop("[mask_and_validate_ast] Internal error: masked event/class not found in the prediction model.",
         call. = FALSE)

  S2 <- prediction_setup$S
  is_correlated <- identical(prediction_setup$residual_structure, "correlated")

  if (is_correlated) {
    # Correlated residual (Commit 6): predicting a masked cell from marginal
    # Phi(mu) alone would ignore every OTHER still-observed class for that
    # same event -- wrong whenever those classes are correlated with the
    # masked one. Instead, condition on that event's other observed classes
    # via the same Gibbs machinery used for DALY profile imputation: rebuild
    # each masked event's full class_cols panel with (a) genuinely untested
    # cells and (b) the cell(s) selected for masking set to NA, and (c)
    # every other class left at its true observed value; then read off
    # P(masked class = R) as the marginal of the resulting joint.
    D <- prediction_setup$D
    class_cols <- prediction_setup$class_cols
    masked_event_ids <- unique(masked_df$event_idx)
    masked_event_rows <- match(masked_event_ids, prediction_setup$event_meta_obs$.event_idx)
    obs_panel_pred <- prediction_setup$obs_ast_mat[masked_event_rows, , drop = FALSE]
    row_in_masked <- match(masked_df$event_idx, masked_event_ids)
    for (i in seq_len(nrow(masked_df))) {
      obs_panel_pred[row_in_masked[i], pred_class_idx[i]] <- NA_real_
    }
    enum_full <- enumerate_binary_profiles(class_cols)
    profile_bin_full <- as.matrix(enum_full[, class_cols, drop = FALSE])
    rows_R_by_class <- lapply(seq_len(D), function(d) which(profile_bin_full[, d] == 1L))

    p_draws_mat <- matrix(NA_real_, nrow(masked_df), S2)
    for (s in seq_len(S2)) {
      mu_all <- prediction_setup$mu_all_for_draw(s)
      mu_sub <- mu_all[masked_event_rows, , drop = FALSE]
      Omega_s <- prediction_setup$Omega_for_draw(s)
      prob_mat <- .gibbs_conditional_profile_probs(
        mu_hp = mu_sub, Omega_hp = Omega_s, obs_panel = obs_panel_pred,
        profile_bin = profile_bin_full, n_burnin = n_gibbs_burnin, n_kept = n_gibbs_kept
      )
      for (i in seq_len(nrow(masked_df))) {
        p_draws_mat[i, s] <- sum(prob_mat[row_in_masked[i], rows_R_by_class[[pred_class_idx[i]]]])
      }
    }
  } else {
    p_draws_mat <- matrix(NA_real_, nrow(masked_df), S2)
    for (s in seq_len(S2)) {
      p_all <- stats::pnorm(prediction_setup$mu_all_for_draw(s))
      p_draws_mat[, s] <- p_all[cbind(pred_event_row, pred_class_idx)]
    }
  }
  predicted_p <- rowMeans(p_draws_mat)

  true_value <- masked_df$value
  eps <- 1e-6
  p_clip <- pmin(pmax(predicted_p, eps), 1 - eps)
  log_score   <- -(true_value * log(p_clip) + (1 - true_value) * log(1 - p_clip))
  brier_score <- (predicted_p - true_value)^2
  predicted_class <- as.integer(predicted_p >= 0.5)
  calibration_group <- as.character(cut(predicted_p, breaks = seq(0, 1, 0.1), include.lowest = TRUE))

  predictions <- tibble::tibble(
    !!upper_re_col := masked_df$hospital,
    !!pathogen_col := masked_df$pathogen,
    event_id_value  = masked_df$event_idx,
    antibiotic_class = masked_df$class_col,
    true_ast_value = true_value,
    predicted_resistance_probability = round(predicted_p, 6L),
    log_score   = round(log_score, 6L),
    brier_score = round(brier_score, 6L),
    predicted_class_at_threshold_0_5 = predicted_class,
    calibration_group = calibration_group,
    validation_mode = validation_mode
  )

  summary_tbl <- tibble::tibble(
    validation_mode = validation_mode,
    n_masked = nrow(predictions),
    n_resistant_masked = sum(true_value == 1),
    n_susceptible_masked = sum(true_value == 0),
    mean_log_score = round(mean(log_score), 6L),
    mean_brier_score = round(mean(brier_score), 6L),
    accuracy_at_threshold_0_5 = round(mean(predicted_class == true_value), 6L),
    auroc = round(.auroc(true_value, predicted_p), 6L),
    fraction_to_mask = fraction_to_mask,
    seed = as.integer(seed)
  )

  list(predictions = predictions, summary = summary_tbl, refit_model = refit_model)
}

# ---------------------------------------------------------------------------
# compute_profile_validation_status()
# ---------------------------------------------------------------------------

#' Summarise Validation Checks into a Profile-Validation Status
#'
#' Combines the outputs of \code{validate_marginal_calibration()},
#' \code{validate_pairwise_calibration()},
#' \code{validate_complete_profile_calibration()}, and
#' \code{mask_and_validate_ast()} into a single \code{profile_validation_status}
#' field. \strong{This is deliberately separate from \code{diagnostic_status}}
#' (see \code{fit_bayesian_multivariate_probit()}): \code{diagnostic_status}
#' describes HMC sampler/model-fitting health (R-hat, ESS, divergences,
#' E-BFMI, tree-depth); \code{profile_validation_status} describes predictive
#' agreement with the observed data and data support for that agreement. A
#' model can converge cleanly (\code{diagnostic_status = "pass"}) while its
#' profile predictions are poorly calibrated, or vice versa, and the two
#' should never be merged into one flag.
#'
#' All thresholds are explicit parameters with documented defaults -- none
#' are silently hard-coded elsewhere. Pass \code{thresholds} to override any
#' subset; unspecified entries keep their default.
#'
#' \strong{In-sample masked validation can never produce a bare \code{"pass"}
#' on its own.} When \code{masked_summary$validation_mode ==
#' "in_sample_no_refit"} (see \code{mask_and_validate_ast()}), the masked
#' cells were already seen during fitting, so predicting them well is not
#' genuine validation. If every threshold is otherwise cleared, the status
#' becomes \code{"in_sample_validation_only"} rather than \code{"pass"} --
#' a real \code{"pass"} requires either \code{validation_mode ==
#' "refit_masked"} (a genuine holdout) or no masked-AST check at all (the
#' other observed-versus-model checks being used purely descriptively).
#'
#' @param marginal_tbl Tibble from \code{validate_marginal_calibration()}, or
#'   \code{NULL} if not run.
#' @param pairwise_tbl Tibble from \code{validate_pairwise_calibration()}, or
#'   \code{NULL} if not run.
#' @param complete_profile_tbl Tibble from
#'   \code{validate_complete_profile_calibration()}, or \code{NULL} if not run.
#' @param masked_summary Tibble (\code{$summary} element) from
#'   \code{mask_and_validate_ast()}, or \code{NULL} if not run.
#' @param thresholds Named list overriding any of: \code{max_mean_abs_error_marginal}
#'   (default \code{0.10}), \code{min_marginal_coverage} (default \code{0.80}),
#'   \code{max_mean_abs_error_pairwise} (default \code{0.10}),
#'   \code{min_pairwise_coverage} (default \code{0.80}),
#'   \code{min_complete_profile_panels_evaluated} (default \code{1L}),
#'   \code{max_masked_brier} (default \code{0.25}),
#'   \code{min_masked_auroc} (default \code{0.6}).
#'
#' @return Named list: \code{status} (one of \code{"pass"},
#'   \code{"in_sample_validation_only"},
#'   \code{"warning_marginal_miscalibration"},
#'   \code{"warning_pairwise_miscalibration"},
#'   \code{"warning_sparse_complete_profiles"},
#'   \code{"fail_masked_ast_validation"}, or \code{"not_evaluated"} when no
#'   inputs were supplied), \code{reasons} (character vector of every
#'   triggered condition -- including \code{"masked_validation_in_sample_only"}
#'   as an informational reason even when a more severe status wins --
#'   the most severe status is returned but all are listed),
#'   and \code{thresholds_used} (the resolved threshold list, for audit trail).
#' @export
compute_profile_validation_status <- function(
    marginal_tbl = NULL,
    pairwise_tbl = NULL,
    complete_profile_tbl = NULL,
    masked_summary = NULL,
    thresholds = list()
) {
  th_defaults <- list(
    max_mean_abs_error_marginal = 0.10,
    min_marginal_coverage       = 0.80,
    max_mean_abs_error_pairwise = 0.10,
    min_pairwise_coverage       = 0.80,
    min_complete_profile_panels_evaluated = 1L,
    max_masked_brier = 0.25,
    min_masked_auroc = 0.6
  )
  th <- utils::modifyList(th_defaults, thresholds)

  if (is.null(marginal_tbl) && is.null(pairwise_tbl) &&
      is.null(complete_profile_tbl) && is.null(masked_summary))
    return(list(status = "not_evaluated", reasons = character(0), thresholds_used = th))

  reasons <- character(0)

  if (!is.null(marginal_tbl) && nrow(marginal_tbl) > 0L) {
    mae      <- mean(marginal_tbl$absolute_error, na.rm = TRUE)
    coverage <- mean(marginal_tbl$interval_contains_observed, na.rm = TRUE)
    if ((!is.na(mae) && mae > th$max_mean_abs_error_marginal) ||
        (!is.na(coverage) && coverage < th$min_marginal_coverage))
      reasons <- c(reasons, "warning_marginal_miscalibration")
  }

  if (!is.null(pairwise_tbl) && nrow(pairwise_tbl) > 0L) {
    mae      <- mean(pairwise_tbl$absolute_error, na.rm = TRUE)
    coverage <- mean(pairwise_tbl$interval_contains_observed, na.rm = TRUE)
    if ((!is.na(mae) && mae > th$max_mean_abs_error_pairwise) ||
        (!is.na(coverage) && coverage < th$min_pairwise_coverage))
      reasons <- c(reasons, "warning_pairwise_miscalibration")
  }

  if (!is.null(complete_profile_tbl) && nrow(complete_profile_tbl) > 0L && "status" %in% names(complete_profile_tbl)) {
    # Count distinct evaluated PANELS (hospital x pathogen), not profile rows.
    evald <- complete_profile_tbl[complete_profile_tbl$status == "evaluated", , drop = FALSE]
    panel_cols <- intersect(c("center_name", "hospital", "pathogen", "profile_class_set"), names(evald))
    n_eval <- if (nrow(evald) == 0L || length(panel_cols) == 0L) 0L
              else nrow(unique(evald[, panel_cols, drop = FALSE]))
    if (n_eval < th$min_complete_profile_panels_evaluated)
      reasons <- c(reasons, "warning_sparse_complete_profiles")
  }

  # A model that predicts a cell after already having seen it during fitting
  # (validation_mode == "in_sample_no_refit") is not genuine validation --
  # this must never be allowed to produce a bare "pass" on its own, even when
  # every threshold is technically cleared. It only reflects real held-out
  # performance when validation_mode == "refit_masked".
  masked_in_sample_only <- FALSE
  if (!is.null(masked_summary) && nrow(masked_summary) > 0L) {
    brier <- masked_summary$mean_brier_score[[1]]
    auroc <- masked_summary$auroc[[1]]
    if ((!is.na(brier) && brier > th$max_masked_brier) ||
        (!is.na(auroc) && auroc < th$min_masked_auroc))
      reasons <- c(reasons, "fail_masked_ast_validation")
    masked_in_sample_only <- identical(masked_summary$validation_mode[[1]], "in_sample_no_refit")
    if (masked_in_sample_only)
      reasons <- c(reasons, "masked_validation_in_sample_only")
  }

  reasons <- unique(reasons)
  status <- if (any(grepl("^fail_", reasons))) {
    reasons[grepl("^fail_", reasons)][1]
  } else if (any(grepl("^warning_", reasons))) {
    reasons[grepl("^warning_", reasons)][1]
  } else if (masked_in_sample_only) {
    "in_sample_validation_only"
  } else {
    "pass"
  }

  list(status = status, reasons = reasons, thresholds_used = th)
}
