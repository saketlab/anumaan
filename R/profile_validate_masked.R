# profile_validate_masked.R
# Held-out/masked-cell validation for a joint Bernoulli fit. See
# profile_fit.R's header for how this fits into the pipeline.
#
# Two modes:
#   - profile_validate_masked()  -- predicts masked cells from the ALREADY-
#                                    fitted model (no refit). Honestly an
#                                    in-sample conditional-prediction
#                                    diagnostic, not a genuine holdout, since
#                                    the masked cells still contributed to
#                                    fitting beta/eta.
#   - profile_validate_heldout() -- masks, REFITS on the masked data via
#                                    profile_fit_joint_bernoulli(), and
#                                    predicts the held-out cells from that
#                                    refit. A genuine holdout.
#
# Predicting a masked cell's resistance probability is an EXACT closed-form
# conditional (profile_compatible_states() restricted further to the masked
# class's R sub-states, normalized by the full compatible set) -- the same
# exact marginalization the Stan model's likelihood already uses for
# untested cells, not a simulated or approximated prediction.


#' (Internal) Select a Reproducible, Eligibility-Preserving Mask
#'
#' Stratified masking over a \code{Y} matrix (any \code{D}), by group x
#' class x observed R/S value; for every group x class cell, at most enough cells
#' are masked to keep tested/resistant/susceptible counts at or above the
#' \code{min_*_after_mask} floors; only observed (non-\code{NA}) cells are
#' candidates; no event is ever left with zero observed classes.
#'
#' @keywords internal
.profile_select_mask <- function(Y, ids, fraction_to_mask, max_mask_per_cell,
                                 min_tested_after_mask, min_resistant_after_mask,
                                 min_susceptible_after_mask, seed) {
  D <- ncol(Y)
  classes <- colnames(Y)
  if (is.null(classes)) classes <- paste0("class", seq_len(D))
  n_observed_per_event <- rowSums(!is.na(Y))

  pool <- lapply(seq_len(D), function(d) {
    vals <- Y[, d]
    obs_idx <- which(!is.na(vals) & n_observed_per_event > 1L)
    if (length(obs_idx) == 0L) return(NULL)
    data.frame(
      row = obs_idx, group = ids$group[obs_idx],
      class_idx = d, class = classes[d], value = vals[obs_idx],
      stringsAsFactors = FALSE
    )
  })
  pool_df <- do.call(rbind, pool)
  if (is.null(pool_df) || nrow(pool_df) == 0L) {
    stop("[profile_validate_masked] No maskable cells (every event has only one observed class).", call. = FALSE)
  }

  .local_seed(as.integer(seed))

  cell_groups <- split(pool_df, list(pool_df$group, pool_df$class), drop = TRUE)
  masked_list <- lapply(cell_groups, function(cell_df) {
    n_total <- nrow(cell_df)
    n_res <- sum(cell_df$value == 1)
    n_sus <- sum(cell_df$value == 0)
    max_mask_total <- max(0L, n_total - min_tested_after_mask)
    max_mask_res <- max(0L, n_res - min_resistant_after_mask)
    max_mask_sus <- max(0L, n_sus - min_susceptible_after_mask)

    target_total <- floor(fraction_to_mask * n_total)
    if (!is.null(max_mask_per_cell)) target_total <- min(target_total, max_mask_per_cell)
    target_total <- min(target_total, max_mask_total)
    if (target_total <= 0L) return(NULL)

    res_df <- cell_df[cell_df$value == 1, , drop = FALSE]
    sus_df <- cell_df[cell_df$value == 0, , drop = FALSE]
    target_res <- min(round(target_total * n_res / n_total), max_mask_res, nrow(res_df))
    target_sus <- min(target_total - target_res, max_mask_sus, nrow(sus_df))
    if (target_res + target_sus < target_total) {
      extra <- min(target_total - target_res - target_sus, max_mask_res - target_res, nrow(res_df) - target_res)
      target_res <- target_res + max(0L, extra)
    }
    pick_res <- if (target_res > 0L) res_df[sample.int(nrow(res_df), target_res), , drop = FALSE] else res_df[0, , drop = FALSE]
    pick_sus <- if (target_sus > 0L) sus_df[sample.int(nrow(sus_df), target_sus), , drop = FALSE] else sus_df[0, , drop = FALSE]
    rbind(pick_res, pick_sus)
  })
  masked_df <- do.call(rbind, masked_list)

  if (is.null(masked_df) || nrow(masked_df) == 0L) {
    stop(
      "[profile_validate_masked] No cells could be masked without breaching the eligibility ",
      "floors -- lower fraction_to_mask/min_*_after_mask or check data support.",
      call. = FALSE
    )
  }
  rownames(masked_df) <- NULL
  masked_df
}


#' (Internal) Predict Masked Cells' Resistance Probability
#'
#' For each masked cell, the event's other classes (genuinely untested, OR
#' masked for this same check) are set to sign \code{0}; the posterior
#' predictive probability of the masked class being resistant is
#' \eqn{P(\text{masked}=R \mid \text{other observed}) = \sum_{s \in
#' \text{compatible} \cap R} p_s \big/ \sum_{s \in \text{compatible}} p_s},
#' using \code{\link{profile_compatible_states}} on the event's own sign
#' vector with the masked class(es) opened up -- exact, not simulated.
#' Averaged over posterior draws.
#'
#' @keywords internal
.profile_predict_masked_cells <- function(fit_result, design, Y, masked_df, n_draws, seed) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for profile_validate_masked()/profile_validate_heldout().", call. = FALSE)
  }
  fit <- if (is.list(fit_result) && "fit" %in% names(fit_result)) fit_result$fit else fit_result
  sd_list <- design$stan_data
  D <- sd_list$D
  K <- sd_list$K
  P <- sd_list$P
  S <- sd_list$S
  B <- sd_list$B
  V <- sd_list$V
  X <- sd_list$X
  states <- design$design$states

  draws_mat <- posterior::as_draws_matrix(fit$draws(variables = c("beta", "eta")))
  n_total <- nrow(draws_mat)
  .local_seed(as.integer(seed))
  n_draws <- min(as.integer(n_draws), n_total)
  draw_idx <- sort(sample.int(n_total, n_draws))
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  beta_cols <- sprintf("beta[%d,%d]", rep(seq_len(K), D), rep(seq_len(D), each = K))
  eta_cols <- if (P > 0L) sprintf("eta[%d]", seq_len(P)) else character(0)

  events_to_predict <- unique(masked_df$row)
  pred_prob <- matrix(NA_real_, nrow(masked_df), n_draws)

  for (i in seq_len(n_draws)) {
    beta_s <- matrix(draws_mat[i, beta_cols], nrow = K, ncol = D)
    eta_s <- if (P > 0L) as.numeric(draws_mat[i, eta_cols]) else numeric(0)
    v_eta <- if (P > 0L) as.vector(V %*% eta_s) else rep(0, S)

    for (e in events_to_predict) {
      masked_here <- which(masked_df$row == e)
      sign_vec <- ifelse(is.na(Y[e, ]), 0L, ifelse(Y[e, ] == 1, 1L, -1L))
      sign_vec[masked_df$class_idx[masked_here]] <- 0L # open up the masked cell(s) too

      mu_e <- as.numeric(X[e, ] %*% beta_s) # length D
      score_e <- as.vector(B %*% mu_e) + v_eta # length S
      score_e <- score_e - max(score_e)
      p_e <- exp(score_e)
      p_e <- p_e / sum(p_e)

      compat_idx <- profile_compatible_states(states, sign_vec)
      denom <- sum(p_e[compat_idx])

      for (mh in masked_here) {
        d <- masked_df$class_idx[mh]
        r_idx <- intersect(compat_idx, which(B[, d] == 1))
        pred_prob[mh, i] <- if (denom > 0) sum(p_e[r_idx]) / denom else NA_real_
      }
    }
  }

  rowMeans(pred_prob, na.rm = TRUE)
}


#' (Internal) Average Precision
#'
#' Area under the recall-weighted precision step function, AP = sum_k
#' (Recall_k - Recall_\{k-1\}) * Precision_k.
#' Depends on resistant prevalence (unlike AUROC); undefined, not zero,
#' with zero positive cases. Matches the standard step-function definition
#' (e.g. scikit-learn's \code{average_precision_score}); ties at a
#' threshold are not specially collapsed -- the telescoping sum is
#' invariant to tie order, so this is not an approximation.
#' @keywords internal
.average_precision <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- p[ok]
  n_pos <- sum(y == 1)
  if (n_pos == 0L) return(NA_real_)
  ord <- order(-p)
  y_sorted <- y[ord]
  tp_cum <- cumsum(y_sorted)
  n_seen <- seq_along(y_sorted)
  precision <- tp_cum / n_seen
  recall <- tp_cum / n_pos
  recall_prev <- c(0, recall[-length(recall)])
  sum((recall - recall_prev) * precision)
}


#' (Internal) Log Loss
#'
#' LogLoss = -1/n * sum[y*log(q) + (1-y)*log(1-q)] --
#' lower is better, the opposite sign convention from this file's own
#' per-row \code{log_score} (kept for backward compatibility: higher is
#' better there). Both are reported; do not confuse the two.
#' @keywords internal
.log_loss <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- pmin(pmax(p[ok], 1e-9), 1 - 1e-9)
  if (length(y) == 0L) return(NA_real_)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}


#' (Internal) Binned Calibration Curve
#'
#' Mean predicted probability versus
#' observed frequency, per bin of predicted probability. Reports every
#' non-empty bin, including its own \code{n} -- a thin bin's point is
#' still shown, not silently smoothed over.
#' @keywords internal
.calibration_bins <- function(y, p, n_bins = 10L) {
  y <- as.numeric(y); p <- as.numeric(p)
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- p[ok]
  if (length(y) == 0L) {
    return(data.frame(bin = integer(0), bin_lower = numeric(0), bin_upper = numeric(0),
                       n = integer(0), mean_predicted = numeric(0), observed_frequency = numeric(0)))
  }
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  bin <- cut(p, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  do.call(rbind, lapply(sort(unique(bin)), function(b) {
    idx <- bin == b
    data.frame(
      bin = b, bin_lower = breaks[b], bin_upper = breaks[b + 1L],
      n = sum(idx), mean_predicted = mean(p[idx]), observed_frequency = mean(y[idx]),
      stringsAsFactors = FALSE
    )
  }))
}


#' (Internal) Calibration Intercept/Slope
#'
#' From a logistic regression of
#' truth on logit(predicted probability) -- slope 1/intercept 0 is perfect
#' calibration. NA when unestimable (one-class truth, separation, glm
#' non-convergence) rather than a fabricated number.
#' @keywords internal
.calibration_slope_intercept <- function(y, p) {
  y <- as.numeric(y); p <- as.numeric(p)
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- pmin(pmax(p[ok], 1e-6), 1 - 1e-6)
  if (length(unique(y)) < 2L) return(c(intercept = NA_real_, slope = NA_real_))
  logit_p <- log(p / (1 - p))
  fit <- tryCatch(stats::glm(y ~ logit_p, family = stats::binomial()), error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(c(intercept = NA_real_, slope = NA_real_))
  co <- stats::coef(fit)
  c(intercept = unname(co[1]), slope = unname(co[2]))
}


#' (Internal) Score Masked-Cell Predictions
#' @keywords internal
.profile_score_mask_predictions <- function(masked_df, predicted_resistance_probability, validation_mode) {
  y <- masked_df$value
  p <- pmin(pmax(predicted_resistance_probability, 1e-9), 1 - 1e-9)
  log_score <- ifelse(y == 1, log(p), log(1 - p))
  brier_score <- (y - p)^2
  predicted_class <- as.integer(p >= 0.5)

  predictions <- data.frame(
    group = masked_df$group, event_row = masked_df$row, class = masked_df$class,
    true_value = y, predicted_resistance_probability = p,
    predicted_class = predicted_class, log_score = log_score, brier_score = brier_score,
    validation_mode = validation_mode,
    stringsAsFactors = FALSE
  )

  cal <- .calibration_slope_intercept(y, p)
  summary <- data.frame(
    n_masked = nrow(predictions),
    n_resistant_masked = sum(y == 1),
    n_susceptible_masked = sum(y == 0),
    mean_log_score = mean(log_score),
    log_loss = .log_loss(y, p),
    mean_brier_score = mean(brier_score),
    accuracy = mean(predicted_class == y),
    auroc = .auroc(y, p),
    average_precision = .average_precision(y, p),
    calibration_intercept = cal[["intercept"]],
    calibration_slope = cal[["slope"]],
    validation_mode = validation_mode,
    stringsAsFactors = FALSE
  )

  list(predictions = predictions, summary = summary, calibration = .calibration_bins(y, p))
}


#' Masked-Cell Validation (In-Sample, No Refit)
#'
#' Masks a reproducible, eligibility-preserving subset of genuinely
#' observed cells in \code{Y} and predicts them from the
#' \strong{already-fitted} model (\code{fit_result}) -- honestly an
#' in-sample conditional-prediction diagnostic
#' (\code{validation_mode = "in_sample_no_refit"}), since the masked cells
#' still contributed to fitting \code{beta}/\code{eta}, not a genuine
#' held-out test. Useful for a quick miscalibrated-subgroup scan; prefer
#' \code{\link{profile_validate_heldout}} for the real holdout.
#'
#' @inheritParams profile_ppc
#' @param fraction_to_mask Numeric in (0, 1]. Target fraction of observed
#'   cells to mask per group x class cell, before the eligibility floor and
#'   \code{max_mask_per_cell} cap. Default \code{0.05}.
#' @param max_mask_per_cell Integer or \code{NULL}. Optional additional cap.
#' @param min_tested_after_mask,min_resistant_after_mask,min_susceptible_after_mask
#'   Integer. Panel-eligibility floors enforced during masking. Default
#'   \code{10L}, \code{3L}, \code{3L}.
#' @param seed Integer. Mask selection AND posterior-draw-subsampling seed.
#'   Default \code{123L}.
#'
#' @return Named list: \code{predictions} (one row per masked cell:
#'   \code{group}, \code{event_row}, \code{class}, \code{true_value},
#'   \code{predicted_resistance_probability}, \code{predicted_class},
#'   \code{log_score}, \code{brier_score}, \code{validation_mode}),
#'   \code{summary} (one-row aggregate: counts, mean log/Brier score,
#'   \code{log_loss} (the LogLoss convention, lower-is-better -- the
#'   NEGATIVE of \code{mean_log_score}, both reported so the sign is never
#'   ambiguous), accuracy, \code{auroc} (\code{NA}, not 0, with one-class
#'   truth), \code{average_precision}, \code{calibration_intercept}/
#'   \code{calibration_slope} (\code{NA} when unestimable)), and
#'   \code{calibration} (binned predicted-vs-observed calibration curve).
#' @export
profile_validate_masked <- function(fit_result, design, Y, ids,
                                    fraction_to_mask = 0.05,
                                    max_mask_per_cell = NULL,
                                    min_tested_after_mask = 10L,
                                    min_resistant_after_mask = 3L,
                                    min_susceptible_after_mask = 3L,
                                    n_draws = 500L,
                                    seed = 123L) {
  Y <- as.matrix(Y)
  masked_df <- .profile_select_mask(
    Y, ids, fraction_to_mask, max_mask_per_cell,
    min_tested_after_mask, min_resistant_after_mask, min_susceptible_after_mask, seed
  )
  message(sprintf(
    "[profile_validate_masked] Masking %d observed cell(s) across %d group x class cell(s).",
    nrow(masked_df), length(unique(paste(masked_df$group, masked_df$class)))
  ))

  predicted <- .profile_predict_masked_cells(fit_result, design, Y, masked_df, n_draws, seed)
  .profile_score_mask_predictions(masked_df, predicted, "in_sample_no_refit")
}


#' Masked-Cell Validation (Genuine Holdout, Refits the Model)
#'
#' Masks a reproducible, eligibility-preserving subset of genuinely
#' observed cells, sets them to \code{NA} in a copy of \code{Y}, rebuilds
#' the Stan design on the masked data via \code{\link{profile_build_design}},
#' refits via \code{\link{profile_fit_joint_bernoulli}}, and predicts the
#' held-out cells from THAT refit
#' (\code{validation_mode = "refit_masked"}). This is the preferred, more
#' expensive procedure once computationally affordable; see
#' \code{\link{profile_validate_masked}} for the cheaper in-sample
#' diagnostic.
#'
#' @inheritParams profile_validate_masked
#' @param X The SAME covariate design matrix used for the original fit
#'   (\code{nrow(X) == nrow(Y)}) -- rebuilt against the masked \code{Y} via
#'   \code{\link{profile_build_design}}.
#' @param output_dir Character. Passed to
#'   \code{\link{profile_fit_joint_bernoulli}} for the refit's artifacts.
#' @param sampler_config Named list. Passed to
#'   \code{\link{profile_fit_joint_bernoulli}} for the refit. Default
#'   \code{list()} (that function's own defaults).
#'
#' @return Named list: \code{predictions}, \code{summary} (same shape as
#'   \code{\link{profile_validate_masked}}, \code{validation_mode =
#'   "refit_masked"}), and \code{refit_result} (the full return value of
#'   \code{\link{profile_fit_joint_bernoulli}} on the masked data, so the
#'   refit's own sampler diagnostics can be checked via
#'   \code{\link{profile_check_sampler}} before trusting these scores).
#' @export
profile_validate_heldout <- function(Y, X, ids, output_dir,
                                     fraction_to_mask = 0.05,
                                     max_mask_per_cell = NULL,
                                     min_tested_after_mask = 10L,
                                     min_resistant_after_mask = 3L,
                                     min_susceptible_after_mask = 3L,
                                     sampler_config = list(),
                                     n_draws = 500L,
                                     seed = 123L) {
  Y <- as.matrix(Y)
  masked_df <- .profile_select_mask(
    Y, ids, fraction_to_mask, max_mask_per_cell,
    min_tested_after_mask, min_resistant_after_mask, min_susceptible_after_mask, seed
  )
  message(sprintf(
    "[profile_validate_heldout] Masking %d observed cell(s) across %d group x class cell(s); refitting...",
    nrow(masked_df), length(unique(paste(masked_df$group, masked_df$class)))
  ))

  Y_masked <- Y
  for (i in seq_len(nrow(masked_df))) {
    Y_masked[masked_df$row[i], masked_df$class_idx[i]] <- NA
  }

  design_masked <- profile_build_design(Y = Y_masked, X = X)
  refit_result <- profile_fit_joint_bernoulli(
    design = design_masked, output_dir = output_dir, sampler_config = sampler_config
  )

  predicted <- .profile_predict_masked_cells(refit_result, design_masked, Y_masked, masked_df, n_draws, seed)
  result <- .profile_score_mask_predictions(masked_df, predicted, "refit_masked")
  result$refit_result <- refit_result
  result
}
