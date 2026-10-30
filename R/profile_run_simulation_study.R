# profile_run_simulation_study.R
# SBC + recovery/robustness simulation harness. See profile_fit.R's
# header for how this fits into the pipeline.
#
# Two DELIBERATELY SEPARATE modes, never conflated:
#   - "sbc": simulate parameters FROM THE DECLARED PRIOR, simulate data
#     FROM THE MODEL given those parameters, refit, check whether each
#     true parameter's RANK among its own posterior draws is uniform
#     across many simulated datasets. Tests the COMPUTATION (is the
#     sampler correctly inverting this exact model), NOT whether the
#     model describes real surveillance data.
#   - "recovery": simulate from a CALLER-DECLARED scenario (fixed true
#     parameters, optionally a declared missingness/selection pattern --
#     which may deliberately misspecify relative to the fitted model),
#     refit, report bias/coverage/width/failure rate. Tests practical
#     estimation quality under a named scenario.
#
# Each simulated dataset is its own real Stan fit via
# profile_fit_joint_bernoulli() -- genuinely expensive. n_sim and
# sampler_config are the caller's choice, never a package default sized
# for "enough" (500 repetitions is not universally sufficient -- Section
# D explicitly works out the MCSE this way: ~0.90 assurance from 500 reps
# has binomial MC SE ~0.0134; size n_sim to the MCSE you actually need).


#' (Internal) One Simulate-Then-Refit Replicate
#' @keywords internal
.profile_simulate_and_refit <- function(design, true_beta, true_eta, missingness_fn,
                                        output_dir, sampler_config, sim_id) {
  sd_list <- design$stan_data
  D <- sd_list$D; K <- sd_list$K; P <- sd_list$P; S <- sd_list$S
  B <- sd_list$B; V <- sd_list$V; X <- sd_list$X
  N <- nrow(X)

  Mu <- X %*% true_beta
  Score <- Mu %*% t(B)
  v_eta <- if (P > 0L) as.vector(V %*% true_eta) else rep(0, S)
  Score <- sweep(Score, 2, v_eta, "+")
  Score <- Score - apply(Score, 1, max)
  Pmat <- exp(Score)
  Pmat <- Pmat / rowSums(Pmat)

  cumP <- t(apply(Pmat, 1, cumsum))
  u <- stats::runif(N)
  sim_state <- rowSums(cumP < u) + 1L
  sim_signs <- B[sim_state, , drop = FALSE]
  Y_sim <- matrix(ifelse(sim_signs == 1, 1, 0), N, D)
  colnames(Y_sim) <- design$design$class_names

  if (!is.null(missingness_fn)) Y_sim <- missingness_fn(Y_sim)

  design_sim <- tryCatch(
    profile_build_design(Y_sim, X, class_names = design$design$class_names, max_order = design$design$max_order,
                          beta_sd = sd_list$beta_sd, interaction_sd = sd_list$eta_sd, grainsize = sd_list$grainsize),
    error = function(e) e
  )
  if (inherits(design_sim, "error")) {
    return(list(sim_id = sim_id, ok = FALSE, failure_reason = paste("design:", conditionMessage(design_sim)), fit = NULL, design_sim = NULL))
  }

  fit_result <- tryCatch(
    profile_fit_joint_bernoulli(design_sim, output_dir = file.path(output_dir, sprintf("sim_%04d", sim_id)), sampler_config = sampler_config),
    error = function(e) e
  )
  if (inherits(fit_result, "error")) {
    return(list(sim_id = sim_id, ok = FALSE, failure_reason = paste("fit:", conditionMessage(fit_result)), fit = NULL, design_sim = NULL))
  }

  list(sim_id = sim_id, ok = TRUE, failure_reason = NA_character_, fit = fit_result, design_sim = design_sim)
}


#' SBC + Recovery/Robustness Simulation Harness
#'
#' @param design A template design from \code{\link{profile_build_design}}
#'   (its \code{X}, \code{B}, \code{V}, \code{beta_sd}, \code{eta_sd},
#'   \code{D}, \code{K}, \code{P}, \code{S} are the structural skeleton
#'   every simulated dataset reuses -- only \code{Y} is re-simulated).
#' @param n_sim Integer. Number of simulated datasets to fit. Each is a
#'   SEPARATE, real \code{\link{profile_fit_joint_bernoulli}} run --
#'   genuinely expensive; size it to your compute budget and to the
#'   Monte Carlo precision you actually need (not a universal default).
#' @param mode Character, \code{"sbc"} or \code{"recovery"}.
#' @param true_params (mode \code{"recovery"} only) list with \code{beta}
#'   (\code{K x D} matrix) and \code{eta} (length-\code{P} vector) --
#'   FIXED true parameters every simulated dataset uses. Required for
#'   \code{"recovery"}; must be \code{NULL} for \code{"sbc"} (which draws
#'   its own true parameters from the declared prior, fresh per dataset).
#' @param missingness_fn Function or \code{NULL}. \code{function(Y)} ->
#'   \code{Y} with some entries set to \code{NA}, applied to each
#'   simulated COMPLETE dataset before fitting -- for declaring a sparse-
#'   mask or selective-testing scenario explicitly. Default \code{NULL}
#'   (every simulated dataset is fully observed).
#' @param sampler_config,output_dir Passed to every simulated fit.
#' @param ci_level Numeric. Default \code{0.95}. Used for
#'   \code{mode = "recovery"}'s coverage check.
#' @param seed Integer.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{mode}}{echoed}
#'     \item{\code{results}}{one row per simulated dataset, INCLUDING
#'       failed ones (never silently dropped from the denominator):
#'       \code{sim_id}, \code{ok}, \code{failure_reason},
#'       \code{sampler_status} (from \code{\link{profile_check_sampler}},
#'       \code{NA} if the fit itself failed). For \code{mode = "sbc"}:
#'       one \code{rank_<param>} column per \code{beta}/\code{eta}
#'       parameter (an integer in \code{0..n_draws_used}, \code{NA} if
#'       not computable) and \code{n_draws_used} (after the ESS-based
#'       thinning mitigation for correlated draws). For
#'       \code{mode = "recovery"}: one \code{bias_<param>},
#'       \code{covered_<param>}, \code{width_<param>} triple per
#'       parameter.}
#'     \item{\code{summary}}{\code{mode = "sbc"}: a uniformity check per
#'       parameter (chi-squared goodness-of-fit of the rank histogram, a
#'       descriptive flag, NOT a formal hypothesis test with a universal
#'       threshold). \code{mode = "recovery"}: mean bias, empirical
#'       coverage, mean width, and failure RATE (over
#'       \code{n_sim}, not over converged fits only) per parameter.}
#'   }
#' @export
profile_run_simulation_study <- function(design, n_sim, mode = c("sbc", "recovery"),
                                         true_params = NULL, missingness_fn = NULL,
                                         sampler_config = list(), output_dir,
                                         ci_level = 0.95, seed = 123L) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required for profile_run_simulation_study().", call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for profile_run_simulation_study().", call. = FALSE)
  }
  mode <- match.arg(mode)
  sd_list <- design$stan_data
  D <- sd_list$D; K <- sd_list$K; P <- sd_list$P
  beta_sd <- sd_list$beta_sd; eta_sd <- sd_list$eta_sd

  if (mode == "sbc" && !is.null(true_params)) {
    stop("true_params must be NULL for mode = \"sbc\" -- true parameters are drawn fresh from the prior for every simulated dataset, not fixed.", call. = FALSE)
  }
  if (mode == "recovery" && is.null(true_params)) {
    stop("true_params (list(beta=, eta=)) is required for mode = \"recovery\".", call. = FALSE)
  }

  .local_seed(as.integer(seed))
  param_names <- c(sprintf("beta[%d,%d]", rep(seq_len(K), D), rep(seq_len(D), each = K)),
                    if (P > 0L) sprintf("eta[%d]", seq_len(P)) else character(0))

  rows <- list()
  for (i in seq_len(n_sim)) {
    if (mode == "sbc") {
      true_beta <- matrix(stats::rnorm(K * D, 0, beta_sd), K, D)
      true_eta <- if (P > 0L) stats::rnorm(P, 0, eta_sd) else numeric(0)
    } else {
      true_beta <- true_params$beta
      true_eta <- if (P > 0L) true_params$eta else numeric(0)
    }

    sim <- .profile_simulate_and_refit(design, true_beta, true_eta, missingness_fn, output_dir, sampler_config, i)

    if (!sim$ok) {
      row <- list(sim_id = i, ok = FALSE, failure_reason = sim$failure_reason, sampler_status = NA_character_)
      if (mode == "sbc") {
        row <- c(row, stats::setNames(as.list(rep(NA_real_, length(param_names))), paste0("rank_", param_names)))
        row$n_draws_used <- NA_integer_
      } else {
        for (p in param_names) {
          row[[paste0("bias_", p)]] <- NA_real_
          row[[paste0("covered_", p)]] <- NA
          row[[paste0("width_", p)]] <- NA_real_
        }
      }
      rows[[i]] <- as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
      next
    }

    sampler_status <- tryCatch(profile_check_sampler(sim$fit)$status, error = function(e) NA_character_)
    draws_mat <- posterior::as_draws_matrix(sim$fit$fit$draws(variables = c("beta", "eta")))

    if (mode == "sbc") {
      true_vec <- c(as.vector(true_beta), true_eta)
      rank_row <- list(sim_id = i, ok = TRUE, failure_reason = NA_character_, sampler_status = sampler_status)
      n_used_report <- NA_integer_
      for (j in seq_along(param_names)) {
        x <- draws_mat[, param_names[j]]
        ess <- tryCatch(posterior::ess_basic(x), error = function(e) length(x))
        n_use <- max(10L, min(length(x), round(ess)))
        thin_idx <- round(seq(1, length(x), length.out = n_use))
        x_thin <- x[thin_idx]
        rank_row[[paste0("rank_", param_names[j])]] <- sum(x_thin < true_vec[j])
        n_used_report <- length(x_thin)
      }
      rank_row$n_draws_used <- n_used_report
      rows[[i]] <- as.data.frame(rank_row, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      true_vec <- c(as.vector(true_beta), true_eta)
      lo_q <- (1 - ci_level) / 2; hi_q <- 1 - lo_q
      rec_row <- list(sim_id = i, ok = TRUE, failure_reason = NA_character_, sampler_status = sampler_status)
      for (j in seq_along(param_names)) {
        x <- draws_mat[, param_names[j]]
        est_mean <- mean(x)
        lo <- stats::quantile(x, lo_q, names = FALSE); hi <- stats::quantile(x, hi_q, names = FALSE)
        rec_row[[paste0("bias_", param_names[j])]] <- est_mean - true_vec[j]
        rec_row[[paste0("covered_", param_names[j])]] <- true_vec[j] >= lo & true_vec[j] <= hi
        rec_row[[paste0("width_", param_names[j])]] <- hi - lo
      }
      rows[[i]] <- as.data.frame(rec_row, stringsAsFactors = FALSE, check.names = FALSE)
    }
  }
  results <- do.call(rbind, rows)
  rownames(results) <- NULL

  if (mode == "sbc") {
    summ_rows <- lapply(param_names, function(p) {
      col <- results[[paste0("rank_", p)]]
      col <- col[!is.na(col)]
      if (length(col) < 10L) {
        return(data.frame(parameter = p, n_valid = length(col), chisq_stat = NA_real_, chisq_p_value = NA_real_, stringsAsFactors = FALSE))
      }
      n_bins <- min(10L, length(col))
      max_rank <- max(col, na.rm = TRUE)
      breaks <- unique(round(seq(0, max_rank, length.out = n_bins + 1L)))
      binned <- cut(col, breaks = breaks, include.lowest = TRUE)
      obs <- table(binned)
      chisq <- suppressWarnings(stats::chisq.test(obs))
      data.frame(parameter = p, n_valid = length(col), chisq_stat = unname(chisq$statistic), chisq_p_value = chisq$p.value, stringsAsFactors = FALSE)
    })
    summary_tbl <- do.call(rbind, summ_rows)
  } else {
    summ_rows <- lapply(param_names, function(p) {
      bias_col <- results[[paste0("bias_", p)]]
      cov_col <- results[[paste0("covered_", p)]]
      width_col <- results[[paste0("width_", p)]]
      data.frame(
        parameter = p,
        n_sim = n_sim, n_failed = sum(!results$ok),
        failure_rate = mean(!results$ok),
        mean_bias = mean(bias_col, na.rm = TRUE),
        empirical_coverage = mean(cov_col, na.rm = TRUE),
        mean_width = mean(width_col, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
    summary_tbl <- do.call(rbind, summ_rows)
  }

  list(mode = mode, results = results, summary = summary_tbl)
}
