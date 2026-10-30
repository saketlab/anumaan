# profile_ppc.R
# Posterior predictive checks for a profile_fit_joint_bernoulli() fit.
# See profile_fit.R's header for how this fits into the pipeline.
#
# Six PPC families, all from ONE genuine joint replication per draw (one
# sampled complete state per event, original observation mask reapplied --
# never independent per-class/per-level Bernoulli sampling, which would
# discard exactly the dependence being checked):
#   - marginal            (per class)
#   - pairwise            (per class pair, all 4 RR/RS/SR/SS cells)
#   - complete_profile    (per full 2^D state, complete events only)
#   - count_groups        (0..D resistance-count distribution, complete events only)
#   - partial_patterns    (sub-profile breakdown per EXACT test mask, not
#                          just co-tested pairs -- the denominator is every
#                          event sharing that exact mask)
#   - covariate_subgroups (same marginal check, stratified by a
#                          prespecified clinical covariate instead of group)


#' Posterior Predictive Checks for a Joint Bernoulli Fit
#'
#' Observed-vs-model comparison at three levels -- marginal (per class),
#' pairwise (per class pair, all four RR/RS/SR/SS cells), and
#' complete-profile (per full \code{2^D} state) -- plus a
#' \code{"not_computable"} status for any group x check cell with too
#' little direct support to evaluate at all. Nothing below the relevant
#' \code{min_*} floor is silently omitted; it is reported with
#' \code{status = "not_computable"} and \code{NA} numeric columns.
#'
#' TWO different quantities are reported for every evaluated cell, and
#' they answer different questions -- do not treat them interchangeably:
#' \describe{
#'   \item{\code{expected_mean}/\code{expected_lower}/\code{expected_upper}
#'     (also exposed as \code{model_mean}/\code{model_lower}/
#'     \code{model_upper} for backward compatibility)}{the posterior mean
#'     and credible interval of the model's EXPECTED probability -- an
#'     exact closed-form sum over state probabilities (the same
#'     softmax-over-states construction \code{\link{profile_build_design}}/
#'     the Stan model use for the likelihood). This interval reflects
#'     PARAMETER uncertainty only: a model that is certain \code{p = 0.5}
#'     reports a near-zero-width interval here, even though the actual
#'     count among replicated events would still vary.}
#'   \item{\code{predictive_mean}/\code{predictive_lower}/
#'     \code{predictive_upper}}{a genuine posterior PREDICTIVE interval,
#'     simulated per posterior draw: for every event, one complete profile
#'     is sampled from that draw's joint state distribution, the event's
#'     ORIGINAL observation mask is reapplied (only classes actually
#'     tested in the real data contribute a replicated value), and the
#'     SAME replicated profile is used to compute the replicated marginal,
#'     pairwise, and complete-profile statistics for that event within
#'     that draw -- so the joint dependence the model was fit to capture
#'     is preserved across all three levels, not simulated independently
#'     per level. This interval reflects BOTH parameter uncertainty and
#'     genuine outcome-level sampling variation, per Stan's own
#'     posterior-predictive-vs-posterior-expectation distinction
#'     (\url{https://mc-stan.org/docs/stan-users-guide/posterior-prediction.html}).
#'   }
#' }
#' \code{interval_contains_observed} (unchanged from earlier versions of
#' this function) checks against the EXPECTED interval, for backward
#' compatibility with \code{\link{profile_release}}'s rollup;
#' \code{predictive_interval_contains_observed} checks against the
#' PREDICTIVE interval -- the two can disagree, and a cell that fails the
#' tighter expected-probability check may still pass the (necessarily
#' wider) predictive check.
#'
#' This also avoids the confirmed defect in the J2-L four-cell-sensitivity
#' experiment's \code{fc_partial_specs()} (grouping comparison rows by the
#' observed outcome pattern itself, which makes "observed proportion" 1.0
#' by construction, not a real check): every "observed" count here --
#' \code{n_tested}/\code{n_resistant}, \code{n_cotested}/\code{n_RR} etc.,
#' empirical state frequency -- is a plain tabulation from \code{Y}
#' against a FIXED category (a class, a pair, a state), never derived by
#' matching each row's own label back to itself.
#'
#' @param fit_result List from \code{\link{profile_fit_joint_bernoulli}}
#'   (reads \code{$fit}).
#' @param design List from \code{\link{profile_build_design}} (reads
#'   \code{$stan_data}: \code{X}, \code{B}, \code{V}).
#' @param Y,ids The SAME \code{Y} matrix and \code{ids} data frame passed
#'   to \code{\link{profile_build_design}} (via
#'   \code{\link{profile_prepare_model_data}}) -- the observed side of
#'   every comparison comes from these, not re-derived from raw data.
#' @param n_draws Integer. Posterior draws used (subsampled if the fit has
#'   more). Default \code{500L}.
#' @param min_tested,min_cotested,min_complete Integer. Below these, a
#'   group x check cell is \code{"not_computable"} rather than evaluated.
#'   Default \code{10L} each.
#' @param min_count_group Integer. Below this many COMPLETE events in a
#'   group, the \code{count_groups} family is \code{"not_computable"} for
#'   that group. Default \code{10L}.
#' @param min_partial_pattern Integer. Below this many events sharing an
#'   EXACT partial test mask (in a group), that (group, mask) combination
#'   is skipped from \code{partial_patterns} entirely (never fabricated as
#'   a zero-support row -- there are too many possible masks to report
#'   every absent one; only masks actually meeting this floor are listed,
#'   with every sub-profile of a LISTED mask reported, including its
#'   zero-count ones). Default \code{10L}.
#' @param covariate_subgroups Named list of vectors/factors, or \code{NULL}
#'   (default). Each vector must have length \code{nrow(Y)}, in the SAME
#'   row order as \code{Y}/\code{ids} -- e.g.
#'   \code{list(specimen_type = ..., age_band = ...)}. When supplied, adds
#'   the \code{covariate_subgroups} table: the SAME marginal-level check as
#'   \code{marginal}, stratified by each named vector's levels instead of
#'   \code{group} (pooled across groups -- these are prespecified clinical
#'   subgroup checks, not another hospital breakdown). Default \code{NULL}
#'   (the returned \code{covariate_subgroups} table has zero rows).
#' @param min_subgroup Integer. Below this many tested events in a
#'   subgroup level x class cell, that cell is \code{"not_computable"}.
#'   Default \code{10L}.
#' @param ci_level Numeric. Default \code{0.95}.
#' @param seed Integer. Draw-subsampling AND replicated-outcome simulation
#'   seed. Default \code{123L}.
#'
#' @return A list of data frames:
#'   \describe{
#'     \item{\code{marginal}}{one row per group x class: \code{status},
#'       \code{n_tested}, \code{n_resistant}, \code{n_susceptible},
#'       \code{observed_resistance}, \code{expected_mean/lower/upper}
#'       (= \code{model_mean/lower/upper}), \code{predictive_mean/lower/upper},
#'       \code{absolute_error}, \code{interval_contains_observed},
#'       \code{predictive_interval_contains_observed}}
#'     \item{\code{pairwise}}{one row per group x class-pair:
#'       \code{status}, \code{n_cotested}, \code{n_RR/n_RS/n_SR/n_SS},
#'       \code{observed_RR_prop} etc., \code{expected_RR_mean/lower/upper}
#'       (= \code{model_RR_mean/lower/upper}, and \code{_RS}/\code{_SR}/
#'       \code{_SS}), \code{predictive_RR_mean/lower/upper},
#'       \code{absolute_error} (RR cell), \code{interval_contains_observed},
#'       \code{predictive_interval_contains_observed}}
#'     \item{\code{complete_profile}}{one row per group x state (one
#'       \code{"not_computable"} row per group below \code{min_complete},
#'       state \code{NA}): \code{status}, \code{n_complete}, \code{state},
#'       \code{empirical_count/frequency}, \code{expected_mean/lower/upper}
#'       (= \code{model_mean/lower/upper}), \code{predictive_mean/lower/upper},
#'       \code{absolute_error}, \code{interval_contains_observed},
#'       \code{predictive_interval_contains_observed}}
#'     \item{\code{count_groups}}{one row per group x resistance-count
#'       \code{0..D} (COMPLETE events only -- total number of resistant
#'       classes): \code{status}, \code{n_complete}, \code{count},
#'       \code{empirical_count/frequency}, \code{expected_mean/lower/upper},
#'       \code{predictive_mean/lower/upper}, \code{absolute_error},
#'       \code{interval_contains_observed},
#'       \code{predictive_interval_contains_observed}}
#'     \item{\code{partial_patterns}}{one row per (group, exact test mask
#'       meeting \code{min_partial_pattern}, sub-profile over the TESTED
#'       classes in that mask) -- the denominator is EVERY event sharing
#'       that EXACT mask, not events co-tested as part of a larger panel
#'       (that is \code{pairwise}'s job) and not limited to one pair: an
#'       event tested on 3 of 5 panel classes contributes to exactly one
#'       mask's sub-profile breakdown here, at whatever order that mask
#'       is. Columns: \code{status}, \code{mask} (tested class names,
#'       "; "-joined), \code{n_mask}, \code{sub_profile} (R/S string over
#'       just the tested classes, in panel order), \code{empirical_count/
#'       frequency}, \code{expected_mean/lower/upper},
#'       \code{predictive_mean/lower/upper}, \code{absolute_error},
#'       \code{interval_contains_observed},
#'       \code{predictive_interval_contains_observed}. EVERY sub-profile of
#'       a listed mask is reported, including zero-count ones.}
#'     \item{\code{covariate_subgroups}}{empty (zero rows) unless
#'       \code{covariate_subgroups} is supplied; one row per (subgroup
#'       name, level, class), same columns as \code{marginal} but
#'       \code{group} replaced by \code{subgroup_name}/\code{subgroup_level}.}
#'   }
#' @export
profile_ppc <- function(fit_result, design, Y, ids,
                        n_draws = 500L,
                        min_tested = 10L, min_cotested = 10L, min_complete = 10L,
                        min_count_group = 10L, min_partial_pattern = 10L,
                        covariate_subgroups = NULL, min_subgroup = 10L,
                        ci_level = 0.95, seed = 123L) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for profile_ppc().", call. = FALSE)
  }
  fit <- if (is.list(fit_result) && "fit" %in% names(fit_result)) fit_result$fit else fit_result
  sd_list <- design$stan_data
  D <- sd_list$D
  K <- sd_list$K
  S <- sd_list$S
  P <- sd_list$P
  B <- sd_list$B
  V <- sd_list$V
  X <- sd_list$X
  N <- nrow(X)
  Y <- as.matrix(Y)
  classes <- colnames(Y)
  if (is.null(classes)) classes <- paste0("class", seq_len(D))
  groups <- unique(ids$group)
  state_labels <- apply(B, 1L, function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))

  draws_mat <- posterior::as_draws_matrix(fit$draws(variables = c("beta", "eta")))
  n_total <- nrow(draws_mat)
  .local_seed(as.integer(seed))
  n_draws <- min(as.integer(n_draws), n_total)
  draw_idx <- sort(sample.int(n_total, n_draws))
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  beta_cols <- sprintf("beta[%d,%d]", rep(seq_len(K), D), rep(seq_len(D), each = K))
  eta_cols <- if (P > 0L) sprintf("eta[%d]", seq_len(P)) else character(0)

  pair_idx <- if (D >= 2L) utils::combn(D, 2L, simplify = FALSE) else list()
  tested_idx_by_class <- lapply(seq_len(D), function(d) which(!is.na(Y[, d])))
  cotested_idx_by_pair <- lapply(pair_idx, function(p) which(!is.na(Y[, p[1]]) & !is.na(Y[, p[2]])))
  complete_idx <- which(rowSums(!is.na(Y)) == D)

  n_groups <- length(groups)
  marg_result <- array(NA_real_, dim = c(n_draws, n_groups, D))
  pair_result <- array(NA_real_, dim = c(n_draws, n_groups, length(pair_idx), 4L))
  profile_result <- array(NA_real_, dim = c(n_draws, n_groups, S))
  # Predictive (replicated-outcome) accumulators -- same shapes as above,
  # but built from one SAMPLED state per event per draw, not the expected
  # probability directly. See the function doc for why these are a
  # genuinely different (and necessarily wider) quantity.
  marg_pred <- array(NA_real_, dim = c(n_draws, n_groups, D))
  pair_pred <- array(NA_real_, dim = c(n_draws, n_groups, length(pair_idx), 4L))
  profile_pred <- array(NA_real_, dim = c(n_draws, n_groups, S))
  group_of_event <- ids$group
  tested_mask <- !is.na(Y) # N x D -- the ORIGINAL observation pattern, reapplied to every replicate
  popcount <- rowSums(B == 1) # length S -- resistance count (0..D) of each state

  # ---- count_groups precomputation: just needs complete_idx, already have it ----
  count_result <- array(NA_real_, dim = c(n_draws, n_groups, D + 1L))
  count_pred <- array(NA_real_, dim = c(n_draws, n_groups, D + 1L))

  # ---- partial_patterns precomputation ----
  # One combo per (group, EXACT test mask) with >= min_partial_pattern
  # events sharing that mask, excluding the all-tested (complete) and
  # all-untested masks (the latter cannot occur -- profile_prepare_model_data()
  # already excludes zero-observed events upstream). For each combo, every
  # sub-profile over just the tested classes gets a compatible-state index
  # set via profile_compatible_states() -- this is already exactly the
  # correct unconditional marginal-over-tested-classes computation (sum
  # Pmat over every state agreeing on the tested positions, regardless of
  # the untested ones), reused rather than re-derived.
  mask_key_of_event <- apply(tested_mask, 1L, function(r) paste(which(r), collapse = ","))
  partial_combos <- list()
  for (gi in seq_along(groups)) {
    g_idx <- which(group_of_event == groups[gi])
    g_idx <- setdiff(g_idx, complete_idx) # complete masks are a separate family
    if (length(g_idx) == 0L) next
    keys <- mask_key_of_event[g_idx]
    for (k in unique(keys)) {
      ev_idx <- g_idx[keys == k]
      if (length(ev_idx) < min_partial_pattern) next
      tested_classes <- as.integer(strsplit(k, ",")[[1]])
      kD <- length(tested_classes)
      sub_signs <- as.matrix(expand.grid(rep(list(c(-1L, 1L)), kD))[, rev(seq_len(kD)), drop = FALSE])
      compat_list <- lapply(seq_len(nrow(sub_signs)), function(r) {
        observed <- rep(0L, D)
        observed[tested_classes] <- sub_signs[r, ]
        profile_compatible_states(B, observed)
      })
      sub_labels <- apply(sub_signs, 1L, function(r) paste(ifelse(r == 1, "R", "S"), collapse = ""))
      partial_combos[[length(partial_combos) + 1]] <- list(
        group = groups[gi], tested_classes = tested_classes, ev_idx = ev_idx,
        compat_list = compat_list, sub_labels = sub_labels
      )
    }
  }
  n_combo <- length(partial_combos)
  partial_result <- lapply(partial_combos, function(cm) matrix(NA_real_, n_draws, length(cm$sub_labels)))
  partial_pred <- lapply(partial_combos, function(cm) matrix(NA_real_, n_draws, length(cm$sub_labels)))

  # ---- covariate_subgroups precomputation ----
  # Pooled across groups (these are prespecified clinical-subgroup checks,
  # not another per-hospital breakdown) -- one level entry per
  # (subgroup_name, level) pair actually present.
  subgroup_levels <- list()
  if (!is.null(covariate_subgroups)) {
    for (nm in names(covariate_subgroups)) {
      v <- covariate_subgroups[[nm]]
      if (length(v) != N) stop(sprintf("covariate_subgroups[['%s']] must have length nrow(Y) = %d", nm, N))
      for (lv in unique(stats::na.omit(v))) {
        subgroup_levels[[length(subgroup_levels) + 1]] <- list(
          name = nm, level = as.character(lv), idx = which(v == lv)
        )
      }
    }
  }
  n_sg <- length(subgroup_levels)
  sg_result <- array(NA_real_, dim = c(n_draws, n_sg, D))
  sg_pred <- array(NA_real_, dim = c(n_draws, n_sg, D))

  .local_seed(as.integer(seed) + 1L) # separate stream from the draw-subsampling seed above

  for (i in seq_len(n_draws)) {
    beta_s <- matrix(draws_mat[i, beta_cols], nrow = K, ncol = D)
    eta_s <- if (P > 0L) as.numeric(draws_mat[i, eta_cols]) else numeric(0)
    Mu <- X %*% beta_s # N x D
    Score <- Mu %*% t(B) # N x S
    v_eta <- if (P > 0L) as.vector(V %*% eta_s) else rep(0, S)
    Score <- sweep(Score, 2, v_eta, "+")
    Score <- Score - apply(Score, 1, max)
    Pmat <- exp(Score)
    Pmat <- Pmat / rowSums(Pmat) # N x S

    marg_prob <- vapply(seq_len(D), function(d) rowSums(Pmat[, B[, d] == 1, drop = FALSE]), numeric(N))

    # ---- ONE replicated complete profile per event, from THIS draw's Pmat ----
    # Vectorized categorical sampling: cumP[e, s] < u[e] recycles u (length N)
    # down cumP's columns correctly since nrow(cumP) == N, giving a per-row
    # (per-event) draw without an explicit loop.
    cumP <- t(apply(Pmat, 1, cumsum))
    u <- stats::runif(N)
    sim_state <- rowSums(cumP < u) + 1L
    sim_signs <- B[sim_state, , drop = FALSE] # N x D, each row the event's simulated full profile
    Y_rep <- matrix(NA_real_, N, D)
    Y_rep[tested_mask] <- ifelse(sim_signs[tested_mask] == 1, 1, 0) # reapply the ORIGINAL mask

    for (gi in seq_along(groups)) {
      g_idx <- which(group_of_event == groups[gi])
      for (d in seq_len(D)) {
        idx <- intersect(g_idx, tested_idx_by_class[[d]])
        if (length(idx) > 0L) {
          marg_result[i, gi, d] <- mean(marg_prob[idx, d])
          marg_pred[i, gi, d] <- mean(Y_rep[idx, d])
        }
      }
      for (pidx in seq_along(pair_idx)) {
        d1 <- pair_idx[[pidx]][1]
        d2 <- pair_idx[[pidx]][2]
        idx <- intersect(g_idx, cotested_idx_by_pair[[pidx]])
        if (length(idx) > 0L) {
          rr <- rowSums(Pmat[idx, B[, d1] == 1 & B[, d2] == 1, drop = FALSE])
          rs <- rowSums(Pmat[idx, B[, d1] == 1 & B[, d2] == -1, drop = FALSE])
          sr <- rowSums(Pmat[idx, B[, d1] == -1 & B[, d2] == 1, drop = FALSE])
          ss <- rowSums(Pmat[idx, B[, d1] == -1 & B[, d2] == -1, drop = FALSE])
          pair_result[i, gi, pidx, ] <- c(mean(rr), mean(rs), mean(sr), mean(ss))
          v1 <- Y_rep[idx, d1]
          v2 <- Y_rep[idx, d2]
          pair_pred[i, gi, pidx, ] <- c(
            mean(v1 == 1 & v2 == 1), mean(v1 == 1 & v2 == 0),
            mean(v1 == 0 & v2 == 1), mean(v1 == 0 & v2 == 0)
          )
        }
      }
      idx <- intersect(g_idx, complete_idx)
      if (length(idx) > 0L) {
        profile_result[i, gi, ] <- colMeans(Pmat[idx, , drop = FALSE])
        labels_rep <- apply(Y_rep[idx, , drop = FALSE], 1L, function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))
        tab_rep <- table(factor(labels_rep, levels = state_labels))
        profile_pred[i, gi, ] <- as.numeric(tab_rep) / length(idx)

        # ---- count_groups: total resistant classes, 0..D, COMPLETE events only ----
        p_by_event_by_count <- vapply(0:D, function(c) rowSums(Pmat[idx, popcount == c, drop = FALSE]), numeric(length(idx)))
        count_result[i, gi, ] <- colMeans(p_by_event_by_count)
        rep_count <- rowSums(Y_rep[idx, , drop = FALSE] == 1)
        count_pred[i, gi, ] <- vapply(0:D, function(c) mean(rep_count == c), numeric(1))
      }
    }

    # ---- partial_patterns: exact test-mask sub-profile breakdown ----
    for (ci in seq_len(n_combo)) {
      cm <- partial_combos[[ci]]
      p_vec <- vapply(cm$compat_list, function(compat) mean(rowSums(Pmat[cm$ev_idx, compat, drop = FALSE])), numeric(1))
      partial_result[[ci]][i, ] <- p_vec
      rep_labels <- apply(Y_rep[cm$ev_idx, cm$tested_classes, drop = FALSE], 1L,
                           function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))
      tab_rep <- table(factor(rep_labels, levels = cm$sub_labels))
      partial_pred[[ci]][i, ] <- as.numeric(tab_rep) / length(cm$ev_idx)
    }

    # ---- covariate_subgroups: same marginal check, stratified by subgroup level instead of group ----
    for (si in seq_len(n_sg)) {
      sg <- subgroup_levels[[si]]
      for (d in seq_len(D)) {
        idx <- intersect(sg$idx, tested_idx_by_class[[d]])
        if (length(idx) > 0L) {
          sg_result[i, si, d] <- mean(marg_prob[idx, d])
          sg_pred[i, si, d] <- mean(Y_rep[idx, d])
        }
      }
    }
  }

  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q

  # ---- marginal ----
  marginal_rows <- list()
  for (gi in seq_along(groups)) {
    g <- groups[gi]
    g_idx <- which(group_of_event == g)
    for (d in seq_len(D)) {
      idx <- intersect(g_idx, tested_idx_by_class[[d]])
      n_tested <- length(idx)
      if (n_tested < min_tested) {
        marginal_rows[[length(marginal_rows) + 1]] <- data.frame(
          group = g, class = classes[d], status = "not_computable",
          n_tested = n_tested, n_resistant = NA_integer_, n_susceptible = NA_integer_,
          observed_resistance = NA_real_,
          expected_mean = NA_real_, expected_lower = NA_real_, expected_upper = NA_real_,
          model_mean = NA_real_, model_lower = NA_real_, model_upper = NA_real_,
          predictive_mean = NA_real_, predictive_lower = NA_real_, predictive_upper = NA_real_,
          absolute_error = NA_real_, interval_contains_observed = NA,
          predictive_interval_contains_observed = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      n_resistant <- sum(Y[idx, d] == 1)
      n_susceptible <- n_tested - n_resistant
      observed_resistance <- n_resistant / n_tested
      model_draws <- marg_result[, gi, d]
      model_mean <- mean(model_draws)
      model_lower <- stats::quantile(model_draws, lo_q, names = FALSE)
      model_upper <- stats::quantile(model_draws, hi_q, names = FALSE)
      pred_draws <- marg_pred[, gi, d]
      pred_mean <- mean(pred_draws)
      pred_lower <- stats::quantile(pred_draws, lo_q, names = FALSE)
      pred_upper <- stats::quantile(pred_draws, hi_q, names = FALSE)
      marginal_rows[[length(marginal_rows) + 1]] <- data.frame(
        group = g, class = classes[d], status = "evaluated",
        n_tested = n_tested, n_resistant = n_resistant, n_susceptible = n_susceptible,
        observed_resistance = observed_resistance,
        expected_mean = model_mean, expected_lower = model_lower, expected_upper = model_upper,
        model_mean = model_mean, model_lower = model_lower, model_upper = model_upper,
        predictive_mean = pred_mean, predictive_lower = pred_lower, predictive_upper = pred_upper,
        absolute_error = abs(observed_resistance - model_mean),
        interval_contains_observed = observed_resistance >= model_lower & observed_resistance <= model_upper,
        predictive_interval_contains_observed = observed_resistance >= pred_lower & observed_resistance <= pred_upper,
        stringsAsFactors = FALSE
      )
    }
  }
  marginal_tbl <- do.call(rbind, marginal_rows)

  # ---- pairwise ----
  pairwise_rows <- list()
  for (gi in seq_along(groups)) {
    g <- groups[gi]
    g_idx <- which(group_of_event == g)
    for (pidx in seq_along(pair_idx)) {
      d1 <- pair_idx[[pidx]][1]
      d2 <- pair_idx[[pidx]][2]
      idx <- intersect(g_idx, cotested_idx_by_pair[[pidx]])
      n_cotested <- length(idx)
      if (n_cotested < min_cotested) {
        pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
          group = g, class_a = classes[d1], class_b = classes[d2], status = "not_computable",
          n_cotested = n_cotested, n_RR = NA_integer_, n_RS = NA_integer_,
          n_SR = NA_integer_, n_SS = NA_integer_,
          observed_RR_prop = NA_real_, observed_RS_prop = NA_real_,
          observed_SR_prop = NA_real_, observed_SS_prop = NA_real_,
          expected_RR_mean = NA_real_, expected_RR_lower = NA_real_, expected_RR_upper = NA_real_,
          expected_RS_mean = NA_real_, expected_SR_mean = NA_real_, expected_SS_mean = NA_real_,
          model_RR_mean = NA_real_, model_RR_lower = NA_real_, model_RR_upper = NA_real_,
          model_RS_mean = NA_real_, model_SR_mean = NA_real_, model_SS_mean = NA_real_,
          predictive_RR_mean = NA_real_, predictive_RR_lower = NA_real_, predictive_RR_upper = NA_real_,
          absolute_error = NA_real_, interval_contains_observed = NA,
          predictive_interval_contains_observed = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      v1 <- Y[idx, d1]
      v2 <- Y[idx, d2]
      n_RR <- sum(v1 == 1 & v2 == 1)
      n_RS <- sum(v1 == 1 & v2 == 0)
      n_SR <- sum(v1 == 0 & v2 == 1)
      n_SS <- sum(v1 == 0 & v2 == 0)
      rr_draws <- pair_result[, gi, pidx, 1]
      model_RR_mean <- mean(rr_draws)
      model_RR_lower <- stats::quantile(rr_draws, lo_q, names = FALSE)
      model_RR_upper <- stats::quantile(rr_draws, hi_q, names = FALSE)
      pred_rr_draws <- pair_pred[, gi, pidx, 1]
      pred_RR_mean <- mean(pred_rr_draws)
      pred_RR_lower <- stats::quantile(pred_rr_draws, lo_q, names = FALSE)
      pred_RR_upper <- stats::quantile(pred_rr_draws, hi_q, names = FALSE)
      observed_RR_prop <- n_RR / n_cotested
      pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
        group = g, class_a = classes[d1], class_b = classes[d2], status = "evaluated",
        n_cotested = n_cotested, n_RR = n_RR, n_RS = n_RS, n_SR = n_SR, n_SS = n_SS,
        observed_RR_prop = observed_RR_prop, observed_RS_prop = n_RS / n_cotested,
        observed_SR_prop = n_SR / n_cotested, observed_SS_prop = n_SS / n_cotested,
        expected_RR_mean = model_RR_mean, expected_RR_lower = model_RR_lower, expected_RR_upper = model_RR_upper,
        expected_RS_mean = mean(pair_result[, gi, pidx, 2]),
        expected_SR_mean = mean(pair_result[, gi, pidx, 3]),
        expected_SS_mean = mean(pair_result[, gi, pidx, 4]),
        model_RR_mean = model_RR_mean, model_RR_lower = model_RR_lower, model_RR_upper = model_RR_upper,
        model_RS_mean = mean(pair_result[, gi, pidx, 2]),
        model_SR_mean = mean(pair_result[, gi, pidx, 3]),
        model_SS_mean = mean(pair_result[, gi, pidx, 4]),
        predictive_RR_mean = pred_RR_mean, predictive_RR_lower = pred_RR_lower, predictive_RR_upper = pred_RR_upper,
        absolute_error = abs(observed_RR_prop - model_RR_mean),
        interval_contains_observed = observed_RR_prop >= model_RR_lower & observed_RR_prop <= model_RR_upper,
        predictive_interval_contains_observed = observed_RR_prop >= pred_RR_lower & observed_RR_prop <= pred_RR_upper,
        stringsAsFactors = FALSE
      )
    }
  }
  pairwise_tbl <- do.call(rbind, pairwise_rows)

  # ---- complete profile ----
  profile_rows <- list()
  for (gi in seq_along(groups)) {
    g <- groups[gi]
    g_idx <- which(group_of_event == g)
    idx <- intersect(g_idx, complete_idx)
    n_complete <- length(idx)
    if (n_complete < min_complete) {
      profile_rows[[length(profile_rows) + 1]] <- data.frame(
        group = g, status = "not_computable", n_complete = n_complete, state = NA_character_,
        empirical_count = NA_integer_, empirical_frequency = NA_real_,
        expected_mean = NA_real_, expected_lower = NA_real_, expected_upper = NA_real_,
        model_mean = NA_real_, model_lower = NA_real_, model_upper = NA_real_,
        predictive_mean = NA_real_, predictive_lower = NA_real_, predictive_upper = NA_real_,
        absolute_error = NA_real_, interval_contains_observed = NA,
        predictive_interval_contains_observed = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    labels <- apply(Y[idx, , drop = FALSE], 1L, function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))
    tab <- table(labels)
    for (s in seq_len(S)) {
      lbl <- state_labels[s]
      emp_count <- if (lbl %in% names(tab)) as.integer(tab[[lbl]]) else 0L
      emp_freq <- emp_count / n_complete
      model_draws <- profile_result[, gi, s]
      model_mean <- mean(model_draws)
      model_lower <- stats::quantile(model_draws, lo_q, names = FALSE)
      model_upper <- stats::quantile(model_draws, hi_q, names = FALSE)
      pred_draws <- profile_pred[, gi, s]
      pred_mean <- mean(pred_draws)
      pred_lower <- stats::quantile(pred_draws, lo_q, names = FALSE)
      pred_upper <- stats::quantile(pred_draws, hi_q, names = FALSE)
      profile_rows[[length(profile_rows) + 1]] <- data.frame(
        group = g, status = "evaluated", n_complete = n_complete, state = lbl,
        empirical_count = emp_count, empirical_frequency = emp_freq,
        expected_mean = model_mean, expected_lower = model_lower, expected_upper = model_upper,
        model_mean = model_mean, model_lower = model_lower, model_upper = model_upper,
        predictive_mean = pred_mean, predictive_lower = pred_lower, predictive_upper = pred_upper,
        absolute_error = abs(emp_freq - model_mean),
        interval_contains_observed = emp_freq >= model_lower & emp_freq <= model_upper,
        predictive_interval_contains_observed = emp_freq >= pred_lower & emp_freq <= pred_upper,
        stringsAsFactors = FALSE
      )
    }
  }
  profile_tbl <- do.call(rbind, profile_rows)

  # ---- count_groups ----
  count_rows <- list()
  for (gi in seq_along(groups)) {
    g <- groups[gi]
    g_idx <- which(group_of_event == g)
    idx <- intersect(g_idx, complete_idx)
    n_complete <- length(idx)
    if (n_complete < min_count_group) {
      for (c in 0:D) {
        count_rows[[length(count_rows) + 1]] <- data.frame(
          group = g, status = "not_computable", n_complete = n_complete, count = c,
          empirical_count = NA_integer_, empirical_frequency = NA_real_,
          expected_mean = NA_real_, expected_lower = NA_real_, expected_upper = NA_real_,
          predictive_mean = NA_real_, predictive_lower = NA_real_, predictive_upper = NA_real_,
          absolute_error = NA_real_, interval_contains_observed = NA,
          predictive_interval_contains_observed = NA, stringsAsFactors = FALSE
        )
      }
      next
    }
    obs_count <- rowSums(Y[idx, , drop = FALSE] == 1)
    for (c in 0:D) {
      emp_count <- sum(obs_count == c)
      emp_freq <- emp_count / n_complete
      model_draws <- count_result[, gi, c + 1L]
      model_mean <- mean(model_draws)
      model_lower <- stats::quantile(model_draws, lo_q, names = FALSE)
      model_upper <- stats::quantile(model_draws, hi_q, names = FALSE)
      pred_draws <- count_pred[, gi, c + 1L]
      pred_mean <- mean(pred_draws)
      pred_lower <- stats::quantile(pred_draws, lo_q, names = FALSE)
      pred_upper <- stats::quantile(pred_draws, hi_q, names = FALSE)
      count_rows[[length(count_rows) + 1]] <- data.frame(
        group = g, status = "evaluated", n_complete = n_complete, count = c,
        empirical_count = emp_count, empirical_frequency = emp_freq,
        expected_mean = model_mean, expected_lower = model_lower, expected_upper = model_upper,
        predictive_mean = pred_mean, predictive_lower = pred_lower, predictive_upper = pred_upper,
        absolute_error = abs(emp_freq - model_mean),
        interval_contains_observed = emp_freq >= model_lower & emp_freq <= model_upper,
        predictive_interval_contains_observed = emp_freq >= pred_lower & emp_freq <= pred_upper,
        stringsAsFactors = FALSE
      )
    }
  }
  count_tbl <- do.call(rbind, count_rows)

  # ---- partial_patterns ----
  partial_rows <- list()
  for (ci in seq_len(n_combo)) {
    cm <- partial_combos[[ci]]
    mask_label <- paste(classes[cm$tested_classes], collapse = "; ")
    n_mask <- length(cm$ev_idx)
    obs_labels <- apply(Y[cm$ev_idx, cm$tested_classes, drop = FALSE], 1L,
                         function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))
    tab_obs <- table(factor(obs_labels, levels = cm$sub_labels))
    for (si in seq_along(cm$sub_labels)) {
      emp_count <- as.integer(tab_obs[si])
      emp_freq <- emp_count / n_mask
      model_draws <- partial_result[[ci]][, si]
      model_mean <- mean(model_draws)
      model_lower <- stats::quantile(model_draws, lo_q, names = FALSE)
      model_upper <- stats::quantile(model_draws, hi_q, names = FALSE)
      pred_draws <- partial_pred[[ci]][, si]
      pred_mean <- mean(pred_draws)
      pred_lower <- stats::quantile(pred_draws, lo_q, names = FALSE)
      pred_upper <- stats::quantile(pred_draws, hi_q, names = FALSE)
      partial_rows[[length(partial_rows) + 1]] <- data.frame(
        group = cm$group, status = "evaluated", mask = mask_label, n_mask = n_mask,
        sub_profile = cm$sub_labels[si],
        empirical_count = emp_count, empirical_frequency = emp_freq,
        expected_mean = model_mean, expected_lower = model_lower, expected_upper = model_upper,
        predictive_mean = pred_mean, predictive_lower = pred_lower, predictive_upper = pred_upper,
        absolute_error = abs(emp_freq - model_mean),
        interval_contains_observed = emp_freq >= model_lower & emp_freq <= model_upper,
        predictive_interval_contains_observed = emp_freq >= pred_lower & emp_freq <= pred_upper,
        stringsAsFactors = FALSE
      )
    }
  }
  partial_tbl <- if (length(partial_rows) > 0L) {
    do.call(rbind, partial_rows)
  } else {
    data.frame(group = character(0), status = character(0), mask = character(0), n_mask = integer(0),
               sub_profile = character(0), empirical_count = integer(0), empirical_frequency = numeric(0),
               expected_mean = numeric(0), expected_lower = numeric(0), expected_upper = numeric(0),
               predictive_mean = numeric(0), predictive_lower = numeric(0), predictive_upper = numeric(0),
               absolute_error = numeric(0), interval_contains_observed = logical(0),
               predictive_interval_contains_observed = logical(0), stringsAsFactors = FALSE)
  }

  # ---- covariate_subgroups ----
  subgroup_rows <- list()
  for (si in seq_len(n_sg)) {
    sg <- subgroup_levels[[si]]
    for (d in seq_len(D)) {
      idx <- intersect(sg$idx, tested_idx_by_class[[d]])
      n_tested <- length(idx)
      if (n_tested < min_subgroup) {
        subgroup_rows[[length(subgroup_rows) + 1]] <- data.frame(
          subgroup_name = sg$name, subgroup_level = sg$level, class = classes[d], status = "not_computable",
          n_tested = n_tested, n_resistant = NA_integer_, n_susceptible = NA_integer_,
          observed_resistance = NA_real_,
          expected_mean = NA_real_, expected_lower = NA_real_, expected_upper = NA_real_,
          predictive_mean = NA_real_, predictive_lower = NA_real_, predictive_upper = NA_real_,
          absolute_error = NA_real_, interval_contains_observed = NA,
          predictive_interval_contains_observed = NA, stringsAsFactors = FALSE
        )
        next
      }
      n_resistant <- sum(Y[idx, d] == 1)
      observed_resistance <- n_resistant / n_tested
      model_draws <- sg_result[, si, d]
      model_mean <- mean(model_draws)
      model_lower <- stats::quantile(model_draws, lo_q, names = FALSE)
      model_upper <- stats::quantile(model_draws, hi_q, names = FALSE)
      pred_draws <- sg_pred[, si, d]
      pred_mean <- mean(pred_draws)
      pred_lower <- stats::quantile(pred_draws, lo_q, names = FALSE)
      pred_upper <- stats::quantile(pred_draws, hi_q, names = FALSE)
      subgroup_rows[[length(subgroup_rows) + 1]] <- data.frame(
        subgroup_name = sg$name, subgroup_level = sg$level, class = classes[d], status = "evaluated",
        n_tested = n_tested, n_resistant = n_resistant, n_susceptible = n_tested - n_resistant,
        observed_resistance = observed_resistance,
        expected_mean = model_mean, expected_lower = model_lower, expected_upper = model_upper,
        predictive_mean = pred_mean, predictive_lower = pred_lower, predictive_upper = pred_upper,
        absolute_error = abs(observed_resistance - model_mean),
        interval_contains_observed = observed_resistance >= model_lower & observed_resistance <= model_upper,
        predictive_interval_contains_observed = observed_resistance >= pred_lower & observed_resistance <= pred_upper,
        stringsAsFactors = FALSE
      )
    }
  }
  subgroup_tbl <- if (length(subgroup_rows) > 0L) {
    do.call(rbind, subgroup_rows)
  } else {
    data.frame(subgroup_name = character(0), subgroup_level = character(0), class = character(0),
               status = character(0), n_tested = integer(0), n_resistant = integer(0),
               n_susceptible = integer(0), observed_resistance = numeric(0),
               expected_mean = numeric(0), expected_lower = numeric(0), expected_upper = numeric(0),
               predictive_mean = numeric(0), predictive_lower = numeric(0), predictive_upper = numeric(0),
               absolute_error = numeric(0), interval_contains_observed = logical(0),
               predictive_interval_contains_observed = logical(0), stringsAsFactors = FALSE)
  }

  list(
    marginal = marginal_tbl, pairwise = pairwise_tbl, complete_profile = profile_tbl,
    count_groups = count_tbl, partial_patterns = partial_tbl, covariate_subgroups = subgroup_tbl
  )
}
