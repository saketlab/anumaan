# profile_aggregate_counts_allocation.R
# Sibling to profile_aggregate_daly.R, for a genuinely different question.
#
# profile_aggregate_for_daly() answers: what PROPORTION of a group is
# EXPECTED to have profile s, as a continuous quantity (mean over draws
# of a per-draw mean probability)?
#
# This file answers: what is the distribution of the REALISED, INTEGER
# completed count N_{gs}, once every genuinely missing profile has
# actually been resolved to one specific state (not left as fractional
# probability mass)? That needs profile_estimate()'s discrete
# "missing_allocation" estimand, whose $allocated_state (n_draws x N
# integer matrix of state indices) is a different shape from
# profile_aggregate_for_daly()'s expected $profile_probability_draws
# (continuous n_draws x N x S array) -- hence a new function rather than
# a patch to that one's existing, documented contract.


#' Aggregate Hospital-Profile Counts Under Discrete Missing-Data Allocation
#'
#' A complete event's compatible-state set is always a single state (see
#' \code{\link{profile_estimate}}'s \code{"missing_allocation"} doc), so
#' it deterministically allocates to its own true state on every draw --
#' there is no need to separately add an "n_observed" term. Passing
#' complete AND partial events together (the same convention
#' \code{\link{profile_aggregate_for_daly}} already uses) reproduces the
#' observed count automatically, with the partial events' allocation
#' uncertainty layered on top.
#'
#' For hospital/group \eqn{h} and profile \eqn{s}, at each posterior draw
#' \eqn{m}:
#' \deqn{N_{hs}^{(m)} = \sum_{e \in h} I(\mathrm{allocated\_state}[m, e] = s)}
#' Collecting \eqn{N_{hs}^{(1)}, \ldots, N_{hs}^{(M)}} gives the realised
#' completed-count distribution for that hospital/profile; this function
#' reports its mean and \code{ci_level} quantile interval. This single
#' procedure automatically mixes parameter uncertainty (the draw-to-draw
#' allocation weights shift because \eqn{\theta^{(m)}} shifts) and
#' allocation uncertainty (even at fixed weights, which specific state
#' gets drawn is itself random) -- the two-term variance decomposition
#' discussed alongside this function's design, obtained here by direct
#' simulation rather than a closed-form formula.
#'
#' @param est List from \code{\link{profile_estimate}} with
#'   \code{estimand = "missing_allocation"} (reads \code{$allocated_state},
#'   \code{$state_labels}).
#' @param ids Data frame with a \code{group} column, same row order/length
#'   as the \code{Y} passed to \code{\link{profile_estimate}} (one row per
#'   event).
#' @param pathogen,cohort_definition Character, echoed into the output for
#'   provenance (same convention as \code{\link{profile_aggregate_for_daly}}).
#' @param sampler_check List or \code{NULL}. Same convention as
#'   \code{\link{profile_aggregate_for_daly}}: reads \code{$status};
#'   \code{"pass"} is the only acceptable value, anything else (including
#'   \code{NULL}) makes every row ineligible with a named reason.
#' @param ci_level Numeric in (0, 1). Quantile interval level for the
#'   realised-count distribution. Default \code{0.95}.
#' @param min_n_events,min_n_draws Integer floors below which a group x
#'   profile row is marked ineligible (same meaning as
#'   \code{\link{profile_aggregate_for_daly}}). Defaults \code{10L},
#'   \code{100L}.
#'
#' @return Data frame, one row per group x profile: \code{N_mean},
#'   \code{N_lower}, \code{N_upper} (the realised completed-count
#'   distribution's mean and \code{ci_level} quantile interval --
#'   integer-valued draws summarised, not a continuous proportion),
#'   \code{n_events}, \code{n_draws}, \code{n_profile_classes},
#'   \code{sampler_status}, \code{computational_checks_passed_YLL}/
#'   \code{computational_checks_passed_YLD} and their
#'   \code{exclusion_reason_*} (same computational-check logic and column
#'   names as \code{\link{profile_aggregate_for_daly}}, for a consistent
#'   downstream filter regardless of which aggregator produced a row --
#'   and see that function's doc for what these flags do NOT establish:
#'   PPC adequacy, held-out validation, missingness sensitivity, death/
#'   episode attribution, double counting, or RR/duration/disability
#'   input availability. This is a computational pre-check, not a
#'   burden-use release decision.), and
#'   \code{estimator = "joint_bernoulli_discrete_allocation"} so the two
#'   aggregators' outputs are never silently confused if concatenated.
#' @export
profile_aggregate_counts_with_allocation <- function(est, ids, pathogen, cohort_definition,
                                                      sampler_check = NULL,
                                                      ci_level = 0.95,
                                                      min_n_events = 10L,
                                                      min_n_draws = 100L) {
  if (!is.list(est) || !all(c("allocated_state", "state_labels") %in% names(est))) {
    stop(
      "est must be the list returned by profile_estimate() with estimand = 'missing_allocation' ",
      "(needs $allocated_state, $state_labels).",
      call. = FALSE
    )
  }
  if (!identical(est$estimand, "missing_allocation")) {
    stop(sprintf(
      "[profile_aggregate_counts_with_allocation] est$estimand is '%s', not 'missing_allocation' -- this function is specifically for the discrete-allocation estimand; use profile_aggregate_for_daly() for continuous probabilities.",
      est$estimand
    ), call. = FALSE)
  }

  allocated <- est$allocated_state # n_draws x N, integer state indices
  state_labels <- est$state_labels
  n_draws_total <- nrow(allocated)
  S <- length(state_labels)
  groups <- unique(ids$group)

  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q

  sampler_status <- if (!is.null(sampler_check)) sampler_check$status else "not_checked"
  sampler_acceptable <- identical(sampler_status, "pass")

  rows <- list()
  for (g in groups) {
    g_idx <- which(ids$group == g)
    n_events_g <- length(g_idx)
    for (s in seq_len(S)) {
      N_draws <- if (n_events_g > 0L) {
        vapply(seq_len(n_draws_total), function(i) sum(allocated[i, g_idx] == s), integer(1))
      } else {
        rep(NA_integer_, n_draws_total)
      }
      n_draws_valid <- sum(!is.na(N_draws))
      N_mean <- mean(N_draws, na.rm = TRUE)
      N_lower <- if (n_draws_valid > 0L) stats::quantile(N_draws, lo_q, na.rm = TRUE, names = FALSE) else NA_real_
      N_upper <- if (n_draws_valid > 0L) stats::quantile(N_draws, hi_q, na.rm = TRUE, names = FALSE) else NA_real_
      if (is.nan(N_mean)) N_mean <- NA_real_

      low_events <- n_events_g < min_n_events
      low_draws <- n_draws_valid < min_n_draws
      eligible <- sampler_acceptable && !low_events && !low_draws

      reason <- if (eligible) {
        NA_character_
      } else if (is.null(sampler_check)) {
        "sampler_not_checked"
      } else if (!sampler_acceptable) {
        "sampler_not_acceptable"
      } else if (n_events_g == 0L) {
        "no_events"
      } else if (low_events) {
        "too_few_events"
      } else {
        "too_few_valid_draws"
      }

      rows[[length(rows) + 1]] <- data.frame(
        group = g, pathogen = pathogen, cohort_definition = cohort_definition,
        profile_delta = state_labels[s],
        N_mean = round(N_mean, 3), N_lower = N_lower, N_upper = N_upper,
        n_events = n_events_g, n_draws = n_draws_valid, n_profile_classes = round(log2(S)),
        sampler_status = sampler_status,
        computational_checks_passed_YLL = eligible, exclusion_reason_YLL = reason,
        computational_checks_passed_YLD = eligible, exclusion_reason_YLD = reason,
        estimator = "joint_bernoulli_discrete_allocation",
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, rows)
  message(sprintf(
    "[profile_aggregate_counts_with_allocation] %d group x profile rows (%d pass computational checks for YLL/YLD -- NOT a burden-use release decision).",
    nrow(result), sum(result$computational_checks_passed_YLL)
  ))
  result
}
