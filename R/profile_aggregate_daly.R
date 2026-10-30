# profile_aggregate_daly.R
# Bridges profile_estimate()'s per-event, per-draw output into the
# hospital x pathogen x profile summary DALY calculation needs. See
# profile_fit.R's header for how this fits into the pipeline.


#' Aggregate Posterior Profile Estimates for DALY Use
#'
#' Bridges \code{\link{profile_estimate}}'s per-event, per-draw
#' state-probability array into the hospital x pathogen x profile summary
#' DALY calculation needs: posterior mean + credible interval, panel
#' composition, event/draw support counts, and explicit
#' \code{eligible_for_YLL}/\code{eligible_for_YLD} flags with reasons.
#'
#' This does NOT split into \code{R_ALL}/\code{R_KNOWN_OUTCOME}/\code{R_NF}
#' cohorts after fitting. In this pipeline, cohort selection (e.g. a Discharged/
#' non-fatal restriction) happens BEFORE fitting, via
#' \code{\link{profile_select_cohort}} -- the fit is already trained on
#' exactly the cohort a DALY estimand wants, never on a broader population
#' split apart afterward. Reporting a post-hoc "R_ALL" alongside an
#' "R_NF" from a model that was only ever fit on non-fatal events would
#' imply a second population the model was never trained on and has no
#' actual information about. Instead, this reports ONE
#' \code{R_mean}/\code{R_lower}/\code{R_upper} per group x profile, for
#' whatever cohort was actually fit, carried in the
#' \code{cohort_definition} column so the number is never mistaken for a
#' different population's. A caller who genuinely needs both an
#' all-events and a non-fatal-only DALY estimand gets both by running
#' \code{\link{profile_select_cohort}} twice and fitting twice -- two real
#' fits, not one fit relabeled twice.
#'
#' A second departure: \code{sampler_acceptable} has no default here. The
#' legacy function defaulted it to \code{TRUE} (opt-out: eligible unless
#' told otherwise). This function requires an actual
#' \code{\link{profile_check_sampler}} result and treats a missing one as
#' not eligible, with \code{exclusion_reason_YLL/YLD = "sampler_not_checked"}
#' -- never silently assumed acceptable.
#'
#' @param est List from \code{\link{profile_estimate}}, called with
#'   \code{estimand = "observed_preserving_conditional"} (the DALY-
#'   appropriate estimand -- never overwrites a tested cell). A warning is
#'   issued if \code{est$estimand} is anything else.
#' @param ids Data frame, the SAME \code{ids} passed to
#'   \code{\link{profile_estimate}} (\code{event_id}, \code{group}).
#' @param pathogen Character, length 1. \code{\link{profile_estimate}}
#'   always operates on one pathogen's panel at a time; this just labels
#'   the output.
#' @param cohort_definition Character, length 1. Free-text description of
#'   which events were fit (e.g. \code{"Discharged only"}, \code{"all
#'   events, no outcome restriction"}) -- echoed into every row.
#' @param sampler_check List or \code{NULL}. Return value of
#'   \code{\link{profile_check_sampler}}. \code{sampler_check$status ==
#'   "pass"} is required for \code{eligible_for_YLL}/\code{eligible_for_YLD}.
#'   \code{NULL} (default): both are always \code{FALSE} with reason
#'   \code{"sampler_not_checked"}.
#' @param ci_level Numeric. Default \code{0.95}.
#' @param min_n_events,min_n_draws Integer. Eligibility floors, same
#'   meaning as the legacy function's. Default \code{10L}/\code{100L}.
#'
#' @return A data frame, one row per group x profile (state):
#'   \code{group}, \code{pathogen}, \code{cohort_definition},
#'   \code{profile_delta} (state label, e.g. \code{"RRSS"}),
#'   \code{R_mean}, \code{R_lower}, \code{R_upper}, \code{n_events},
#'   \code{n_draws}, \code{n_profile_classes}, \code{sampler_status},
#'   \code{eligible_for_YLL}, \code{exclusion_reason_YLL},
#'   \code{eligible_for_YLD}, \code{exclusion_reason_YLD} (identical to
#'   the YLL columns here -- there is only one cohort to be eligible
#'   for), \code{estimator} (\code{"joint_bernoulli"}).
#' @export
profile_aggregate_for_daly <- function(est, ids, pathogen, cohort_definition,
                                       sampler_check = NULL,
                                       ci_level = 0.95,
                                       min_n_events = 10L,
                                       min_n_draws = 100L) {
  if (!is.list(est) || !all(c("profile_probability_draws", "state_labels") %in% names(est))) {
    stop(
      "est must be the list returned by profile_estimate() with estimand != 'missing_allocation' ",
      "(needs $profile_probability_draws, $state_labels)."
    )
  }
  if (!identical(est$estimand, "observed_preserving_conditional")) {
    warning(sprintf(
      "[profile_aggregate_for_daly] est$estimand is '%s', not 'observed_preserving_conditional' -- ",
      est$estimand
    ), "DALY use should never overwrite a tested cell with an unconditional model prediction.",
    call. = FALSE)
  }

  draws <- est$profile_probability_draws # n_draws x N x S
  state_labels <- est$state_labels
  n_draws_total <- dim(draws)[1]
  S <- dim(draws)[3]
  D <- round(log2(S))
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
      draw_means <- if (n_events_g > 0L) {
        vapply(seq_len(n_draws_total), function(i) mean(draws[i, g_idx, s]), numeric(1))
      } else {
        rep(NA_real_, n_draws_total)
      }
      n_draws_valid <- sum(!is.na(draw_means))
      R_mean <- mean(draw_means, na.rm = TRUE)
      R_lower <- if (n_draws_valid > 0L) stats::quantile(draw_means, lo_q, na.rm = TRUE, names = FALSE) else NA_real_
      R_upper <- if (n_draws_valid > 0L) stats::quantile(draw_means, hi_q, na.rm = TRUE, names = FALSE) else NA_real_
      if (is.nan(R_mean)) R_mean <- NA_real_

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
        R_mean = round(R_mean, 6), R_lower = round(R_lower, 6), R_upper = round(R_upper, 6),
        n_events = n_events_g, n_draws = n_draws_valid, n_profile_classes = D,
        sampler_status = sampler_status,
        eligible_for_YLL = eligible, exclusion_reason_YLL = reason,
        eligible_for_YLD = eligible, exclusion_reason_YLD = reason,
        estimator = "joint_bernoulli",
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, rows)
  message(sprintf(
    "[profile_aggregate_for_daly] %d group x profile rows (%d eligible for YLL/YLD).",
    nrow(result), sum(result$eligible_for_YLL)
  ))
  result
}
