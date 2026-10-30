# profile_estimate.R
# Final per-event resistance-profile estimates from a profile_fit_joint_bernoulli()
# fit. See profile_fit.R's header for how this fits into the pipeline.
#
# This is the function the broader resistance-profiling pipeline design
# insisted could not infer its own answer: every downstream estimand
# (unconditional model prediction, observed-preserving conditional
# completion, or a discrete missing-data allocation) is a real, explicit
# argument, never a default this function picks for the caller.


#' Estimate Resistance Profiles From a Fitted Joint Bernoulli Model
#'
#' Produces per-event joint-profile estimates under one of three
#' EXPLICITLY DIFFERENT estimands -- \code{estimand} has no default, by
#' design:
#'
#' \describe{
#'   \item{\code{"unconditional"}}{The model's own implied profile
#'     distribution from an event's covariates alone, ignoring whatever
#'     was actually tested for that event. A population-level baseline,
#'     never a per-patient "best guess" when AST results exist for that
#'     patient.}
#'   \item{\code{"observed_preserving_conditional"}}{Never overwrites a
#'     tested cell: distributes probability mass only over the states
#'     compatible with the event's own observed R/S pattern (via
#'     \code{\link{profile_compatible_states}}) -- exactly the same exact
#'     marginalization the Stan model's own likelihood already uses for
#'     untested classes, not an approximation. A fully-tested event
#'     collapses to 100% on its one observed state; a fully-untested event
#'     reduces exactly to the unconditional distribution (every state is
#'     "compatible" when nothing was tested). This is the estimand
#'     DALY-style burden calculations need: never claim to know more about
#'     a patient's resistance pattern than was actually measured.}
#'   \item{\code{"missing_allocation"}}{Like
#'     \code{"observed_preserving_conditional"}, but instead of a
#'     continuous probability over compatible states, DRAWS one specific
#'     state per event per posterior draw from that same conditional
#'     distribution -- carrying allocation uncertainty forward as a set of
#'     discrete posterior draws, for downstream calculations that need an
#'     actual profile assignment per patient (e.g. counting how many
#'     patients fall in profile RRRR) with Monte Carlo uncertainty across
#'     draws, rather than fractional probability mass.}
#' }
#'
#' @param fit_result List from \code{\link{profile_fit_joint_bernoulli}}
#'   (reads \code{$fit}), or a raw \code{cmdstanr} \code{CmdStanMCMC}.
#' @param design List from \code{\link{profile_build_design}} (reads
#'   \code{$stan_data} and \code{$design$states}).
#' @param Y,ids The SAME \code{Y} matrix and \code{ids} data frame passed
#'   to \code{\link{profile_build_design}}.
#' @param estimand Character, one of \code{"unconditional"},
#'   \code{"observed_preserving_conditional"}, \code{"missing_allocation"}.
#'   No default -- the caller must choose.
#' @param n_draws Integer. Posterior draws used (subsampled if the fit has
#'   more). Default \code{500L}. Note: for \code{"unconditional"}/
#'   \code{"observed_preserving_conditional"}, the returned
#'   \code{profile_probability_draws} array is \code{n_draws x
#'   nrow(Y) x 2^D} -- memory scales with all three, so reduce
#'   \code{n_draws} for a very large \code{Y} with a large panel.
#' @param ci_level Numeric in (0, 1). Credible-interval level for
#'   \code{profile_probability_lower}/\code{upper} and
#'   \code{class_probability_lower}/\code{upper}. Default \code{0.95}
#'   (giving the \code{[Q_.025, Q_.975]} marginal CrI). Ignored for
#'   \code{estimand = "missing_allocation"} (its uncertainty is already
#'   the full set of discrete per-draw allocations, not a further
#'   interval). These are MARGINAL per-event/per-class/per-state CrIs,
#'   each computed by aggregating WITHIN a draw first (summing state
#'   probabilities into a class margin, say) and taking quantiles ACROSS
#'   draws second -- never the reverse, and never a simultaneous region
#'   across multiple events/classes/states at once.
#' @param seed Integer. Default \code{123L}.
#'
#' @return A list with \code{estimand} (echoes the argument, for
#'   provenance) and \code{state_labels} (character vector, length
#'   \code{2^D}, e.g. \code{"RRSS"}), plus:
#'   \describe{
#'     \item{\code{"unconditional"} / \code{"observed_preserving_conditional"}}{
#'       \code{profile_probability}: matrix, events x \code{2^D} states,
#'       the posterior-mean probability of each state per event;
#'       \code{profile_probability_lower}/\code{profile_probability_upper}:
#'       same shape, the \code{ci_level} marginal CrI per event x state;
#'       \code{profile_probability_draws}: the full \code{n_draws x
#'       events x 2^D} array; \code{class_probability}/
#'       \code{class_probability_lower}/\code{class_probability_upper}:
#'       matrix, events x \code{D} classes -- the per-class resistance
#'       margin (states summed WITHIN each draw, then meaned/quantiled
#'       across draws), with its own CrI, not derived by summing the
#'       per-state CrIs (that would NOT be a valid interval for the sum));
#'       \code{class_probability_draws}: the full \code{n_draws x events x
#'       D} array behind \code{class_probability} (e.g. for feeding
#'       \code{\link{profile_check_sampler}}'s \code{derived_estimates}
#'       MCSE check).}
#'     \item{\code{"missing_allocation"}}{\code{allocated_state}: integer
#'       matrix, \code{n_draws x events}, each entry a state index into
#'       \code{state_labels}.}
#'   }
#' @export
profile_estimate <- function(fit_result, design, Y, ids,
                             estimand,
                             n_draws = 500L,
                             ci_level = 0.95,
                             seed = 123L) {
  estimand <- match.arg(estimand, c("unconditional", "observed_preserving_conditional", "missing_allocation"))
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for profile_estimate().", call. = FALSE)
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
  N <- nrow(X)
  Y <- as.matrix(Y)
  state_labels <- apply(B, 1L, function(row) paste(ifelse(row == 1, "R", "S"), collapse = ""))

  draws_mat <- posterior::as_draws_matrix(fit$draws(variables = c("beta", "eta")))
  n_total <- nrow(draws_mat)
  .local_seed(as.integer(seed))
  n_draws <- min(as.integer(n_draws), n_total)
  draw_idx <- sort(sample.int(n_total, n_draws))
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  beta_cols <- sprintf("beta[%d,%d]", rep(seq_len(K), D), rep(seq_len(D), each = K))
  eta_cols <- if (P > 0L) sprintf("eta[%d]", seq_len(P)) else character(0)

  needs_conditioning <- estimand != "unconditional"
  sign_mat <- matrix(0L, N, D)
  sign_mat[Y == 1] <- 1L
  sign_mat[Y == 0] <- -1L
  # Compatible-state sets depend only on the event's own observed pattern,
  # never on the posterior draw -- computed once here, reused across every
  # draw below, instead of re-deriving it N * n_draws times.
  compat_list <- if (needs_conditioning) {
    lapply(seq_len(N), function(e) profile_compatible_states(states, sign_mat[e, ]))
  } else {
    NULL
  }

  if (estimand == "missing_allocation") {
    allocated_state <- matrix(NA_integer_, n_draws, N)
  } else {
    profile_probability_draws <- array(NA_real_, dim = c(n_draws, N, S))
  }

  for (i in seq_len(n_draws)) {
    beta_s <- matrix(draws_mat[i, beta_cols], nrow = K, ncol = D)
    eta_s <- if (P > 0L) as.numeric(draws_mat[i, eta_cols]) else numeric(0)
    Mu <- X %*% beta_s # N x D
    Score <- Mu %*% t(B) # N x S
    v_eta <- if (P > 0L) as.vector(V %*% eta_s) else rep(0, S)
    Score <- sweep(Score, 2, v_eta, "+")
    Score <- Score - apply(Score, 1, max)
    Pmat <- exp(Score)
    Pmat <- Pmat / rowSums(Pmat) # N x S, unconditional per event

    if (estimand == "unconditional") {
      profile_probability_draws[i, , ] <- Pmat
      next
    }

    Pcond <- matrix(0, N, S)
    for (e in seq_len(N)) {
      idx <- compat_list[[e]]
      denom <- sum(Pmat[e, idx])
      if (denom > 0) Pcond[e, idx] <- Pmat[e, idx] / denom
    }

    if (estimand == "observed_preserving_conditional") {
      profile_probability_draws[i, , ] <- Pcond
    } else {
      for (e in seq_len(N)) {
        allocated_state[i, e] <- sample.int(S, 1L, prob = Pcond[e, ])
      }
    }
  }

  if (estimand == "missing_allocation") {
    list(estimand = estimand, state_labels = state_labels, allocated_state = allocated_state)
  } else {
    profile_probability <- apply(profile_probability_draws, c(2, 3), mean)
    dimnames(profile_probability) <- list(NULL, state_labels)

    lo_q <- (1 - ci_level) / 2
    hi_q <- 1 - lo_q
    profile_probability_lower <- apply(profile_probability_draws, c(2, 3), stats::quantile, probs = lo_q, names = FALSE)
    profile_probability_upper <- apply(profile_probability_draws, c(2, 3), stats::quantile, probs = hi_q, names = FALSE)
    dimnames(profile_probability_lower) <- dimnames(profile_probability_upper) <- list(NULL, state_labels)

    # Class margins: states summed WITHIN each draw first (a per-draw
    # sum is still a valid probability), THEN mean/quantile across draws
    # -- never quantiles-of-states summed afterward, which would not be a
    # correct interval for the sum.
    class_probability_draws <- array(NA_real_, dim = c(n_draws, N, D))
    for (d in seq_len(D)) {
      class_probability_draws[, , d] <- apply(profile_probability_draws[, , B[, d] == 1, drop = FALSE], c(1, 2), sum)
    }
    class_probability <- apply(class_probability_draws, c(2, 3), mean)
    class_probability_lower <- apply(class_probability_draws, c(2, 3), stats::quantile, probs = lo_q, names = FALSE)
    class_probability_upper <- apply(class_probability_draws, c(2, 3), stats::quantile, probs = hi_q, names = FALSE)
    class_names <- colnames(Y)
    if (is.null(class_names)) class_names <- paste0("class", seq_len(D))
    dimnames(class_probability) <- dimnames(class_probability_lower) <- dimnames(class_probability_upper) <-
      list(NULL, class_names)

    list(
      estimand = estimand, state_labels = state_labels,
      profile_probability = profile_probability,
      profile_probability_lower = profile_probability_lower,
      profile_probability_upper = profile_probability_upper,
      profile_probability_draws = profile_probability_draws,
      class_probability = class_probability,
      class_probability_lower = class_probability_lower,
      class_probability_upper = class_probability_upper,
      class_probability_draws = class_probability_draws
    )
  }
}
