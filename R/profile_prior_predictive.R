# profile_prior_predictive.R
# Prior predictive checks for a profile_build_design() design, BEFORE any
# fitting. See profile_fit.R's header for how this fits into the
# pipeline. Uses the EXACT SAME prediction kernel profile_ppc()/
# profile_estimate() use for the posterior (softmax-over-states, B/V
# design) -- just with parameters drawn from the declared priors
# (beta_sd/eta_sd, already recorded in design$stan_data) instead of from
# posterior draws. No Stan run, no cmdstanr dependency, no fit_result.


#' Prior Predictive Checks for a Joint Bernoulli Design
#'
#' Draws \code{beta}/\code{eta} from the declared priors
#' (\code{design$stan_data$beta_sd}/\code{eta_sd} -- i.e. whatever
#' \code{\link{profile_build_design}} was actually called with, not a
#' separately-specified prior), computes the implied joint-state
#' distribution at each event's REAL covariates via the SAME prediction
#' kernel the Stan model/\code{\link{profile_ppc}}/
#' \code{\link{profile_estimate}} use, and summarizes it by group:
#' per-class resistance margins, selected pair probabilities, the two most
#' extreme profiles (all-resistant, all-susceptible), the resistance-count
#' distribution, and a concentration/entropy summary. This is a check on
#' the PRIOR alone -- it never looks at \code{Y}, and a prior that looks
#' implausible here is a prior problem, not evidence the data disagrees
#' with it.
#'
#' Inspect whether the prior puts most of its mass on near-deterministic
#' profiles (entropy near 0) or implausibly strong cross-class
#' associations (e.g. all-resistant/all-susceptible probability far above
#' what independence at the marginal rates would give) across most
#' hospitals BEFORE fitting. When increasing \code{D} or the interaction
#' order, re-run this: many individually modest coefficients can combine
#' into an extreme combined score, so a prior that looked fine at a
#' smaller \code{D} is not guaranteed to still look fine at a larger one
#' -- rescaling by a fixed formula without re-inspecting the induced
#' distribution is exactly the mistake this function exists to catch.
#'
#' @param design List from \code{\link{profile_build_design}} (reads
#'   \code{$stan_data}: \code{X}, \code{B}, \code{V}, \code{beta_sd},
#'   \code{eta_sd}; and \code{$design}: \code{class_names}).
#' @param ids Data frame, one row per event (matching \code{nrow(X)}),
#'   with a \code{group} column.
#' @param n_draws Integer. Prior draws. Default \code{1000L} (prior draws
#'   are cheap -- no Stan sampling -- so this can be large).
#' @param ci_level Numeric. Default \code{0.95}.
#' @param seed Integer. Default \code{123L}.
#'
#' @return A list of data frames, each one row per group (x the family's
#'   own extra key):
#'   \describe{
#'     \item{\code{class_margins}}{group x class: \code{mean},
#'       \code{lower}, \code{upper} prior-predictive resistance
#'       probability.}
#'     \item{\code{pair_probabilities}}{group x class-pair x cell
#'       (\code{RR}/\code{RS}/\code{SR}/\code{SS}): \code{mean},
#'       \code{lower}, \code{upper}.}
#'     \item{\code{extreme_profiles}}{group: \code{all_resistant_mean/
#'       lower/upper}, \code{all_susceptible_mean/lower/upper}, and
#'       \code{independence_all_resistant}/\code{independence_all_susceptible}
#'       -- what the all-R/all-S probability WOULD be under independence
#'       at this group's OWN mean marginal rates (the product of
#'       \code{class_margins}), so a prior-induced value far above the
#'       independence baseline signals strong prior-induced association,
#'       not just "a small number".}
#'     \item{\code{count_distribution}}{group x \code{count} (\code{0..D}
#'       resistant classes): \code{mean}, \code{lower}, \code{upper}.}
#'     \item{\code{concentration}}{group: \code{mean_entropy} (mean,
#'       across events and draws, of the per-event Shannon entropy of the
#'       \code{2^D}-state distribution, in nats), \code{max_entropy}
#'       (\eqn{\log(2^D)}, the uniform-distribution ceiling), and
#'       \code{entropy_ratio} (\code{mean_entropy / max_entropy} -- near 0
#'       means near-deterministic profiles; near 1 means close to the
#'       uniform/no-information prior).}
#'   }
#' @export
profile_prior_predictive <- function(design, ids, n_draws = 1000L, ci_level = 0.95, seed = 123L) {
  sd_list <- design$stan_data
  D <- sd_list$D
  K <- sd_list$K
  P <- sd_list$P
  S <- sd_list$S
  B <- sd_list$B
  V <- sd_list$V
  X <- sd_list$X
  beta_sd <- sd_list$beta_sd
  eta_sd <- sd_list$eta_sd
  N <- nrow(X)
  classes <- design$design$class_names
  if (is.null(classes)) classes <- paste0("class", seq_len(D))
  groups <- unique(ids$group)
  n_groups <- length(groups)
  group_of_event <- ids$group
  popcount <- rowSums(B == 1)

  pair_idx <- if (D >= 2L) utils::combn(D, 2L, simplify = FALSE) else list()

  .local_seed(as.integer(seed))

  class_margin_draws <- array(NA_real_, dim = c(n_draws, n_groups, D))
  pair_draws <- array(NA_real_, dim = c(n_draws, n_groups, length(pair_idx), 4L))
  allR_draws <- array(NA_real_, dim = c(n_draws, n_groups))
  allS_draws <- array(NA_real_, dim = c(n_draws, n_groups))
  count_draws <- array(NA_real_, dim = c(n_draws, n_groups, D + 1L))
  entropy_draws <- array(NA_real_, dim = c(n_draws, n_groups))

  all_R_state <- which(popcount == D)
  all_S_state <- which(popcount == 0L)

  for (i in seq_len(n_draws)) {
    beta_s <- matrix(stats::rnorm(K * D, mean = 0, sd = beta_sd), nrow = K, ncol = D)
    eta_s <- if (P > 0L) stats::rnorm(P, mean = 0, sd = eta_sd) else numeric(0)
    Mu <- X %*% beta_s # N x D
    Score <- Mu %*% t(B) # N x S
    v_eta <- if (P > 0L) as.vector(V %*% eta_s) else rep(0, S)
    Score <- sweep(Score, 2, v_eta, "+")
    Score <- Score - apply(Score, 1, max)
    Pmat <- exp(Score)
    Pmat <- Pmat / rowSums(Pmat) # N x S

    marg_prob <- vapply(seq_len(D), function(d) rowSums(Pmat[, B[, d] == 1, drop = FALSE]), numeric(N))
    count_prob <- vapply(0:D, function(c) rowSums(Pmat[, popcount == c, drop = FALSE]), numeric(N))
    entropy_e <- -rowSums(Pmat * log(pmax(Pmat, 1e-300)))

    for (gi in seq_along(groups)) {
      g_idx <- which(group_of_event == groups[gi])
      class_margin_draws[i, gi, ] <- colMeans(marg_prob[g_idx, , drop = FALSE])
      allR_draws[i, gi] <- mean(Pmat[g_idx, all_R_state])
      allS_draws[i, gi] <- mean(Pmat[g_idx, all_S_state])
      count_draws[i, gi, ] <- colMeans(count_prob[g_idx, , drop = FALSE])
      entropy_draws[i, gi] <- mean(entropy_e[g_idx])
      for (pidx in seq_along(pair_idx)) {
        d1 <- pair_idx[[pidx]][1]
        d2 <- pair_idx[[pidx]][2]
        rr <- mean(rowSums(Pmat[g_idx, B[, d1] == 1 & B[, d2] == 1, drop = FALSE]))
        rs <- mean(rowSums(Pmat[g_idx, B[, d1] == 1 & B[, d2] == -1, drop = FALSE]))
        sr <- mean(rowSums(Pmat[g_idx, B[, d1] == -1 & B[, d2] == 1, drop = FALSE]))
        ss <- mean(rowSums(Pmat[g_idx, B[, d1] == -1 & B[, d2] == -1, drop = FALSE]))
        pair_draws[i, gi, pidx, ] <- c(rr, rs, sr, ss)
      }
    }
  }

  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q
  .summ <- function(x) c(mean = mean(x), lower = stats::quantile(x, lo_q, names = FALSE), upper = stats::quantile(x, hi_q, names = FALSE))

  class_margins <- do.call(rbind, lapply(seq_along(groups), function(gi) {
    do.call(rbind, lapply(seq_len(D), function(d) {
      s <- .summ(class_margin_draws[, gi, d])
      data.frame(group = groups[gi], class = classes[d], mean = s[["mean"]], lower = s[["lower"]], upper = s[["upper"]], stringsAsFactors = FALSE)
    }))
  }))

  pair_probabilities <- if (length(pair_idx) > 0L) {
    do.call(rbind, lapply(seq_along(groups), function(gi) {
      do.call(rbind, lapply(seq_along(pair_idx), function(pidx) {
        d1 <- pair_idx[[pidx]][1]; d2 <- pair_idx[[pidx]][2]
        do.call(rbind, lapply(seq_along(c("RR", "RS", "SR", "SS")), function(ci) {
          cell <- c("RR", "RS", "SR", "SS")[ci]
          s <- .summ(pair_draws[, gi, pidx, ci])
          data.frame(group = groups[gi], class_a = classes[d1], class_b = classes[d2], cell = cell,
                     mean = s[["mean"]], lower = s[["lower"]], upper = s[["upper"]], stringsAsFactors = FALSE)
        }))
      }))
    }))
  } else {
    data.frame(group = character(0), class_a = character(0), class_b = character(0), cell = character(0),
               mean = numeric(0), lower = numeric(0), upper = numeric(0), stringsAsFactors = FALSE)
  }

  extreme_profiles <- do.call(rbind, lapply(seq_along(groups), function(gi) {
    sR <- .summ(allR_draws[, gi])
    sS <- .summ(allS_draws[, gi])
    marg_mean <- vapply(seq_len(D), function(d) mean(class_margin_draws[, gi, d]), numeric(1))
    indep_R <- prod(marg_mean)
    indep_S <- prod(1 - marg_mean)
    data.frame(
      group = groups[gi],
      all_resistant_mean = sR[["mean"]], all_resistant_lower = sR[["lower"]], all_resistant_upper = sR[["upper"]],
      all_susceptible_mean = sS[["mean"]], all_susceptible_lower = sS[["lower"]], all_susceptible_upper = sS[["upper"]],
      independence_all_resistant = indep_R, independence_all_susceptible = indep_S,
      stringsAsFactors = FALSE
    )
  }))

  count_distribution <- do.call(rbind, lapply(seq_along(groups), function(gi) {
    do.call(rbind, lapply(0:D, function(c) {
      s <- .summ(count_draws[, gi, c + 1L])
      data.frame(group = groups[gi], count = c, mean = s[["mean"]], lower = s[["lower"]], upper = s[["upper"]], stringsAsFactors = FALSE)
    }))
  }))

  max_entropy <- log(S)
  concentration <- do.call(rbind, lapply(seq_along(groups), function(gi) {
    me <- mean(entropy_draws[, gi])
    data.frame(group = groups[gi], mean_entropy = me, max_entropy = max_entropy,
               entropy_ratio = me / max_entropy, stringsAsFactors = FALSE)
  }))

  list(
    class_margins = class_margins, pair_probabilities = pair_probabilities,
    extreme_profiles = extreme_profiles, count_distribution = count_distribution,
    concentration = concentration
  )
}
