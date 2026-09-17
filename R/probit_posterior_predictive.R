# probit_posterior_predictive.R
#
# Posterior predictive checking for a fitted Bayesian multivariate probit
# model (fit_bayesian_multivariate_probit()). This answers a third, distinct
# question from the other two families of model assessment already in this
# package:
#
#   diagnostic_status           -- "did the HMC/NUTS sampler behave well?"
#                                   (computed inline in
#                                   fit_bayesian_multivariate_probit(); see
#                                   daly_resistance_profiles.R)
#   profile_validation_status   -- "do fitted/implied cell probabilities
#                                   agree with observed rates?" (marginal/
#                                   pairwise/complete-profile calibration and
#                                   masked-AST refit; see
#                                   daly_resistance_validation.R)
#   posterior_predictive_status -- "can the fitted model generate WHOLE
#                                   REPLICATED DATASETS that resemble the one
#                                   actually observed?" (this file)
#
# A model can pass every marginal-cell calibration check while still failing
# to reproduce joint multidrug-resistance burden, hospital-to-hospital
# heterogeneity, or within-admission/within-patient clustering -- exactly the
# kind of misfit a posterior predictive check is designed to catch. See also
# prior_predictive_status (probit_prior_predictive.R) and mixed predictive
# checks (probit_mixed_predictive.R).
#
# These status fields are NEVER merged into one generic "model_status"
# anywhere in this package. A downstream decision layer (in
# anumaan-analysis) may combine them for model selection / DALY eligibility;
# anumaan itself only ever exposes them separately.
#
# Methodological reference: Stan User's Guide, "Posterior and Prior
# Predictive Checks":
# https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html
#
#   y_rep ~ p(y_rep | y) = integral p(y_rep | theta) p(theta | y) d theta
#
# CRITICAL DISTINCTION from the DALY profile-completion machinery
# (compute_event_profile_probabilities() / the internal
# .gibbs_conditional_profile_probs()): those functions sample the MISSING
# cells of an event's resistance profile CONDITIONAL on that event's own
# OBSERVED AST outcomes (obs_panel) -- an "observed-plus-imputed" estimand,
# correct for DALY use. Posterior predictive replication draws ALL classes'
# outcomes unconditionally, Y_rep ~ P(Y | theta), for a theta drawn from the
# posterior p(theta | Y) -- observed data influences Y_rep only THROUGH
# theta, never by conditioning the replicate itself on the actually-observed
# values. The functions in this file therefore never call
# .gibbs_conditional_profile_probs() (see
# test-probit-posterior-predictive.R for an explicit regression test).
#
# "ppc_tail_probability" columns produced by compute_probit_ppc_statistics()
# are posterior-predictive tail-probability-like quantities
# (mean(T_rep >= T_obs)), NOT classical calibrated hypothesis-test p-values.
# No object in this file is ever named literally "p_value".
#
# "top_observed_profiles" vs. "clinical_priority_profiles": the
# complete_profile statistic family's per-profile check
# (profile_top_observed_frequency) is a purely DATA-DOMINANT summary -- the
# top-3 most frequent profiles actually observed -- and deliberately makes no
# claim of clinical importance. A hospital's top-3 observed profiles may or
# may not include, say, carbapenem/3GC/fluoroquinolone co-resistance patterns
# that matter clinically. A future clinical_priority_profiles family, driven
# by an explicit burden-framework/expert-rule table (not yet implemented),
# would be ADDITIVE to this one, not a replacement -- the two answer
# different questions ("what does the data show most often?" vs. "how well
# does the model reproduce the patterns clinicians care about?").
#
# compute_posterior_predictive_status()'s default thresholds (tail_warning =
# 0.025, tail_severe = 0.005, max_fraction_core_extreme = 0.20) are INITIAL
# AMR-PROJECT DEFAULTS, not universal truths -- they are a reasonable
# starting point validated only against this file's synthetic recovery
# scenarios (test-probit-predictive-synthetic-recovery.R), and should be
# recalibrated once real multi-pathogen experiment results are available.

# ---------------------------------------------------------------------------
# Internal: shared posterior-draw setup for predictive-check functions
# ---------------------------------------------------------------------------

#' Shared posterior-draw setup for posterior predictive simulation
#'
#' Mirrors the draw-subsampling / mu-reconstruction pattern of
#' \code{.probit_validation_draws_setup()} in \code{daly_resistance_validation.R}
#' (same \code{posterior::as_draws_matrix()} subsampling, same \code{.arr()}
#' unpacking, same \code{re_contribution()}-based mu reconstruction), but is
#' intentionally a SEPARATE helper rather than a shared/reused one: it needs
#' the raw \code{L_Omega} Cholesky factor (not the reconstructed/regularised
#' \code{Omega} matrix the validation module needs for \code{pbivnorm}), and
#' \code{daly_resistance_validation.R} currently has zero test coverage, so
#' keeping predictive-check code isolated from it avoids introducing any risk
#' of silently changing existing validation behaviour.
#' @keywords internal
.probit_predictive_draws_setup <- function(fitted_model, n_states, seed) {
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Package 'posterior' is required (installed with cmdstanr).", call. = FALSE)

  set.seed(as.integer(seed))

  draws        <- fitted_model$draws
  class_cols   <- fitted_model$class_cols
  event_meta   <- fitted_model$event_metadata
  upper_re_col <- fitted_model$upper_re_col
  pathogen_col <- fitted_model$pathogen_col
  re_prep      <- fitted_model$random_effects_prep
  residual_structure <- .null_default(fitted_model$residual_structure, "identity")

  D <- length(fitted_model$index_maps$class_levels)
  X_event_sim  <- if (!is.null(fitted_model$X_event)) fitted_model$X_event else fitted_model$X_design
  event_re_idx <- fitted_model$event_re_idx    # N_events x R flattened level index
  K <- ncol(X_event_sim)

  draws_mat <- posterior::as_draws_matrix(draws)
  n_total   <- nrow(draws_mat)
  S <- min(as.integer(n_states), n_total)
  if (is.na(S) || S < 1L)
    stop("`n_states` must resolve to a positive integer.", call. = FALSE)
  draw_idx <- if (S < n_total) sort(sample.int(n_total, S)) else seq_len(n_total)
  draws_mat <- draws_mat[draw_idx, , drop = FALSE]

  .arr <- function(prefix, d1, d2) {
    cols <- as.vector(outer(seq_len(d1), seq_len(d2),
                            function(a, b) sprintf("%s[%d,%d]", prefix, a, b)))
    array(draws_mat[, cols, drop = FALSE], dim = c(S, d1, d2))
  }
  beta_arr    <- .arr("beta", K, D)
  re_eff_arr  <- if (re_prep$R > 0L) .arr("re_effect", D, re_prep$total_re_levels) else NULL
  L_omega_arr <- if (identical(residual_structure, "correlated")) .arr("L_Omega", D, D) else NULL

  stopifnot(".event_idx" %in% names(event_meta))
  has_obs <- event_meta$.event_idx %in% fitted_model$data_long$ev_idx
  event_meta_obs <- event_meta[has_obs, , drop = FALSE]
  N_ev <- nrow(event_meta_obs)
  ev_row_idx <- event_meta_obs$.event_idx
  X_event <- X_event_sim[ev_row_idx, , drop = FALSE]
  flat_re_idx_obs <- event_re_idx[ev_row_idx, , drop = FALSE]   # N_ev x R

  obs_ast_mat <- as.matrix(event_meta_obs[, class_cols, drop = FALSE])
  storage.mode(obs_ast_mat) <- "double"
  # Observation mask: TRUE where the real AST matrix has a tested (non-NA)
  # value for that event/class. Part 3's "preserve_observation_mask" applies
  # this to the COMPLETE generated Y_rep after the fact -- testing-coverage
  # differences must not be misread as resistance-model misfit.
  obs_mask <- !is.na(obs_ast_mat)

  # Generic random-effect contribution: sums over an arbitrary number of
  # declared blocks via the shared re_contribution() helper -- never
  # hardcodes hospital/admission/patient.
  mu_all_for_draw <- function(s) {
    beta_s   <- matrix(beta_arr[s, , ],   nrow = K, ncol = D)
    re_term <- if (re_prep$R > 0L) {
      re_eff_s <- matrix(re_eff_arr[s, , ], nrow = D, ncol = re_prep$total_re_levels)
      re_contribution(re_eff_s, flat_re_idx_obs)
    } else matrix(0, nrow = N_ev, ncol = D)
    (X_event %*% beta_s) + re_term
  }
  # Raw stored Cholesky factor of the residual correlation matrix (a Stan
  # cholesky_factor_corr[D] parameter -- already a valid Cholesky factor, no
  # tcrossprod()/re-Cholesky reconstruction needed for generative sampling).
  L_omega_for_draw <- if (identical(residual_structure, "correlated")) {
    function(s) matrix(L_omega_arr[s, , ], nrow = D, ncol = D)
  } else NULL

  list(
    S = S, N_ev = N_ev, D = D, K = K,
    class_cols = class_cols,
    event_meta_obs = event_meta_obs,
    obs_ast_mat = obs_ast_mat, obs_mask = obs_mask,
    flat_re_idx_obs = flat_re_idx_obs,
    random_effects_prep = re_prep,
    upper_re_col = upper_re_col, pathogen_col = pathogen_col,
    residual_structure = residual_structure,
    eligibility_report = fitted_model$eligibility_report,
    mu_all_for_draw = mu_all_for_draw,
    L_omega_for_draw = L_omega_for_draw
  )
}

# ---------------------------------------------------------------------------
# Internal: unconditional generative draws (identity / correlated residual)
# ---------------------------------------------------------------------------

#' Identity-residual posterior predictive draw: Y_rep ~ Bernoulli(Phi(mu))
#'
#' Distributionally identical to drawing Z ~ N(mu, 1) and thresholding at 0,
#' but avoids materialising the unneeded N_events x D latent Z matrix and
#' mirrors what \code{validate_marginal_calibration()} already does with
#' \code{stats::pnorm(mu)} for the identity case.
#' @keywords internal
.ppc_generate_identity <- function(mu_s) {
  p <- stats::pnorm(as.vector(mu_s))
  matrix(stats::rbinom(length(p), 1L, p), nrow = nrow(mu_s), ncol = ncol(mu_s))
}

#' Correlated-residual posterior predictive draw: Z_rep ~ MVN(mu, Omega),
#' Y_rep = I(Z_rep > 0)
#'
#' Reuses the state's already-stored Cholesky factor \code{L_omega_s}
#' directly (one Cholesky per state, per Part 20's performance requirement)
#' via a single D x N_events matrix multiply, vectorised across every event
#' in that state -- no per-event loop.
#' @keywords internal
.ppc_generate_correlated <- function(mu_s, L_omega_s) {
  N_ev <- nrow(mu_s)
  D    <- ncol(mu_s)
  eps <- matrix(stats::rnorm(D * N_ev), nrow = D, ncol = N_ev)
  z <- mu_s + t(L_omega_s %*% eps)
  matrix(as.integer(z > 0), nrow = N_ev, ncol = D)
}

#' Deterministic per-state seed, a function of (seed, s) only
#'
#' Ensures a given state's replicate depends ONLY on (seed, s, theta_s) --
#' never on call order, on interleaving with other generator objects built
#' from the same or a different seed, or on any other RNG consumption
#' elsewhere in the calling session between generator construction and the
#' first \code{generate_state(s)} call. Small differences in the seed
#' passed to \code{set.seed()} produce uncorrelated Mersenne Twister
#' streams, so a simple affine combination is sufficient here (this is not
#' a cryptographic hash, just a reproducible index-to-seed map).
#' @keywords internal
.ppc_state_seed <- function(seed, s) {
  as.integer((as.numeric(seed) * 1000003 + as.numeric(s)) %% 2147483647L)
}

#' Run \code{expr} under a temporary, locally-seeded RNG state, restoring the
#' caller's global RNG state (\code{.Random.seed}) afterwards -- so a
#' predictive-check generator never disturbs any OTHER stochastic code
#' running in the same session before or after it is used.
#' @keywords internal
.ppc_with_local_seed <- function(seed, expr) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  force(expr)
}

# ---------------------------------------------------------------------------
# simulate_probit_posterior_predictive()
# ---------------------------------------------------------------------------

#' Simulate posterior predictive replicate datasets from a fitted probit model
#'
#' For each of \code{n_states} saved posterior states \eqn{s} and each fitted
#' event \eqn{e}, reconstructs
#' \deqn{\mu_e^{(s)} = X_e \beta^{(s)} + \sum_r u_{r,g_r(e)}^{(s)}}
#' using the existing generic random-effect machinery
#' (\code{\link{prepare_random_effects}}/\code{\link{re_contribution}} --
#' never hardcoded block names), then draws a COMPLETE replicate outcome
#' matrix \eqn{Y_{rep}^{(s)}} unconditionally from \eqn{P(Y \mid \mu^{(s)},
#' \Omega^{(s)})}: \code{Bernoulli(Phi(mu))} per class for identity residual
#' models, or \code{Z ~ MVN(mu, Omega); Y = I(Z > 0)} for correlated residual
#' models (see Stan User's Guide,
#' \url{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}).
#'
#' This function NEVER calls the internal DALY conditional-completion helper
#' (\code{.gibbs_conditional_profile_probs()}) -- that function samples
#' missing cells conditional on an event's own observed AST signs, which is a
#' fundamentally different (and, for this purpose, wrong) generative
#' mechanism. See the file header of \code{probit_posterior_predictive.R} and
#' \code{test-probit-posterior-predictive.R} for the explicit distinction and
#' regression test.
#'
#' By default (\code{return_replicates = FALSE}), no \code{S x N_events x D}
#' array is ever materialised -- this function returns a lightweight
#' generator object (class \code{"amr_ppc_draws"}) that
#' \code{\link{compute_probit_ppc_statistics}} drives one state at a time,
#' streaming, to bound peak memory at O(N_events x D) rather than
#' O(S x N_events x D). Each state's replicate is generated under a
#' temporary, locally-seeded RNG state derived deterministically from
#' \code{(seed, s)} alone (restoring the caller's global RNG state
#' afterwards) -- so a given state's draw depends only on \code{seed},
#' \code{s}, and that state's posterior parameters, never on call order, on
#' interleaving with any other generator object, or on any other RNG use
#' elsewhere in the calling session. Draws are cached per state, so calling
#' the generator repeatedly for the same \code{s} is idempotent.
#'
#' \code{return_replicates = TRUE} additionally materialises the full
#' replicate arrays and is intended for small debugging datasets only; a
#' safety ceiling refuses to silently run this on a large realistic dataset
#' (reduce \code{n_states} or leave \code{return_replicates = FALSE}
#' instead).
#'
#' @param fitted_model A fitted model object as returned by
#'   \code{\link{fit_bayesian_multivariate_probit}} (or the equivalent
#'   \code{fit_light} object saved by the analysis layer with \code{$fit}
#'   stripped -- this function never touches \code{$fit}, only stored draws
#'   and design matrices).
#' @param n_states Integer; number of saved posterior states to use (a
#'   random subsample, seeded, if fewer than the total number of saved
#'   draws). Default 500; use 1000-2000 for shortlisted models.
#' @param seed Integer seed for state subsampling and stochastic generation.
#' @param preserve_observation_mask Logical, default \code{TRUE}. If
#'   \code{TRUE}, the returned generator's \code{Y_rep} slot masks the
#'   complete generated replicate to exactly the same tested/untested cell
#'   pattern as the real observed AST matrix, so that discrepancy statistics
#'   computed downstream compare like with like. The COMPLETE (unmasked)
#'   replicate is always separately available as \code{Y_rep_complete} for
#'   statistics that genuinely need it (e.g. complete-profile statistics
#'   restricted to adequately-supported fully-observed subsets) -- the two
#'   are never compared without explicit qualification.
#' @param return_replicates Logical, default \code{FALSE}. If \code{TRUE},
#'   materialises full \code{S x N_events x D} replicate arrays
#'   (\code{Y_rep_array}, masked per \code{preserve_observation_mask}, and
#'   \code{Y_rep_complete_array}, always unmasked) -- small debug datasets
#'   only; refuses above a fixed element-count ceiling.
#' @param replicate_random_effects,random_effect_blocks_to_replicate Reserved
#'   for API compatibility with the (distinct) mixed predictive check. Not
#'   implemented here -- passing a non-default value errors with a pointer to
#'   \code{\link{simulate_probit_mixed_predictive}} rather than silently
#'   being ignored.
#' @return An object of class \code{"amr_ppc_draws"}: a list with
#'   \code{setup} (internal draw-setup object), \code{generate_state}
#'   (a function of one integer argument \code{s} returning
#'   \code{list(Y_rep_complete, Y_rep)} for that state), \code{n_states_used},
#'   \code{seed_used}, \code{preserve_observation_mask}, \code{residual_structure},
#'   \code{return_replicates}, and (only when \code{return_replicates = TRUE})
#'   \code{Y_rep_array}/\code{Y_rep_complete_array}.
#' @references
#' Stan Development Team. "Posterior and Prior Predictive Checks." Stan
#' User's Guide. \url{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}
#' @export
simulate_probit_posterior_predictive <- function(
    fitted_model,
    n_states = 500L,
    seed = 123L,
    preserve_observation_mask = TRUE,
    return_replicates = FALSE,
    replicate_random_effects = FALSE,
    random_effect_blocks_to_replicate = NULL
) {
  if (isTRUE(replicate_random_effects) || !is.null(random_effect_blocks_to_replicate)) {
    stop(
      "replicate_random_effects/random_effect_blocks_to_replicate are not implemented in ",
      "simulate_probit_posterior_predictive() -- generating NEW random-effect levels while ",
      "retaining posterior hyperparameters is a MIXED predictive check, a distinct concept ",
      "from a standard posterior predictive check (which conditions on the fitted random ",
      "effects exactly as estimated). Use simulate_probit_mixed_predictive() instead.",
      call. = FALSE
    )
  }

  setup <- .probit_predictive_draws_setup(fitted_model, n_states, seed)

  generate_state <- local({
    cache <- new.env(parent = emptyenv())
    function(s) {
      key <- as.character(as.integer(s))
      cached <- cache[[key]]
      if (!is.null(cached)) return(cached)

      val <- .ppc_with_local_seed(.ppc_state_seed(seed, s), {
        mu_s <- setup$mu_all_for_draw(s)
        Y_complete <- if (identical(setup$residual_structure, "correlated")) {
          .ppc_generate_correlated(mu_s, setup$L_omega_for_draw(s))
        } else {
          .ppc_generate_identity(mu_s)
        }
        Y_masked <- Y_complete
        Y_masked[!setup$obs_mask] <- NA_integer_
        list(
          Y_rep_complete = Y_complete,
          Y_rep = if (isTRUE(preserve_observation_mask)) Y_masked else Y_complete
        )
      })
      cache[[key]] <- val
      val
    }
  })

  out <- list(
    setup = setup,
    generate_state = generate_state,
    n_states_used = setup$S,
    seed_used = as.integer(seed),
    preserve_observation_mask = isTRUE(preserve_observation_mask),
    residual_structure = setup$residual_structure,
    return_replicates = isTRUE(return_replicates),
    Y_rep_array = NULL,
    Y_rep_complete_array = NULL
  )

  if (isTRUE(return_replicates)) {
    n_elements <- as.double(setup$S) * setup$N_ev * setup$D
    max_elements <- 5e7
    if (n_elements > max_elements) {
      stop(sprintf(
        paste(
          "return_replicates = TRUE would materialize %.0f array elements",
          "(n_states=%d x N_events=%d x D=%d), above the safety ceiling of %.0f.",
          "This path is intended for small debug datasets only -- reduce n_states,",
          "or leave return_replicates = FALSE and use compute_probit_ppc_statistics()",
          "for the compact streaming summary instead."
        ),
        n_elements, setup$S, setup$N_ev, setup$D, max_elements
      ), call. = FALSE)
    }
    Y_rep_array          <- array(NA_integer_, dim = c(setup$S, setup$N_ev, setup$D))
    Y_rep_complete_array <- array(NA_integer_, dim = c(setup$S, setup$N_ev, setup$D))
    for (s in seq_len(setup$S)) {
      st <- generate_state(s)
      Y_rep_complete_array[s, , ] <- st$Y_rep_complete
      Y_rep_array[s, , ] <- st$Y_rep
    }
    out$Y_rep_array <- Y_rep_array
    out$Y_rep_complete_array <- Y_rep_complete_array
  }

  class(out) <- "amr_ppc_draws"
  out
}

# ---------------------------------------------------------------------------
# compute_probit_ppc_statistics(): AMR-specific discrepancy statistics
# ---------------------------------------------------------------------------
#
# Architecture: a one-time "spec-building" pass (.ppc_build_spec()) resolves,
# for every requested statistic family, a fixed list of "items" -- one per
# (statistic_name, stratum) pair -- each carrying its observed value, support
# status, and a `compute(Y_rep, Y_rep_complete)` closure. The actual S-state
# loop (in compute_probit_ppc_statistics() itself) is then a single pass that
# calls generate_state(s) ONCE per state and feeds the resulting replicate
# into every item's closure, preallocating one S x n_items numeric matrix --
# this is what bounds the state loop to O(N_events x D) transient memory per
# state (Part 20) rather than re-deriving mu/Y_rep per statistic per state.
#
# Every closure created inside a `for` loop below is wrapped in local({...})
# so it captures a FRESH copy of any loop variable it needs at the time the
# item is constructed -- R's `for` loops reuse a single frame across
# iterations, so an unwrapped closure would otherwise see only the LAST
# iteration's values when it is actually invoked later, in the state loop.

#' @keywords internal
.ppc_quantile <- function(x, p) unname(stats::quantile(x, probs = p, na.rm = TRUE, type = 7))

#' Reduce one item's (observed_value, T_rep across states) to the Part 7
#' output schema: statistic_name, stratum, observed_value, replicated_mean/
#' sd/q025/q50/q975, ppc_tail_probability, ppc_two_sided, n_replications,
#' support_status. ppc_tail_probability/ppc_two_sided are posterior
#' predictive TAIL-PROBABILITY-LIKE quantities. The reported one-sided value
#' is the inclusive upper tail, \code{(1 + sum(T_rep >= T_obs)) / (B + 1)};
#' the two-sided value is \code{min(1, 2 * min(lower_tail, upper_tail))}, with
#' the lower tail calculated analogously. These finite-replicate corrected
#' summaries are NOT classical calibrated p-values.
#' @keywords internal
.ppc_summarize_statistic <- function(statistic_name, stratum, observed_value, T_rep,
                                      support_status, ci_level) {
  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q
  T_rep <- T_rep[!is.na(T_rep)]
  n_rep <- length(T_rep)

  if (!identical(support_status, "supported") || is.na(observed_value) || n_rep == 0L) {
    return(tibble::tibble(
      statistic_name = statistic_name, stratum = stratum,
      observed_value = observed_value,
      replicated_mean = NA_real_, replicated_sd = NA_real_,
      replicated_q025 = NA_real_, replicated_q50 = NA_real_, replicated_q975 = NA_real_,
      ppc_tail_probability = NA_real_, ppc_two_sided = NA_real_,
      n_replications = n_rep, support_status = support_status
    ))
  }

  lower_tail <- (1 + sum(T_rep <= observed_value)) / (n_rep + 1)
  upper_tail <- (1 + sum(T_rep >= observed_value)) / (n_rep + 1)
  two_sided <- min(1, 2 * min(lower_tail, upper_tail))

  tibble::tibble(
    statistic_name = statistic_name, stratum = stratum,
    observed_value = observed_value,
    replicated_mean = mean(T_rep), replicated_sd = stats::sd(T_rep),
    replicated_q025 = .ppc_quantile(T_rep, lo_q),
    replicated_q50  = .ppc_quantile(T_rep, 0.5),
    replicated_q975 = .ppc_quantile(T_rep, hi_q),
    ppc_tail_probability = upper_tail, ppc_two_sided = two_sided,
    n_replications = n_rep, support_status = support_status
  )
}

# --- Part 5A: marginal resistance -------------------------------------------

#' @keywords internal
.ppc_items_marginal <- function(setup) {
  items <- list()
  D <- setup$D; class_cols <- setup$class_cols
  obs <- setup$obs_ast_mat; mask <- setup$obs_mask

  # Global (pooled, dataset-wide) per-class marginal -- support thresholds
  # mirror .validate_panel_eligibility()'s own defaults (30/5/5).
  for (d in seq_len(D)) {
    cc <- class_cols[d]
    tested_idx <- which(mask[, d])
    n_tested <- length(tested_idx)
    n_res <- sum(obs[tested_idx, d] == 1)
    n_sus <- sum(obs[tested_idx, d] == 0)
    supported <- n_tested >= 30L && n_res >= 5L && n_sus >= 5L
    obs_val <- if (n_tested > 0L) n_res / n_tested else NA_real_
    items[[length(items) + 1L]] <- local({
      idx <- tested_idx; dd <- d
      list(
        statistic_name = "marginal_resistance", stratum = sprintf("class:%s", cc),
        observed_value = obs_val,
        support_status = if (supported) "supported" else "insufficient_support",
        compute = function(Y_rep, Y_rep_complete) mean(Y_rep[idx, dd] == 1)
      )
    })
  }

  # Hospital x pathogen x class, reusing eligibility_report$marginal exactly
  # as validate_marginal_calibration() does.
  elig <- setup$eligibility_report$marginal
  if (!is.null(elig) && nrow(elig) > 0L) {
    upper_re_col <- setup$upper_re_col; pathogen_col <- setup$pathogen_col
    event_meta <- setup$event_meta_obs
    elig_rows <- elig[elig$eligible, , drop = FALSE]
    for (r in seq_len(nrow(elig_rows))) {
      h <- elig_rows[[upper_re_col]][r]; k <- elig_rows[[pathogen_col]][r]
      cc <- elig_rows$antibiotic_class[r]
      d <- match(cc, class_cols)
      if (is.na(d)) next
      sub_idx <- which(event_meta[[upper_re_col]] == h & event_meta[[pathogen_col]] == k & mask[, d])
      n_tested <- length(sub_idx)
      if (n_tested == 0L) next
      n_res <- sum(obs[sub_idx, d] == 1)
      obs_val <- n_res / n_tested
      items[[length(items) + 1L]] <- local({
        idx <- sub_idx; dd <- d
        list(
          statistic_name = "marginal_resistance",
          stratum = sprintf("hospital:%s|class:%s", h, cc),
          observed_value = obs_val, support_status = "supported",
          compute = function(Y_rep, Y_rep_complete) mean(Y_rep[idx, dd] == 1)
        )
      })
    }
  }
  items
}

# --- Part 5B: resistant classes per event -----------------------------------

#' @keywords internal
.ppc_items_resistant_count <- function(setup, min_complete_profile_events) {
  items <- list()
  N_ev <- setup$N_ev
  obs <- setup$obs_ast_mat; mask <- setup$obs_mask

  tested_list <- lapply(seq_len(N_ev), function(i) which(mask[i, ]))
  C_obs <- vapply(seq_len(N_ev), function(i) {
    idx <- tested_list[[i]]
    if (length(idx) == 0L) return(NA_real_)
    sum(obs[i, idx] == 1)
  }, numeric(1L))
  has_tested <- !is.na(C_obs)

  summary_fns <- list(
    resistant_count_mean     = function(x) mean(x),
    resistant_count_variance = function(x) stats::var(x),
    resistant_count_median   = function(x) stats::median(x),
    resistant_count_p90      = function(x) .ppc_quantile(x, 0.9),
    resistant_count_max      = function(x) max(x)
  )
  for (nm in names(summary_fns)) {
    fn <- summary_fns[[nm]]
    items[[length(items) + 1L]] <- local({
      fn2 <- fn; tl <- tested_list; ht <- has_tested
      list(
        statistic_name = nm, stratum = "all_events",
        observed_value = fn2(C_obs[ht]),
        support_status = if (sum(ht) >= 2L) "supported" else "insufficient_support",
        compute = function(Y_rep, Y_rep_complete) {
          Crep <- vapply(which(ht), function(i) {
            idx <- tl[[i]]
            sum(Y_rep[i, idx] == 1, na.rm = TRUE)
          }, numeric(1L))
          fn2(Crep)
        }
      )
    })
  }

  # Complete-panel per hospital x pathogen: full C = 0..Dp count distribution
  # (the "true full count distribution" for fully-observed events only).
  hp_pairs <- unique(setup$event_meta_obs[, c(setup$upper_re_col, setup$pathogen_col), drop = FALSE])
  for (r in seq_len(nrow(hp_pairs))) {
    h <- hp_pairs[[setup$upper_re_col]][r]; k <- hp_pairs[[setup$pathogen_col]][r]
    panel <- .resolve_profile_class_panel(
      class_cols = setup$class_cols, hospital = h, pathogen = k,
      eligibility_report = setup$eligibility_report,
      upper_re_col = setup$upper_re_col, pathogen_col = setup$pathogen_col,
      residual_structure = setup$residual_structure
    )
    panel_classes <- panel$classes
    Dp <- length(panel_classes)
    if (Dp == 0L) next
    d_idx <- match(panel_classes, setup$class_cols)
    hp_idx <- which(setup$event_meta_obs[[setup$upper_re_col]] == h &
                     setup$event_meta_obs[[setup$pathogen_col]] == k)
    complete_idx <- hp_idx[rowSums(mask[hp_idx, d_idx, drop = FALSE]) == Dp]
    n_complete <- length(complete_idx)
    supported <- n_complete >= min_complete_profile_events
    C_obs_complete <- if (n_complete > 0L) rowSums(obs[complete_idx, d_idx, drop = FALSE] == 1) else integer(0)
    for (cnt in 0:Dp) {
      obs_val <- if (supported) mean(C_obs_complete == cnt) else NA_real_
      items[[length(items) + 1L]] <- local({
        idx <- complete_idx; di <- d_idx; cc <- cnt
        list(
          statistic_name = "resistant_count_proportion",
          stratum = sprintf("hospital:%s|pathogen:%s|panel_size:%d|C=%d", h, k, Dp, cc),
          observed_value = obs_val,
          support_status = if (supported) "supported" else "insufficient_support",
          compute = function(Y_rep, Y_rep_complete) {
            if (length(idx) == 0L) return(NA_real_)
            Crep <- rowSums(Y_rep[idx, di, drop = FALSE] == 1, na.rm = TRUE)
            mean(Crep == cc)
          }
        )
      })
    }
  }
  items
}

# --- Part 5C: pairwise co-resistance -----------------------------------------

#' @keywords internal
.ppc_items_pairwise <- function(setup) {
  items <- list()
  D <- setup$D
  if (D < 2L) return(items)
  class_cols <- setup$class_cols
  obs <- setup$obs_ast_mat; mask <- setup$obs_mask
  pw <- setup$eligibility_report$pairwise
  if (is.null(pw) || nrow(pw) == 0L) return(items)
  upper_re_col <- setup$upper_re_col
  event_meta <- setup$event_meta_obs
  sufficient_rows <- pw[pw$sufficient, , drop = FALSE]

  for (r in seq_len(nrow(sufficient_rows))) {
    h <- sufficient_rows[[upper_re_col]][r]
    c1 <- sufficient_rows$class_1[r]; c2 <- sufficient_rows$class_2[r]
    d1 <- match(c1, class_cols); d2 <- match(c2, class_cols)
    if (is.na(d1) || is.na(d2)) next
    sub_idx <- which(event_meta[[upper_re_col]] == h)
    cotested_idx <- sub_idx[mask[sub_idx, d1] & mask[sub_idx, d2]]
    n_cotested <- length(cotested_idx)
    if (n_cotested == 0L) next
    v1 <- obs[cotested_idx, d1]; v2 <- obs[cotested_idx, d2]
    obs_props <- c(
      RR = mean(v1 == 1 & v2 == 1), RS = mean(v1 == 1 & v2 == 0),
      SR = mean(v1 == 0 & v2 == 1), SS = mean(v1 == 0 & v2 == 0)
    )
    for (lab in c("RR", "RS", "SR", "SS")) {
      items[[length(items) + 1L]] <- local({
        idx <- cotested_idx; dd1 <- d1; dd2 <- d2; l <- lab
        list(
          statistic_name = sprintf("pairwise_%s", l),
          stratum = sprintf("hospital:%s|%s_x_%s", h, c1, c2),
          observed_value = unname(obs_props[l]), support_status = "supported",
          compute = function(Y_rep, Y_rep_complete) {
            a <- Y_rep[idx, dd1]; b <- Y_rep[idx, dd2]
            switch(l,
              RR = mean(a == 1 & b == 1), RS = mean(a == 1 & b == 0),
              SR = mean(a == 0 & b == 1), SS = mean(a == 0 & b == 0)
            )
          }
        )
      })
    }
  }
  items
}

# --- Part 5D: complete-profile statistics ------------------------------------

#' @keywords internal
.ppc_items_complete_profile <- function(setup, min_complete_profile_events) {
  items <- list()
  hp_pairs <- unique(setup$event_meta_obs[, c(setup$upper_re_col, setup$pathogen_col), drop = FALSE])
  obs <- setup$obs_ast_mat; mask <- setup$obs_mask
  label_of <- function(mat_slice) {
    apply(mat_slice, 1L, function(row) paste0(ifelse(row == 1, "R", "S"), collapse = ""))
  }
  shannon <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  simpson <- function(p) sum(p^2)

  for (r in seq_len(nrow(hp_pairs))) {
    h <- hp_pairs[[setup$upper_re_col]][r]; k <- hp_pairs[[setup$pathogen_col]][r]
    panel <- .resolve_profile_class_panel(
      class_cols = setup$class_cols, hospital = h, pathogen = k,
      eligibility_report = setup$eligibility_report,
      upper_re_col = setup$upper_re_col, pathogen_col = setup$pathogen_col,
      residual_structure = setup$residual_structure
    )
    panel_classes <- panel$classes
    Dp <- length(panel_classes)
    if (Dp == 0L) next
    d_idx <- match(panel_classes, setup$class_cols)
    hp_idx <- which(setup$event_meta_obs[[setup$upper_re_col]] == h &
                     setup$event_meta_obs[[setup$pathogen_col]] == k)
    complete_idx <- hp_idx[rowSums(mask[hp_idx, d_idx, drop = FALSE]) == Dp]
    n_complete <- length(complete_idx)
    supported <- n_complete >= min_complete_profile_events
    strat_base <- sprintf("hospital:%s|pathogen:%s", h, k)

    if (!supported) {
      for (nm in c("profile_all_susceptible_frequency", "profile_all_resistant_frequency",
                   "profile_most_common_frequency", "profile_n_distinct",
                   "profile_shannon_entropy", "profile_simpson_concentration")) {
        items[[length(items) + 1L]] <- local({
          nm2 <- nm; sb <- strat_base
          list(statistic_name = nm2, stratum = sb, observed_value = NA_real_,
               support_status = "insufficient_support",
               compute = function(Y_rep, Y_rep_complete) NA_real_)
        })
      }
      next
    }

    obs_labels <- label_of(obs[complete_idx, d_idx, drop = FALSE])
    obs_tab <- table(obs_labels)
    all_s_label <- paste(rep("S", Dp), collapse = "")
    all_r_label <- paste(rep("R", Dp), collapse = "")
    mode_label <- names(obs_tab)[which.max(obs_tab)]
    obs_freq <- as.numeric(obs_tab) / n_complete

    freq_of <- function(labels, lab) mean(labels == lab)

    # Item builder: every free variable this closure needs (idx, di, and the
    # fixed comparison label/reducer) is snapshotted into a FRESH local()
    # environment at construction time, so later hp_pairs iterations cannot
    # retroactively change what an earlier item's compute() sees.
    mk_item <- function(statistic_name, observed_value, reducer) {
      local({
        idx <- complete_idx; di <- d_idx; red <- reducer
        list(
          statistic_name = statistic_name, stratum = strat_base,
          observed_value = observed_value, support_status = "supported",
          compute = function(Y_rep, Y_rep_complete) {
            rep_labels <- label_of(Y_rep[idx, di, drop = FALSE])
            red(rep_labels)
          }
        )
      })
    }

    items[[length(items) + 1L]] <- mk_item(
      "profile_all_susceptible_frequency", freq_of(obs_labels, all_s_label),
      local({ lab <- all_s_label; function(labels) mean(labels == lab) })
    )
    items[[length(items) + 1L]] <- mk_item(
      "profile_all_resistant_frequency", freq_of(obs_labels, all_r_label),
      local({ lab <- all_r_label; function(labels) mean(labels == lab) })
    )
    items[[length(items) + 1L]] <- mk_item(
      "profile_most_common_frequency", freq_of(obs_labels, mode_label),
      local({ lab <- mode_label; function(labels) mean(labels == lab) })
    )
    items[[length(items) + 1L]] <- mk_item(
      "profile_n_distinct", length(unique(obs_labels)),
      function(labels) length(unique(labels))
    )
    items[[length(items) + 1L]] <- mk_item(
      "profile_shannon_entropy", shannon(obs_freq),
      function(labels) { tb <- table(labels); shannon(as.numeric(tb) / length(labels)) }
    )
    items[[length(items) + 1L]] <- mk_item(
      "profile_simpson_concentration", simpson(obs_freq),
      function(labels) { tb <- table(labels); simpson(as.numeric(tb) / length(labels)) }
    )

    # Top OBSERVED profiles (data-dominant, not a clinical-importance claim):
    # the top-3 most frequent profiles actually seen in the data, purely as a
    # descriptive summary. This is deliberately NOT labelled "clinically
    # important" -- clinical priority (e.g. profiles involving resistance to
    # carbapenems/3GCs/fluoroquinolones) is a separate, currently unimplemented
    # concept that would need to come from a burden framework or expert rule
    # table, not be inferred from observed frequency. A future
    # `clinical_priority_profiles` statistic family, driven by such a table,
    # would be additive to (not a replacement for) this one.
    top_labels <- names(sort(obs_tab, decreasing = TRUE))[seq_len(min(3L, length(obs_tab)))]
    for (lab in top_labels) {
      items[[length(items) + 1L]] <- mk_item(
        "profile_top_observed_frequency", freq_of(obs_labels, lab),
        local({ l <- lab; function(labels) mean(labels == l) })
      )
      items[[length(items)]]$stratum <- sprintf("%s|profile:%s", strat_base, lab)
    }
  }
  items
}

# --- Part 5E: hospital heterogeneity -----------------------------------------

#' @keywords internal
.ppc_items_hospital_heterogeneity <- function(setup, min_hospital_support) {
  items <- list()
  D <- setup$D; class_cols <- setup$class_cols
  obs <- setup$obs_ast_mat; mask <- setup$obs_mask
  upper_re_col <- setup$upper_re_col
  event_meta <- setup$event_meta_obs
  hospitals <- sort(unique(event_meta[[upper_re_col]]))

  spread_fns <- list(
    hospital_heterogeneity_sd    = function(x) stats::sd(x),
    hospital_heterogeneity_iqr   = function(x) unname(diff(stats::quantile(x, c(0.25, 0.75), type = 7))),
    hospital_heterogeneity_mad   = function(x) stats::mad(x),
    hospital_heterogeneity_range = function(x) max(x) - min(x)
  )

  for (d in seq_len(D)) {
    cc <- class_cols[d]
    hosp_idx_list <- list()
    for (h in hospitals) {
      idx <- which(event_meta[[upper_re_col]] == h & mask[, d])
      if (length(idx) >= min_hospital_support) hosp_idx_list[[h]] <- idx
    }
    supported <- length(hosp_idx_list) >= 2L
    obs_props <- if (supported) {
      vapply(hosp_idx_list, function(idx) mean(obs[idx, d] == 1), numeric(1L))
    } else numeric(0)

    for (nm in names(spread_fns)) {
      fn <- spread_fns[[nm]]
      items[[length(items) + 1L]] <- local({
        idx_list <- hosp_idx_list; dd <- d; fn2 <- fn
        list(
          statistic_name = nm, stratum = sprintf("class:%s", cc),
          observed_value = if (supported) fn2(obs_props) else NA_real_,
          support_status = if (supported) "supported" else "insufficient_support",
          compute = function(Y_rep, Y_rep_complete) {
            if (length(idx_list) < 2L) return(NA_real_)
            props <- vapply(idx_list, function(ix) mean(Y_rep[ix, dd] == 1), numeric(1L))
            fn2(props)
          }
        )
      })
    }
  }
  items
}

# --- Part 6: random-effect-specific PPC statistics (admission/patient) ------
#
# anumaan never hardcodes "hospital"/"admission"/"patient" as DECLARED BLOCK
# NAMES -- it only recognises those three strings as ROLE values inside the
# optional `random_effect_roles` argument. Each role maps to the name of a
# GROUPING COLUMN in `fitted_model$event_metadata` (which retains every
# original column, not just those declared as random-effect blocks). This is
# DELIBERATELY decoupled from whether that column is actually modelled as a
# random-effect block in this particular fit: a clustering discrepancy
# statistic's whole purpose is to ask "do this fit's REPLICATED datasets
# reproduce the within-group similarity the OBSERVED data shows for this
# grouping?", and that question is equally meaningful -- arguably more
# useful -- for a grouping the model does NOT declare an RE block for (e.g.
# to detect exactly the omitted-clustering failure mode a fit lacking an
# admission random effect would produce). When the named column happens to
# also be a declared RE block's group_col, the fit's own random effect for
# it will generally help the replicate reproduce that clustering; when it is
# not declared at all, nothing in the generative process ties same-group
# events together, and the replicate's within-group similarity will tend to
# collapse toward the between-group baseline -- precisely the discrepancy
# this statistic is designed to surface. The analysis layer
# (anumaan-analysis) is responsible for populating `random_effect_roles`
# from its own experiment configuration/data dictionary.

#' @keywords internal
.ppc_role_column_name <- function(role, random_effect_roles) {
  if (is.null(random_effect_roles)) return(NULL)
  nm <- names(random_effect_roles)[random_effect_roles == role]
  if (length(nm) == 0L) return(NULL)
  nm[1L]
}

#' @keywords internal
.ppc_role_eligible_fallback <- function(cluster_id, min_groups_with_2plus = 20L) {
  tab <- table(cluster_id)
  sum(tab >= 2L) >= min_groups_with_2plus
}

#' @keywords internal
.ppc_role_eligible <- function(role, random_effect_eligibility_audit) {
  if (is.null(random_effect_eligibility_audit)) return(NA)
  aud <- random_effect_eligibility_audit
  if (!all(c("block", "status") %in% names(aud))) return(NA)
  row <- aud[aud$block == role, , drop = FALSE]
  if (nrow(row) == 0L) return(NA)
  isTRUE(row$status[1L] %in% c("eligible_primary", "eligible_sensitivity_only"))
}

#' @keywords internal
.ppc_cluster_unsupported_items <- function(role, stratum, reason) {
  nms <- c("cluster_within_group_same_class_agreement",
           "cluster_within_group_hamming_distance",
           "cluster_within_between_agreement_diff")
  lapply(nms, function(nm) {
    list(statistic_name = nm, stratum = stratum, observed_value = NA_real_,
         support_status = reason, compute = function(Y_rep, Y_rep_complete) NA_real_)
  })
}

#' Within-/between-cluster same-class agreement and Hamming distance for an
#' arbitrary event grouping (used for both admission and patient clustering).
#'
#' \code{cluster_within_group_same_class_agreement}: pooled agreement rate
#' (fraction of within-group event-pair x class comparisons, restricted to
#' classes BOTH events tested, where the two values agree).
#' \code{cluster_within_group_hamming_distance}: mean, over within-group
#' event pairs, of that PAIR's own mismatch rate across its jointly-tested
#' classes (a per-pair-normalised Hamming distance -- distinct from the
#' pooled agreement statistic above, which does not average per pair first).
#' \code{cluster_within_between_agreement_diff}: within-group agreement minus
#' agreement computed over a deterministic random sample of between-group
#' event pairs (sampling all possible between-group pairs is infeasible for
#' large N; the sample size scales with the number of within-group pairs).
#' @keywords internal
.ppc_cluster_pair_stats <- function(setup, role, cluster_id, restrict_diff_id = NULL,
                                     min_group_events = 2L, max_between_pairs = 5000L) {
  N_ev <- setup$N_ev
  obs <- setup$obs_ast_mat; mask <- setup$obs_mask
  stratum <- sprintf("role:%s", role)

  grp <- split(seq_len(N_ev), cluster_id)
  if (!is.null(restrict_diff_id)) {
    grp <- grp[vapply(grp, function(ix) length(unique(restrict_diff_id[ix])) >= 2L, logical(1L))]
  } else {
    grp <- grp[vapply(grp, length, integer(1L)) >= min_group_events]
  }
  if (length(grp) == 0L) return(.ppc_cluster_unsupported_items(role, stratum, "insufficient_support"))

  pair_list <- list()
  for (g in grp) {
    if (length(g) < 2L) next
    combs <- utils::combn(g, 2L, simplify = FALSE)
    if (!is.null(restrict_diff_id)) {
      combs <- Filter(function(pr) restrict_diff_id[pr[1]] != restrict_diff_id[pr[2]], combs)
    }
    pair_list <- c(pair_list, combs)
  }
  if (length(pair_list) == 0L) return(.ppc_cluster_unsupported_items(role, stratum, "insufficient_support"))
  within_i <- vapply(pair_list, `[`, integer(1L), 1L)
  within_j <- vapply(pair_list, `[`, integer(1L), 2L)

  all_idx <- unique(unlist(grp))
  n_within <- length(within_i)
  n_between_target <- min(max_between_pairs, n_within * 5L)
  between_pairs <- .ppc_with_local_seed(999L, {
    bi <- integer(0); bj <- integer(0)
    attempts <- 0L
    max_attempts <- max(n_between_target * 20L, 200L)
    while (length(bi) < n_between_target && attempts < max_attempts && length(all_idx) >= 2L) {
      attempts <- attempts + 1L
      pr <- sample(all_idx, 2L)
      if (cluster_id[pr[1]] == cluster_id[pr[2]]) next
      bi <- c(bi, pr[1]); bj <- c(bj, pr[2])
    }
    list(bi = bi, bj = bj)
  })
  between_i <- between_pairs$bi; between_j <- between_pairs$bj

  .pair_agreement <- function(Y, idx_i, idx_j) {
    mat_i <- Y[idx_i, , drop = FALSE]; mat_j <- Y[idx_j, , drop = FALSE]
    mask_i <- mask[idx_i, , drop = FALSE]; mask_j <- mask[idx_j, , drop = FALSE]
    both <- mask_i & mask_j
    n_both <- rowSums(both)
    agree <- rowSums((mat_i == mat_j) & both)
    list(agree = agree, n_both = n_both)
  }

  compute_agreement_hamming <- function(Y) {
    aw <- .pair_agreement(Y, within_i, within_j)
    valid_w <- aw$n_both > 0L
    agreement_within <- if (any(valid_w)) sum(aw$agree[valid_w]) / sum(aw$n_both[valid_w]) else NA_real_
    hamming_within <- if (any(valid_w)) mean((aw$n_both[valid_w] - aw$agree[valid_w]) / aw$n_both[valid_w]) else NA_real_

    agreement_between <- NA_real_
    if (length(between_i) > 0L) {
      ab <- .pair_agreement(Y, between_i, between_j)
      valid_b <- ab$n_both > 0L
      agreement_between <- if (any(valid_b)) sum(ab$agree[valid_b]) / sum(ab$n_both[valid_b]) else NA_real_
    }
    list(agreement_within = agreement_within, hamming_within = hamming_within,
         diff = agreement_within - agreement_between)
  }

  obs_res <- compute_agreement_hamming(obs)

  list(
    list(statistic_name = "cluster_within_group_same_class_agreement", stratum = stratum,
         observed_value = obs_res$agreement_within, support_status = "supported",
         compute = function(Y_rep, Y_rep_complete) compute_agreement_hamming(Y_rep)$agreement_within),
    list(statistic_name = "cluster_within_group_hamming_distance", stratum = stratum,
         observed_value = obs_res$hamming_within, support_status = "supported",
         compute = function(Y_rep, Y_rep_complete) compute_agreement_hamming(Y_rep)$hamming_within),
    list(statistic_name = "cluster_within_between_agreement_diff", stratum = stratum,
         observed_value = obs_res$diff, support_status = "supported",
         compute = function(Y_rep, Y_rep_complete) compute_agreement_hamming(Y_rep)$diff)
  )
}

#' @keywords internal
.ppc_items_role_clustering <- function(setup, role, random_effect_roles, eligibility_audit) {
  stratum <- sprintf("role:%s", role)
  col_name <- .ppc_role_column_name(role, random_effect_roles)
  if (is.null(col_name) || !(col_name %in% names(setup$event_meta_obs)))
    return(.ppc_cluster_unsupported_items(role, stratum, "role_not_mapped"))

  cluster_id <- setup$event_meta_obs[[col_name]]

  eligible <- .ppc_role_eligible(role, eligibility_audit)
  if (is.na(eligible)) {
    # No audit supplied: fall back to a generic, self-contained structural
    # check (>=20 groups with >=2 observations) so anumaan does not REQUIRE
    # the analysis layer's audit format to function.
    eligible <- .ppc_role_eligible_fallback(cluster_id)
  }
  if (!isTRUE(eligible)) return(.ppc_cluster_unsupported_items(role, stratum, "insufficient_support"))

  restrict_diff_id <- NULL
  if (identical(role, "patient")) {
    admission_col <- .ppc_role_column_name("admission", random_effect_roles)
    if (!is.null(admission_col) && admission_col %in% names(setup$event_meta_obs)) {
      restrict_diff_id <- setup$event_meta_obs[[admission_col]]
    }
  }

  .ppc_cluster_pair_stats(setup, role, cluster_id, restrict_diff_id = restrict_diff_id)
}

# --- spec assembly + main entry point ----------------------------------------

#' @keywords internal
.ppc_build_spec <- function(setup, statistics, random_effect_roles,
                             random_effect_eligibility_audit,
                             min_hospital_support, min_complete_profile_events) {
  items <- list()
  if ("marginal" %in% statistics)
    items <- c(items, .ppc_items_marginal(setup))
  if ("resistant_count" %in% statistics)
    items <- c(items, .ppc_items_resistant_count(setup, min_complete_profile_events))
  if ("pairwise" %in% statistics)
    items <- c(items, .ppc_items_pairwise(setup))
  if ("complete_profile" %in% statistics)
    items <- c(items, .ppc_items_complete_profile(setup, min_complete_profile_events))
  if ("hospital_heterogeneity" %in% statistics)
    items <- c(items, .ppc_items_hospital_heterogeneity(setup, min_hospital_support))
  if ("admission_clustering" %in% statistics)
    items <- c(items, .ppc_items_role_clustering(setup, "admission", random_effect_roles,
                                                  random_effect_eligibility_audit))
  if ("patient_clustering" %in% statistics)
    items <- c(items, .ppc_items_role_clustering(setup, "patient", random_effect_roles,
                                                  random_effect_eligibility_audit))
  items
}

#' Compute AMR-specific posterior predictive discrepancy statistics
#'
#' Drives \code{ppc_draws$generate_state(s)} once per posterior state (see
#' \code{\link{simulate_probit_posterior_predictive}}), streaming each
#' state's complete and observation-masked replicate through every requested
#' discrepancy statistic and accumulating only the compact per-state scalar
#' contributions -- never an \code{S x N_events x D} array. Statistic
#' families (Stan User's Guide-style posterior predictive checks, adapted to
#' AMR resistance data):
#' \describe{
#'   \item{marginal}{Per-class resistance rate, pooled and per eligible
#'     hospital x pathogen x class cell (reusing
#'     \code{fitted_model$eligibility_report$marginal}).}
#'   \item{resistant_count}{Number of resistant classes per event, among
#'     whatever was tested for that event, summarised (mean/variance/median/
#'     p90/max); plus, for fully-observed eligible panels, the full count
#'     distribution.}
#'   \item{pairwise}{RR/RS/SR/SS co-resistance proportions for every
#'     sufficiently co-tested class pair (reusing
#'     \code{fitted_model$eligibility_report$pairwise}).}
#'   \item{complete_profile}{All-susceptible/all-resistant/most-common-profile
#'     frequency, distinct-profile count, Shannon entropy, Simpson
#'     concentration, and the top-3 most frequent OBSERVED profiles
#'     (\code{profile_top_observed_frequency} -- a purely data-dominant
#'     descriptive summary, deliberately NOT a claim of clinical importance;
#'     see the file header of \code{probit_posterior_predictive.R} for the
#'     distinction and the deferred \code{clinical_priority_profiles}
#'     extension point), restricted to adequately-supported fully-observed
#'     hospital x pathogen panels (\code{\link{enumerate_binary_profiles}}-
#'     style class panel via \code{.resolve_profile_class_panel()}).}
#'   \item{hospital_heterogeneity}{SD/IQR/MAD/range of per-hospital
#'     resistance proportions among hospitals meeting
#'     \code{min_hospital_support}, per class.}
#'   \item{admission_clustering, patient_clustering}{Within-/between-group
#'     same-class agreement and Hamming distance, only when the relevant role
#'     is declared eligible (see \code{random_effect_roles} and
#'     \code{random_effect_eligibility_audit} below); skipped cleanly
#'     (\code{support_status != "supported"}), never as a failure, when
#'     unsupported.}
#' }
#'
#' @param fitted_model_or_ppc_draws Either a fitted model object (as for
#'   \code{\link{simulate_probit_posterior_predictive}}) or an
#'   \code{"amr_ppc_draws"} object already returned by that function -- the
#'   latter avoids re-deriving posterior draw setup when statistics are
#'   recomputed for the same simulated replicates.
#' @param n_states,seed,preserve_observation_mask See
#'   \code{\link{simulate_probit_posterior_predictive}}. Ignored when
#'   \code{fitted_model_or_ppc_draws} is already an \code{"amr_ppc_draws"}
#'   object (its own settings are used).
#' @param statistics Character vector of statistic families to compute.
#'   Default: all of \code{"marginal"}, \code{"resistant_count"},
#'   \code{"pairwise"}, \code{"complete_profile"},
#'   \code{"hospital_heterogeneity"}, \code{"admission_clustering"},
#'   \code{"patient_clustering"}.
#' @param random_effect_roles Optional named character vector/list mapping
#'   the name of a GROUPING COLUMN present in
#'   \code{fitted_model$event_metadata} (the NAME of each entry) to one of
#'   \code{"hospital"}, \code{"admission"}, \code{"patient"} (the VALUE -- a
#'   closed set of ROLE labels, never an assumption about what any
#'   particular experiment's columns are actually named), e.g.
#'   \code{c(admission_id = "admission")}. Deliberately decoupled from
#'   whether that column is declared as a random-effect block in this
#'   particular fit -- e.g. tagging a column that is NOT modelled as a
#'   random intercept with role \code{"admission"} is exactly how a fit
#'   lacking an admission random effect gets checked for whether it fails to
#'   reproduce observed within-admission clustering (see Stan User's Guide-
#'   style synthetic recovery tests in
#'   \code{test-probit-predictive-synthetic-recovery.R}). \code{NULL}
#'   (default) means no admission/patient-role clustering statistics can run
#'   (reported with \code{support_status = "role_not_mapped"}).
#' @param random_effect_eligibility_audit Optional data frame with columns
#'   \code{block} (role: \code{"hospital"}/\code{"admission"}/\code{"patient"})
#'   and \code{status} (\code{"eligible_primary"}/\code{"eligible_sensitivity_only"}/
#'   \code{"not_eligible"}) -- e.g. the pathogen-specific random-effect
#'   eligibility audit produced by the analysis layer. When \code{NULL}, a
#'   generic structural fallback check is used instead (see
#'   \code{\link{simulate_probit_posterior_predictive}}'s file header for the
#'   general never-hardcode-eligibility-in-the-package philosophy this
#'   mirrors).
#' @param min_hospital_support Integer; minimum tested events for a hospital
#'   to enter the hospital-heterogeneity comparison. Default 30.
#' @param min_complete_profile_events Integer; minimum fully-observed events
#'   for a hospital x pathogen panel to enter complete-profile / full
#'   resistant-count-distribution statistics. Default 30.
#' @param ci_level Numeric; determines the \code{replicated_q025}/
#'   \code{replicated_q975} quantile levels (\code{(1-ci_level)/2} and its
#'   complement). Default 0.95.
#' @return A tibble with columns \code{statistic_name}, \code{stratum},
#'   \code{observed_value}, \code{replicated_mean}, \code{replicated_sd},
#'   \code{replicated_q025}, \code{replicated_q50}, \code{replicated_q975},
#'   \code{ppc_tail_probability}, \code{ppc_two_sided}, \code{n_replications},
#'   \code{support_status}.
#' @references
#' Stan Development Team. "Posterior and Prior Predictive Checks." Stan
#' User's Guide. \url{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}
#' @export
compute_probit_ppc_statistics <- function(
    fitted_model_or_ppc_draws,
    n_states = 500L,
    seed = 123L,
    preserve_observation_mask = TRUE,
    statistics = c("marginal", "resistant_count", "pairwise", "complete_profile",
                   "hospital_heterogeneity", "admission_clustering", "patient_clustering"),
    random_effect_roles = NULL,
    random_effect_eligibility_audit = NULL,
    min_hospital_support = 30L,
    min_complete_profile_events = 30L,
    ci_level = 0.95
) {
  ppc_draws <- if (inherits(fitted_model_or_ppc_draws, "amr_ppc_draws")) {
    fitted_model_or_ppc_draws
  } else {
    simulate_probit_posterior_predictive(
      fitted_model_or_ppc_draws, n_states = n_states, seed = seed,
      preserve_observation_mask = preserve_observation_mask, return_replicates = FALSE
    )
  }
  setup <- ppc_draws$setup
  S <- ppc_draws$n_states_used

  spec <- .ppc_build_spec(setup, statistics, random_effect_roles,
                           random_effect_eligibility_audit,
                           min_hospital_support, min_complete_profile_events)

  if (length(spec) == 0L) {
    warning("[compute_probit_ppc_statistics] No statistics were requested or computable.", call. = FALSE)
    return(tibble::tibble())
  }

  n_stat <- length(spec)
  T_rep_mat <- matrix(NA_real_, nrow = S, ncol = n_stat)
  needs_state <- vapply(spec, function(it) identical(it$support_status, "supported"), logical(1L))

  for (s in seq_len(S)) {
    st <- ppc_draws$generate_state(s)
    for (i in which(needs_state)) {
      T_rep_mat[s, i] <- spec[[i]]$compute(st$Y_rep, st$Y_rep_complete)
    }
  }

  rows <- vector("list", n_stat)
  for (i in seq_len(n_stat)) {
    it <- spec[[i]]
    rows[[i]] <- .ppc_summarize_statistic(
      it$statistic_name, it$stratum, it$observed_value,
      T_rep_mat[, i], it$support_status, ci_level
    )
  }
  dplyr::bind_rows(rows)
}

# ---------------------------------------------------------------------------
# Faceted PPC plots (Phase A redesign): replace the single generic
# .ppc_dotplot_page() -- which put every hospital x class x pair x profile
# combination on one un-faceted, unparsed `stratum` string axis -- with
# purpose-built plots per statistic family, each parsing `stratum` back into
# real columns (hospital, class, class pair, C, panel_size) and faceting on
# the dimension that has too many levels for one axis. .ppc_dotplot_page()
# itself is UNCHANGED and still used by families not yet redesigned
# (resistant_count_* summary stats, complete-profile statistics, cluster
# statistics) -- see the phased plan this was built under.
# ---------------------------------------------------------------------------

#' @keywords internal
.ppc_common_theme <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 7)
    )
}

#' @keywords internal
.ppc_forest_by_hospital <- function(df, facet_var, title, subtitle, x_label,
                                    facet_ncol = NULL, caption = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df$hospital_display <- factor(hospital_display_label(df$hospital),
                                levels = sort(unique(hospital_display_label(df$hospital)), decreasing = TRUE))
  ggplot2::ggplot(df, ggplot2::aes(y = .data$hospital_display)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$replicated_q025, xmax = .data$replicated_q975),
                           height = 0, colour = "#4292C6", alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$replicated_mean), colour = "#4292C6", size = 1.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed_value), colour = "red", shape = 4,
                        size = 2, stroke = 1) +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_var)), ncol = facet_ncol) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = "Hospital", caption = caption) +
    .ppc_common_theme() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8, face = "italic"))
}

#' Marginal resistance PPC, faceted by class instead of one hospital x class
#' forest plot -- keeps the same statistic (marginal_resistance,
#' hospital-level rows only; the pooled class-level rows are a separate,
#' coarser stratum not shown here) and observed/replicated convention (red
#' cross = observed, blue point + 95\% interval = replicated).
#' @keywords internal
.ppc_plot_marginal <- function(st, title_base) {
  df <- st[st$statistic_name == "marginal_resistance" &
             grepl("^hospital:", st$stratum) &
             st$support_status == "supported", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  m <- regmatches(df$stratum, regexec("^hospital:(.+)\\|class:(.+)$", df$stratum))
  df$hospital <- vapply(m, `[`, character(1L), 2L)
  df$class <- vapply(m, `[`, character(1L), 3L)
  df$class_display <- class_display_label(df$class)
  .ppc_forest_by_hospital(
    df, "class_display",
    paste(title_base, "-- Posterior Predictive Check: Marginal Resistance by Hospital"),
    "Red cross = observed resistance proportion | blue point + bar = posterior-predictive mean and 95% interval",
    "Resistance proportion"
  )
}

#' Resistant-class-count distribution (P(C=c)) for fully-observed panels,
#' faceted by hospital (panel size shown in the facet label) so C=4 for a
#' 4-class panel is never plotted alongside C=5 for a 5-class panel as if
#' they were the same quantity.
#' @keywords internal
.ppc_plot_resistant_count_distribution <- function(st, title_base) {
  df <- st[st$statistic_name == "resistant_count_proportion" &
             st$support_status == "supported", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  m <- regmatches(df$stratum, regexec("^hospital:(.+)\\|pathogen:(.+)\\|panel_size:(\\d+)\\|C=(\\d+)$", df$stratum))
  df$hospital <- vapply(m, `[`, character(1L), 2L)
  df$panel_size <- as.integer(vapply(m, `[`, character(1L), 4L))
  df$C <- as.integer(vapply(m, `[`, character(1L), 5L))
  df$facet_label <- sprintf("%s (panel size = %d)", hospital_display_label(df$hospital), df$panel_size)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$C)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$replicated_q025, ymax = .data$replicated_q975),
                           width = 0.15, colour = "#4292C6", alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$replicated_mean), colour = "#4292C6", size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed_value), colour = "red", shape = 4,
                        size = 2.2, stroke = 1.1) +
    ggplot2::facet_wrap(~facet_label) +
    ggplot2::scale_x_continuous(breaks = scales_int_breaks) +
    ggplot2::labs(
      title = paste(title_base, "-- Posterior Predictive Check: Resistant-Class Count Distribution"),
      subtitle = paste("P(C = c) among fully-observed-panel events, per hospital.",
                       "Red cross = observed | blue point + bar = posterior-predictive mean and 95% interval."),
      x = "Number of resistant classes (C)", y = "Proportion of fully-observed events"
    ) +
    .ppc_common_theme()
}
#' @keywords internal
scales_int_breaks <- function(x) unique(round(pretty(x)))

#' Pairwise co-resistance PPC. Returns a NAMED LIST of four plots (RR, RS,
#' SR, SS) instead of one page -- the previous single forest plot mixed all
#' four cell types with hospital x class-pair labels on one axis, which was
#' the single worst page in the report. Each returned plot facets by class
#' pair (short labels) with hospital on the y-axis.
#' @keywords internal
.ppc_plot_pairwise <- function(st, title_base) {
  df <- st[grepl("^pairwise_", st$statistic_name) & st$support_status == "supported", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  m <- regmatches(df$stratum, regexec("^hospital:(.+)\\|(.+)_x_(.+)$", df$stratum))
  df$hospital <- vapply(m, `[`, character(1L), 2L)
  df$class_1 <- vapply(m, `[`, character(1L), 3L)
  df$class_2 <- vapply(m, `[`, character(1L), 4L)
  df$class_pair <- class_pair_label(df$class_1, df$class_2)
  df$cell <- sub("^pairwise_", "", df$statistic_name)

  cell_labels <- c(RR = "Both resistant (R,R)", RS = "First resistant, second susceptible (R,S)",
                   SR = "First susceptible, second resistant (S,R)", SS = "Both susceptible (S,S)")
  plots <- list()
  for (cell in c("RR", "RS", "SR", "SS")) {
    sub_df <- df[df$cell == cell, , drop = FALSE]
    if (nrow(sub_df) == 0L) next
    plots[[cell]] <- .ppc_forest_by_hospital(
      sub_df, "class_pair",
      paste(title_base, sprintf("-- Posterior Predictive Check: Pairwise Co-resistance (%s)", cell)),
      sprintf("%s | Red cross = observed proportion among co-tested events | blue point + bar = posterior-predictive mean and 95%% interval",
              cell_labels[[cell]]),
      sprintf("Observed/replicated %s proportion", cell)
    )
  }
  plots
}

#' Resistant-class-count summary statistics (Phase B redesign): mean,
#' median, 90th percentile and maximum are counts of resistant classes and
#' share one axis; variance is a dispersion statistic on a different scale
#' (roughly the square of a count) and must not share that axis, even
#' though the old .ppc_dotplot_page() call put all five on one row range
#' labelled "Value". Split into two plots: Panel A (location/tail summaries)
#' and Panel B (variance). There is no hospital/class dimension here (each
#' statistic is a single dataset-wide row, stratum == "all_events"), so
#' neither panel needs faceting.
#' @keywords internal
.ppc_plot_resistant_count_summary <- function(st, title_base) {
  df <- st[grepl("^resistant_count_", st$statistic_name) & st$stratum == "all_events" &
             st$support_status == "supported", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)

  metric_labels_a <- c(resistant_count_mean = "Mean", resistant_count_median = "Median",
                       resistant_count_p90 = "90th percentile", resistant_count_max = "Maximum")

  plots <- list()

  df_a <- df[df$statistic_name %in% names(metric_labels_a), , drop = FALSE]
  if (nrow(df_a) > 0L) {
    df_a$metric_display <- factor(metric_labels_a[df_a$statistic_name],
                                  levels = rev(unname(metric_labels_a)))
    plots$location <- ggplot2::ggplot(df_a, ggplot2::aes(y = .data$metric_display)) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$replicated_q025, xmax = .data$replicated_q975),
                             height = 0, colour = "#4292C6", alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(x = .data$replicated_mean), colour = "#4292C6", size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = .data$observed_value), colour = "red", shape = 4,
                          size = 2.4, stroke = 1.2) +
      ggplot2::expand_limits(x = 0) +
      ggplot2::labs(
        title = paste(title_base, "-- Posterior Predictive Check: Resistant Classes per Event (location)"),
        subtitle = "Red cross = observed | blue point + bar = posterior-predictive mean and 95% interval.",
        x = "Number of resistant classes per event", y = "Summary statistic"
      ) + .ppc_common_theme()
  }

  df_b <- df[df$statistic_name == "resistant_count_variance", , drop = FALSE]
  if (nrow(df_b) > 0L) {
    df_b$metric_display <- "Variance"
    plots$dispersion <- ggplot2::ggplot(df_b, ggplot2::aes(y = .data$metric_display)) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$replicated_q025, xmax = .data$replicated_q975),
                             height = 0, colour = "#4292C6", alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(x = .data$replicated_mean), colour = "#4292C6", size = 2) +
      ggplot2::geom_point(ggplot2::aes(x = .data$observed_value), colour = "red", shape = 4,
                          size = 2.4, stroke = 1.2) +
      ggplot2::expand_limits(x = 0) +
      ggplot2::labs(
        title = paste(title_base, "-- Posterior Predictive Check: Resistant Classes per Event (dispersion)"),
        subtitle = paste("Red cross = observed | blue point + bar = posterior-predictive mean and 95% interval.",
                         "NOT on the same scale as the location/tail-summary page (variance is not a count)."),
        x = "Variance of resistant-class count", y = NULL
      ) + .ppc_common_theme() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }

  plots
}

#' Between-hospital heterogeneity, faceted by spread metric (SD/IQR/MAD/
#' range) instead of concatenating class + metric into one axis. These four
#' metrics are not numerically interchangeable, but they share the same
#' underlying units (spread of a resistance PROPORTION), so faceting them
#' side by side (rather than requiring four separate pages) is appropriate
#' here -- unlike the complete-profile page, where probability, count, and
#' entropy are genuinely different units.
#' @keywords internal
.ppc_plot_hospital_heterogeneity <- function(st, title_base) {
  df <- st[grepl("^hospital_heterogeneity_", st$statistic_name) &
             st$support_status == "supported", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  m <- regmatches(df$stratum, regexec("^class:(.+)$", df$stratum))
  df$class <- vapply(m, `[`, character(1L), 2L)
  df$class_display <- class_display_label(df$class)
  metric_labels <- c(sd = "Standard deviation", iqr = "Interquartile range",
                     mad = "Median absolute deviation", range = "Range")
  df$metric <- sub("^hospital_heterogeneity_", "", df$statistic_name)
  df$metric_display <- factor(metric_labels[df$metric], levels = unname(metric_labels))
  df$class_display <- factor(df$class_display, levels = sort(unique(df$class_display), decreasing = TRUE))

  ggplot2::ggplot(df, ggplot2::aes(y = .data$class_display)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$replicated_q025, xmax = .data$replicated_q975),
                           height = 0, colour = "#4292C6", alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$replicated_mean), colour = "#4292C6", size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$observed_value), colour = "red", shape = 4,
                        size = 2.2, stroke = 1.1) +
    ggplot2::facet_wrap(~metric_display, scales = "free_x") +
    ggplot2::labs(
      title = paste(title_base, "-- Posterior Predictive Check: Between-Hospital Heterogeneity"),
      subtitle = paste("Observed vs. replicated spread of hospital-specific resistance proportions, by class.",
                       "Red cross = observed | blue point + bar = posterior-predictive mean and 95% interval."),
      x = "Spread of hospital-specific resistance proportion", y = "Antimicrobial class"
    ) +
    .ppc_common_theme()
}

#' Complete-profile statistics, split by measurement unit (Phase A does the
#' minimum split needed so the axis is never labelled "Value" across
#' incompatible units; a fuller redesign -- worst-deviation labelling,
#' rare-vs-common profile zoom -- is a later phase). Group A (0-1 scale:
#' frequencies + Simpson concentration) is faceted by hospital; group B
#' (profile count) and group C (entropy) are compact dotplots, since there
#' is only one statistic each per hospital, not a hospital x metric grid.
#' @keywords internal
.ppc_plot_complete_profile_statistics <- function(st, title_base) {
  df <- st[grepl("^profile_", st$statistic_name) &
             st$support_status == "supported" &
             st$statistic_name != "profile_top_observed_frequency", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  m <- regmatches(df$stratum, regexec("^hospital:(.+)\\|pathogen:(.+)$", df$stratum))
  df$hospital <- vapply(m, `[`, character(1L), 2L)

  prob_names <- c("profile_all_susceptible_frequency", "profile_all_resistant_frequency",
                  "profile_most_common_frequency", "profile_simpson_concentration")
  prob_labels <- c(profile_all_susceptible_frequency = "All modelled classes susceptible",
                   profile_all_resistant_frequency = "All modelled classes resistant",
                   profile_most_common_frequency = "Most-common observed profile",
                   profile_simpson_concentration = "Simpson concentration")

  plots <- list()

  df_a <- df[df$statistic_name %in% prob_names, , drop = FALSE]
  if (nrow(df_a) > 0L) {
    df_a$metric_display <- factor(prob_labels[df_a$statistic_name], levels = unname(prob_labels))
    plots$probability <- .ppc_forest_by_hospital(
      df_a, "metric_display",
      paste(title_base, "-- Posterior Predictive Check: Complete-Profile Frequencies"),
      "Fully-observed-panel events only. Red cross = observed | blue point + bar = posterior-predictive mean and 95% interval.",
      "Proportion / concentration",
      caption = paste("\"All modelled classes resistant\" is resistance to only the classes in THIS fit's panel,",
                      "not necessarily every clinically relevant class (e.g. Colistin/Reserve-category agents",
                      "may be excluded from the modelled panel) -- do not read this as \"pan-resistant\".")
    )
  }

  df_b <- df[df$statistic_name == "profile_n_distinct", , drop = FALSE]
  if (nrow(df_b) > 0L) {
    df_b$hospital_display <- factor(hospital_display_label(df_b$hospital),
                                    levels = sort(unique(hospital_display_label(df_b$hospital)), decreasing = TRUE))
    plots$n_distinct <- ggplot2::ggplot(df_b, ggplot2::aes(y = .data$hospital_display)) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$replicated_q025, xmax = .data$replicated_q975),
                             height = 0, colour = "#4292C6", alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(x = .data$replicated_mean), colour = "#4292C6", size = 1.8) +
      ggplot2::geom_point(ggplot2::aes(x = .data$observed_value), colour = "red", shape = 4,
                          size = 2.2, stroke = 1.1) +
      ggplot2::labs(
        title = paste(title_base, "-- Posterior Predictive Check: Number of Distinct Resistance Profiles"),
        subtitle = "Count of distinct complete profiles observed vs. replicated, per hospital -- NOT on the same scale as the frequency/entropy pages.",
        x = "Number of distinct resistance profiles", y = "Hospital"
      ) + .ppc_common_theme()
  }

  df_c <- df[df$statistic_name == "profile_shannon_entropy", , drop = FALSE]
  if (nrow(df_c) > 0L) {
    df_c$hospital_display <- factor(hospital_display_label(df_c$hospital),
                                    levels = sort(unique(hospital_display_label(df_c$hospital)), decreasing = TRUE))
    plots$entropy <- ggplot2::ggplot(df_c, ggplot2::aes(y = .data$hospital_display)) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$replicated_q025, xmax = .data$replicated_q975),
                             height = 0, colour = "#4292C6", alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(x = .data$replicated_mean), colour = "#4292C6", size = 1.8) +
      ggplot2::geom_point(ggplot2::aes(x = .data$observed_value), colour = "red", shape = 4,
                          size = 2.2, stroke = 1.1) +
      ggplot2::labs(
        title = paste(title_base, "-- Posterior Predictive Check: Complete-Profile Entropy"),
        subtitle = "Shannon entropy of the observed complete-profile distribution vs. replicated -- NOT on the same scale as the frequency/count pages.",
        x = "Profile entropy (nats)", y = "Hospital"
      ) + .ppc_common_theme()
  }

  plots
}

# ---------------------------------------------------------------------------
# plot_probit_posterior_predictive_checks(): compact multi-page PDF
# ---------------------------------------------------------------------------

#' @keywords internal
.ppc_dotplot_page <- function(df, title, subtitle, value_label = "Value") {
  df <- df[!is.na(df$observed_value) & df$support_status == "supported", , drop = FALSE]
  if (nrow(df) == 0L) return(NULL)
  df <- df[order(df$stratum), , drop = FALSE]
  df$stratum <- factor(df$stratum, levels = unique(df$stratum))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$stratum)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$replicated_q025, ymax = .data$replicated_q975),
                            width = 0.2, colour = "#4292C6") +
    ggplot2::geom_point(ggplot2::aes(y = .data$replicated_mean), colour = "#4292C6", size = 1.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$observed_value), colour = "red", shape = 4,
                        size = 2.2, stroke = 1.1) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = value_label) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                   axis.text.y = ggplot2::element_text(size = if (nrow(df) > 25L) 6 else 8))
}

#' @keywords internal
.ppc_small_multiples_plot <- function(ppc_replicates, title_base, n_small_multiples, seed) {
  setup <- ppc_replicates$setup
  Y_arr <- ppc_replicates$Y_rep_array   # S x N_ev x D, masked
  S <- dim(Y_arr)[1L]

  count_per_event <- function(Y) {
    mask <- setup$obs_mask
    vapply(seq_len(nrow(Y)), function(i) {
      idx <- which(mask[i, ])
      if (length(idx) == 0L) return(NA_real_)
      sum(Y[i, idx] == 1, na.rm = TRUE)
    }, numeric(1L))
  }

  obs_counts <- count_per_event(setup$obs_ast_mat)
  set.seed(as.integer(seed))
  n_pick <- min(as.integer(n_small_multiples), S)
  picked <- sort(sample.int(S, n_pick))

  rows <- list(data.frame(panel = "Observed", count = obs_counts[!is.na(obs_counts)]))
  for (s in picked) {
    counts_s <- count_per_event(Y_arr[s, , , drop = FALSE][1L, , ])
    rows[[length(rows) + 1L]] <- data.frame(panel = sprintf("Replicate %d", s),
                                             count = counts_s[!is.na(counts_s)])
  }
  df <- do.call(rbind, rows)
  df$panel <- factor(df$panel, levels = c("Observed", sprintf("Replicate %d", picked)))
  df$is_observed <- df$panel == "Observed"

  ggplot2::ggplot(df, ggplot2::aes(x = .data$count, fill = .data$is_observed)) +
    ggplot2::geom_histogram(binwidth = 1, colour = "white") +
    ggplot2::facet_wrap(~ .data$panel, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "#4292C6"), guide = "none") +
    ggplot2::labs(
      title = paste(title_base, "-- Small Multiples: Resistant Classes per Event"),
      subtitle = sprintf(
        "Observed panel (red) vs. %d randomly selected posterior replications (seed=%d)",
        n_pick, seed),
      x = "Number of resistant classes (among tested)", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
}

#' Plot posterior predictive checks for a fitted probit model
#'
#' Compact multi-page PDF: observed-vs-replicated marginal resistance,
#' number-resistant-per-event distribution, observed-vs-replicated pairwise
#' RR/RS/SR/SS, complete-profile concentration/entropy, hospital
#' heterogeneity, admission/patient clustering, and (when
#' \code{ppc_replicates} with materialised replicate arrays is supplied) a
#' small-multiples page comparing the observed data against
#' \code{n_small_multiples} randomly selected posterior replications. A
#' good posterior predictive plot makes it visually apparent whether the
#' observed data look typical among the replicated datasets.
#'
#' @param ppc_statistics Tibble returned by
#'   \code{\link{compute_probit_ppc_statistics}}.
#' @param ppc_replicates Optional \code{"amr_ppc_draws"} object returned by
#'   \code{\link{simulate_probit_posterior_predictive}} with
#'   \code{return_replicates = TRUE} -- enables the small-multiples page.
#'   \code{NULL} (default) skips that page.
#' @param output_pdf_path Path to the PDF to write.
#' @param title_base Character; prefixed to every page's title (e.g.
#'   \code{"<experiment_id>\\n<pathogen>"}).
#' @param n_small_multiples Integer; number of replicate states to show on
#'   the small-multiples page. Default 10.
#' @param small_multiple_seed Integer seed for selecting which replicate
#'   states appear on the small-multiples page. Default 123.
#' @return Invisibly returns the path to the saved PDF, or \code{NULL} if
#'   skipped.
#' @export
plot_probit_posterior_predictive_checks <- function(
    ppc_statistics,
    ppc_replicates = NULL,
    output_pdf_path,
    title_base,
    n_small_multiples = 10L,
    small_multiple_seed = 123L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[plot_probit_posterior_predictive_checks] 'ggplot2' not installed -- skipping.")
    return(invisible(NULL))
  }
  if (is.null(ppc_statistics) || nrow(ppc_statistics) == 0L) {
    message("[plot_probit_posterior_predictive_checks] ppc_statistics is empty -- skipping.")
    return(invisible(NULL))
  }

  dir.create(dirname(output_pdf_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_pdf_path, width = 14, height = 8, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  .try_plot <- function(expr, label) {
    tryCatch(print(expr), error = function(e)
      message(sprintf("  [WARN] %s: %s", label, conditionMessage(e))))
  }

  st <- ppc_statistics

  # 1. Marginal resistance -- faceted by class (Phase A redesign); replaces
  # the old single forest plot that put ~hospital*class rows on one axis.
  p <- .ppc_plot_marginal(st, title_base)
  if (!is.null(p)) .try_plot(p, "marginal resistance")

  # 2. Number resistant per event (summary stats) -- Phase B redesign:
  # split into location/tail summaries (mean/median/p90/max, a count) and
  # dispersion (variance, not a count) -- previously mixed on one "Value"
  # axis.
  rc_plots <- .ppc_plot_resistant_count_summary(st, title_base)
  if (!is.null(rc_plots)) {
    for (nm in names(rc_plots)) .try_plot(rc_plots[[nm]], sprintf("resistant count summary (%s)", nm))
  }
  # 3. Resistant-count distribution (complete panels) -- faceted by hospital
  # (Phase A redesign), panel size shown in the facet label so different
  # panel widths are never plotted as though C meant the same thing.
  p <- .ppc_plot_resistant_count_distribution(st, title_base)
  if (!is.null(p)) .try_plot(p, "resistant count distribution")

  # 4. Pairwise co-resistance -- FOUR separate pages (RR/RS/SR/SS), each
  # faceted by class pair (Phase A redesign); replaces the single worst page
  # in the old report (hundreds of hospital x pair x cell labels on one axis).
  pw_plots <- .ppc_plot_pairwise(st, title_base)
  if (!is.null(pw_plots)) {
    for (cell in names(pw_plots)) .try_plot(pw_plots[[cell]], sprintf("pairwise co-resistance (%s)", cell))
  }

  # 5. Complete-profile statistics -- split by measurement unit (Phase A
  # minimum split: probability/concentration vs. distinct-profile count vs.
  # entropy are never on the same axis; full worst-deviation labelling and
  # rare-profile zoom are a later phase), faceted by hospital within the
  # probability group.
  prof_plots <- .ppc_plot_complete_profile_statistics(st, title_base)
  if (!is.null(prof_plots)) {
    for (nm in names(prof_plots)) .try_plot(prof_plots[[nm]], sprintf("complete profile statistics (%s)", nm))
  }

  # 5c. Hospital heterogeneity -- faceted by spread metric (Phase A
  # redesign); replaces the old class+metric-concatenated single axis.
  p <- .ppc_plot_hospital_heterogeneity(st, title_base)
  if (!is.null(p)) .try_plot(p, "hospital heterogeneity")

  # 5b. Cluster (admission/patient) statistics
  cl <- st[grepl("^cluster_", st$statistic_name), , drop = FALSE]
  if (nrow(cl) > 0L) {
    cl$stratum <- paste(cl$stratum, cl$statistic_name, sep = "|")
    p <- .ppc_dotplot_page(cl, paste(title_base, "-- Admission/Patient Clustering"),
                            "Within-/between-group same-class agreement and Hamming distance",
                            "Value")
    if (!is.null(p)) .try_plot(p, "cluster statistics")
  }

  # 6. Small multiples: observed vs n_small_multiples random replicate states
  if (!is.null(ppc_replicates) && inherits(ppc_replicates, "amr_ppc_draws") &&
      !is.null(ppc_replicates$Y_rep_array)) {
    .try_plot(
      .ppc_small_multiples_plot(ppc_replicates, title_base, n_small_multiples, small_multiple_seed),
      "small multiples"
    )
  }

  message(sprintf("[plot_probit_posterior_predictive_checks] Done -- %s", output_pdf_path))
  invisible(output_pdf_path)
}

# ---------------------------------------------------------------------------
# compute_posterior_predictive_status()
# ---------------------------------------------------------------------------

#' Classify overall posterior predictive fit from a discrepancy-statistics table
#'
#' Does NOT classify from a single tail probability. Groups supported
#' statistics into families (\code{marginal}, \code{pairwise}, \code{profile}
#' [complete-profile and resistant-count statistics -- both address joint
#' multidrug-resistance burden], \code{hospital_heterogeneity}, \code{cluster}
#' [admission and patient clustering]) and flags a family when either (a) any
#' statistic in it is SEVERE (\code{ppc_tail_probability < tail_severe} or
#' \code{> 1 - tail_severe}), or (b) the FRACTION of that family's supported
#' statistics that are extreme exceeds \code{max_fraction_core_extreme}.
#' Thresholds are configurable, not universal hardcoded truths, and are
#' applied as discrepancy FLAGS, not frequentist rejection tests.
#'
#' Statistics with essentially ZERO posterior-predictive variance
#' (\code{replicated_sd} at or near 0 -- e.g. a discrete statistic that
#' lands on the same value in every replicated state, such as
#' \code{max(resistant_count)} pinned at its ceiling for a small class
#' panel) are excluded from extreme/severe classification (though still
#' returned in \code{ppc_statistics} with their true, unmodified
#' \code{ppc_tail_probability}): \code{mean(T_rep >= T_obs)} is exactly 1
#' whenever every replicate exactly EQUALS the observed value, which is a
#' well-known degeneracy of posterior-predictive tail probabilities for
#' near-constant discrete statistics, not evidence of misfit -- a
#' discrepancy statistic that cannot vary under the posterior predictive
#' distribution carries no discriminating information about model fit.
#'
#' @param ppc_statistics Tibble returned by
#'   \code{\link{compute_probit_ppc_statistics}}.
#' @param thresholds Named list overriding any of the defaults:
#'   \code{tail_warning = 0.025}, \code{tail_severe = 0.005},
#'   \code{max_fraction_core_extreme = 0.20}. These are INITIAL AMR-PROJECT
#'   DEFAULTS, not universal statistical truths -- validated only against
#'   this package's own synthetic recovery scenarios
#'   (\code{test-probit-predictive-synthetic-recovery.R}). Expect to
#'   recalibrate them once real multi-pathogen experiment results are
#'   available; override via this argument rather than editing the defaults
#'   in place, so the thresholds actually used for any given run remain
#'   visible in \code{$thresholds_used}.
#' @return List with \code{status} (one of \code{"pass"},
#'   \code{"warning_marginal_ppc"}, \code{"warning_pairwise_ppc"},
#'   \code{"warning_profile_ppc"}, \code{"warning_hospital_heterogeneity_ppc"},
#'   \code{"warning_cluster_ppc"}, \code{"fail_major_ppc_misfit"},
#'   \code{"insufficient_ppc_support"}), \code{reasons} (character vector of
#'   every triggered family-level flag), \code{thresholds_used}, and
#'   \code{family_status} (per-family n/n_extreme/n_severe/fraction_extreme).
#' @export
compute_posterior_predictive_status <- function(ppc_statistics, thresholds = list()) {
  th_defaults <- list(
    tail_warning = 0.025,
    tail_severe = 0.005,
    max_fraction_core_extreme = 0.20
  )
  th <- utils::modifyList(th_defaults, thresholds)

  if (is.null(ppc_statistics) || nrow(ppc_statistics) == 0L)
    return(list(status = "insufficient_ppc_support", reasons = character(0),
                thresholds_used = th, family_status = list()))

  supported <- ppc_statistics[!is.na(ppc_statistics$ppc_tail_probability), , drop = FALSE]
  if (nrow(supported) == 0L)
    return(list(status = "insufficient_ppc_support", reasons = character(0),
                thresholds_used = th, family_status = list()))

  # Zero-(near-zero-)variance statistics carry no discriminating information
  # for extreme/severe classification -- see roxygen note above. Excluded
  # from both the flagging masks AND the family denominator, not merely
  # forced to "not extreme", so they don't dilute a family's fraction_extreme.
  informative <- is.na(supported$replicated_sd) | supported$replicated_sd > 1e-8
  supported <- supported[informative, , drop = FALSE]
  if (nrow(supported) == 0L)
    return(list(status = "insufficient_ppc_support", reasons = character(0),
                thresholds_used = th, family_status = list()))

  family_of_one <- function(nm) {
    if (grepl("^marginal_", nm)) "marginal"
    else if (grepl("^pairwise_", nm)) "pairwise"
    else if (grepl("^profile_", nm) || grepl("^resistant_count_", nm)) "profile"
    else if (grepl("^hospital_heterogeneity_", nm)) "hospital_heterogeneity"
    else if (grepl("^cluster_", nm)) "cluster"
    else "other"
  }
  fam <- vapply(supported$statistic_name, family_of_one, character(1L))

  is_extreme <- supported$ppc_tail_probability < th$tail_warning |
                supported$ppc_tail_probability > (1 - th$tail_warning)
  is_severe  <- supported$ppc_tail_probability < th$tail_severe |
                supported$ppc_tail_probability > (1 - th$tail_severe)

  reasons <- character(0)
  family_status <- list()
  for (f in unique(fam)) {
    idx <- fam == f
    n_f <- sum(idx)
    n_extreme_f <- sum(is_extreme[idx])
    n_severe_f  <- sum(is_severe[idx])
    frac_extreme_f <- n_extreme_f / n_f
    fam_status <- "ok"
    if (n_severe_f > 0L) {
      fam_status <- sprintf("warning_%s_ppc", f)
      reasons <- c(reasons, "fail_major_ppc_misfit", fam_status)
    } else if (frac_extreme_f > th$max_fraction_core_extreme) {
      fam_status <- sprintf("warning_%s_ppc", f)
      reasons <- c(reasons, fam_status)
    }
    family_status[[f]] <- list(n = n_f, n_extreme = n_extreme_f, n_severe = n_severe_f,
                                fraction_extreme = frac_extreme_f, status = fam_status)
  }

  reasons <- unique(reasons)
  status <- if (any(grepl("^fail_", reasons))) {
    "fail_major_ppc_misfit"
  } else if (any(grepl("^warning_", reasons))) {
    reasons[grepl("^warning_", reasons)][1]
  } else {
    "pass"
  }

  list(status = status, reasons = reasons, thresholds_used = th, family_status = family_status)
}
