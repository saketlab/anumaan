# profile_sensitivity.R
# Sensitivity contrasts across alternative joint Bernoulli fits. See
# profile_fit.R's header for how this fits into the pipeline.
#
# Directly operationalizes the identifiability critique this pipeline was
# built against: "a prior can produce a finite posterior interval even
# when direct information is weak -- that does not establish data-based
# identification... shared interactions can help a hospital with few
# complete panels, but that hospital's estimate must then be labelled as
# relying on cross-hospital structural assumptions." This function makes
# that labelling a number, not a caveat in a methods section: refit once
# under the normal interaction-prior schedule and once under a much
# tighter one (effectively forcing near-independence across classes), and
# any group whose estimate barely moves between the two is directly
# data-supported; any group whose estimate moves a lot is relying on the
# shared eta terms, i.e. cross-group structural borrowing.


#' Named Sensitivity-Scenario Catalog
#'
#' The 10 named scenarios a scientific interpretation of this family's
#' output should check, as a plain reference table -- NOT something this
#' function runs. Running any of them means fitting the compared models
#' explicitly (one \code{\link{profile_fit_joint_bernoulli}} call per
#' side, under whatever changed assumption the scenario names) and
#' passing the results to \code{\link{profile_sensitivity}}; this
#' function never searches for a narrow interval or picks a model for the
#' caller, matching \code{\link{profile_sensitivity}}'s own standing rule.
#'
#' @return A data frame, one row per scenario: \code{scenario_id},
#'   \code{name}, \code{comparable} (what must stay the same between the
#'   two sides for the contrast to mean anything), \code{tests} (what the
#'   comparison actually tells you).
#' @export
profile_sensitivity_scenarios <- function() {
  data.frame(
    scenario_id = c(
      "prior_scale_by_order", "pair_vs_richer_interactions", "shared_vs_hospital_specific_interactions",
      "all_vs_complete_case", "full_vs_supported_roster", "ast_classification_policy",
      "event_dependence_handling", "selective_testing_scenarios", "panel_size", "covariate_functional_form"
    ),
    name = c(
      "Prior scales by interaction order",
      "Pair-only versus richer interactions",
      "Shared versus selected hospital-specific interaction deviations",
      "ALL versus complete-case training",
      "Full hospital roster versus supported roster",
      "Alternative AST/class-collapse policies",
      "One event per patient or dependence-aware model",
      "Plausible selective-testing scenarios",
      "Four versus five/seven class panel",
      "Covariate functional form and case-mix target"
    ),
    comparable = c(
      "Same model target, cohort and panel",
      "Same target and panel",
      "Same target; explicit alternative model",
      "Same prediction target and estimand",
      "Compare common represented hospital targets",
      "Explicit target/policy changes recorded",
      "Explicit revised sampling target",
      "State the selection/completion perturbation",
      "Marginalize extended predictions to the common panel",
      "Same scientific contrast"
    ),
    tests = c(
      "Prior dependence, especially for sparse hospitals",
      "Need for higher-order structure versus instability",
      "Robustness of the common-interaction assumption",
      "Dependence on incomplete data and selection assumptions",
      "Influence of the roster; no unseen-hospital fixed-effect prediction",
      "Sensitivity to measurement and classification",
      "Sensitivity to repeated-event independence",
      "Nonignorability robustness -- cannot claim missingness corrected without a model",
      "Comparable probabilities without confusing a new estimand with improved precision",
      "Adjustment and extrapolation sensitivity"
    ),
    stringsAsFactors = FALSE
  )
}


#' Sensitivity Contrast Across Alternative Model Fits
#'
#' Compares \code{\link{profile_estimate}} output across two or more
#' already-completed fits of the SAME design (same panel, same
#' covariates, same \code{Y}/\code{ids}) made under different modelling
#' assumptions -- different interaction-prior scales (\code{eta_sd} via
#' \code{profile_build_design()}'s \code{interaction_sd}), different
#' \code{max_order}, different missingness handling, etc.
#'
#' This function does NOT run any fits itself, by design, matching this
#' package's standing rule that it never silently chooses a model for the
#' caller: the caller fits each alternative explicitly (one
#' \code{\link{profile_fit_joint_bernoulli}} call per candidate
#' assumption, each its own tracked run with its own sampler diagnostics
#' to check via \code{\link{profile_check_sampler}} before trusting it in
#' this comparison) and passes the resulting \code{fit_results} in.
#'
#' @param fit_results Named list, length \code{>= 2}, each element a
#'   \code{\link{profile_fit_joint_bernoulli}} return value (or raw
#'   \code{CmdStanMCMC}) fitted against the SAME \code{design}/\code{Y}/
#'   \code{ids} under a different assumption. Names become this
#'   comparison's columns (e.g. \code{list(default_priors = fit1,
#'   independence_priors = fit2)}).
#' @param design,Y,ids Same as \code{\link{profile_estimate}} -- the SAME
#'   values for every fit (only the posterior draws differ between fits;
#'   \code{B}/\code{V}/\code{X} are structural and do not depend on the
#'   prior that was used when fitting).
#' @param reference Character or \code{NULL}. Name of the
#'   \code{fit_results} entry everything else is contrasted against.
#'   Default \code{NULL}: the first name in \code{fit_results}.
#' @param estimand Passed to \code{\link{profile_estimate}} for every fit.
#'   Default \code{"observed_preserving_conditional"}.
#' @param level Character, \code{"marginal"} or \code{"complete_profile"}.
#'   \code{"marginal"} collapses \code{profile_estimate()}'s per-state
#'   probabilities to per-class resistance probability first (easier to
#'   read across many classes); \code{"complete_profile"} contrasts the
#'   full per-state probabilities directly.
#' @param n_draws,seed Passed to \code{\link{profile_estimate}} for every
#'   fit.
#'
#' @return A data frame, one row per group x class
#'   (\code{level = "marginal"}) or group x state
#'   (\code{level = "complete_profile"}): each fit's value (one column per
#'   name in \code{fit_results}), one \code{<name>_minus_<reference>}
#'   SIGNED difference column per non-reference fit (direction matters --
#'   a model that moves every group's estimate up under the alternative
#'   assumption is a different finding from one that moves some up and
#'   some down by the same average magnitude, which \code{sensitivity_contrast}
#'   alone cannot distinguish), and \code{sensitivity_contrast} -- the
#'   maximum ABSOLUTE difference from the reference across the other fits.
#'   This single number answers "how much does this group's estimate
#'   move" but not "in which direction" -- use the signed columns for that.
#' @export
profile_sensitivity <- function(fit_results, design, Y, ids,
                                reference = NULL,
                                estimand = "observed_preserving_conditional",
                                level = c("marginal", "complete_profile"),
                                n_draws = 500L,
                                seed = 123L) {
  level <- match.arg(level)
  if (length(fit_results) < 2L || is.null(names(fit_results)) || any(names(fit_results) == "")) {
    stop("fit_results must be a named list of length >= 2 (one entry per alternative assumption).")
  }
  if (is.null(reference)) reference <- names(fit_results)[1]
  if (!reference %in% names(fit_results)) {
    stop(sprintf("reference '%s' not found in names(fit_results)", reference))
  }

  Y <- as.matrix(Y)
  groups <- unique(ids$group)
  fit_names <- names(fit_results)

  estimates <- lapply(fit_results, function(fr) {
    profile_estimate(fr, design, Y, ids, estimand = estimand, n_draws = n_draws, seed = seed)
  })

  if (level == "marginal") {
    D <- ncol(Y)
    classes <- colnames(Y)
    if (is.null(classes)) classes <- paste0("class", seq_len(D))
    B <- design$stan_data$B

    marg_by_fit <- lapply(estimates, function(est) {
      pp <- est$profile_probability # N x S
      vapply(seq_len(D), function(d) rowSums(pp[, B[, d] == 1, drop = FALSE]), numeric(nrow(pp)))
    })

    other_names <- fit_names[fit_names != reference]
    rows <- list()
    for (g in groups) {
      g_idx <- which(ids$group == g)
      for (d in seq_len(D)) {
        vals <- vapply(fit_names, function(nm) mean(marg_by_fit[[nm]][g_idx, d]), numeric(1))
        other_vals <- vals[other_names]
        row <- c(list(group = g, class = classes[d]), as.list(vals))
        for (nm in other_names) row[[paste0(nm, "_minus_", reference)]] <- vals[[nm]] - vals[[reference]]
        row$sensitivity_contrast <- max(abs(other_vals - vals[[reference]]))
        rows[[length(rows) + 1]] <- row
      }
    }
  } else {
    state_labels <- estimates[[1]]$state_labels
    S <- length(state_labels)
    pp_by_fit <- lapply(estimates, function(est) est$profile_probability)

    other_names <- fit_names[fit_names != reference]
    rows <- list()
    for (g in groups) {
      g_idx <- which(ids$group == g)
      for (s in seq_len(S)) {
        vals <- vapply(fit_names, function(nm) mean(pp_by_fit[[nm]][g_idx, s]), numeric(1))
        other_vals <- vals[other_names]
        row <- c(list(group = g, state = state_labels[s]), as.list(vals))
        for (nm in other_names) row[[paste0(nm, "_minus_", reference)]] <- vals[[nm]] - vals[[reference]]
        row$sensitivity_contrast <- max(abs(other_vals - vals[[reference]]))
        rows[[length(rows) + 1]] <- row
      }
    }
  }

  do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
}
