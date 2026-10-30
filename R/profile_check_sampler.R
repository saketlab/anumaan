# profile_check_sampler.R
# Sampler-validity verdict + diagnostic plots for a profile_fit_joint_bernoulli()
# fit. See profile_fit.R's header for how this fits into the pipeline.


#' Check Sampler Validity and Render Diagnostic Plots
#'
#' The numeric half of "did this Stan run actually work": checks the
#' standard HMC/NUTS validity criteria (divergences, max-treedepth hits,
#' E-BFMI, split-Rhat, bulk/tail ESS) against widely used thresholds from
#' the Stan/\code{posterior} documentation, and renders the three standard
#' diagnostic plots (trace, rank, energy) via \code{bayesplot} -- the same
#' \code{mcmc_trace()}/\code{mcmc_rank_overlay()} functions already used
#' elsewhere in \code{anumaan} (\code{plot_worst_offender_diagnostics.R}),
#' plus \code{mcmc_nuts_energy()} for the energy check that file doesn't
#' currently render.
#'
#' This is deliberately upstream of any scientific judgment about the
#' MODEL (calibration, held-out performance, which classes/panels to
#' trust) -- those are \code{profile_ppc}, \code{profile_validate_masked},
#' \code{profile_validate_heldout}, and \code{profile_sensitivity}
#' (elsewhere in this family). A sampler can pass every check here and
#' still fit a badly miscalibrated model; this function only answers
#' whether the SAMPLER did its job, which is a precondition for trusting
#' anything downstream of it.
#'
#' @param fit_result The list returned by
#'   \code{\link{profile_fit_joint_bernoulli}} (reads \code{$fit},
#'   \code{$summary}, \code{$diagnostics}), or a raw \code{cmdstanr}
#'   \code{CmdStanMCMC} object (summary/diagnostics are then computed
#'   fresh from it).
#' @param rhat_max Numeric. Split-Rhat above this value fails a parameter.
#'   Default \code{1.01} (Vehtari et al. 2021's revised threshold).
#' @param ess_min Integer. Bulk or tail ESS below this fails a parameter.
#'   Default \code{400L} (the Stan-documented minimum for stable 95%
#'   interval estimates).
#' @param ebfmi_min Numeric. E-BFMI below this, per chain, fails that
#'   chain. Default \code{0.2} (Betancourt 2017's threshold).
#' @param variables Character vector or \code{NULL}. Which parameters'
#'   Rhat/ESS to check and plot. Default \code{NULL}: every parameter in
#'   the summary except \code{lp__}.
#' @param output_dir Character or \code{NULL}. If supplied, writes
#'   \code{sampler_check.rds} (status/reasons/checks, not the plot
#'   objects) and each plot as its own PDF into this directory.
#' @param ci_level Numeric in (0, 1). Credible-interval level the
#'   quantile-MCSE checks (both \code{mcse} and \code{derived_mcse}) are
#'   computed at. Default \code{0.95} (MCSE of the 0.025 and 0.975
#'   quantiles -- \code{posterior::mcse_quantile()}'s \code{ess_tail()}
#'   alone summarizes a DIFFERENT combination of tails and does not
#'   replace checking the actual interval endpoints you report).
#' @param derived_estimates Named list or \code{NULL}. Each element a
#'   numeric matrix of posterior draws (\code{n_draws x n_quantities}) for
#'   a DERIVED estimand not among \code{variables} -- e.g. per-hospital
#'   resistance probabilities from \code{\link{profile_estimate}}'s
#'   \code{class_probability_draws}/\code{profile_probability_draws},
#'   reshaped to one column per hospital x class/state. When supplied,
#'   mean and quantile MCSE are computed per column and returned as
#'   \code{derived_mcse} -- raw sampler diagnostics on \code{beta}/
#'   \code{eta} do not by themselves say whether a DERIVED probability is
#'   precise enough to report. Default \code{NULL} (skipped). NOTE: unlike
#'   \code{mcse} (computed on the fit's own chain-aware draws),
#'   \code{derived_estimates} columns are treated as one pooled sequence
#'   -- \code{\link{profile_estimate}} does not track which chain each
#'   subsampled draw came from, so no between-chain correction is applied
#'   here; this MCSE is a reasonable, slightly optimistic approximation,
#'   not the exact chain-aware quantity \code{mcse} reports.
#' @param mcse_budget Numeric or \code{NULL}. The MCSE budget each
#'   \code{derived_estimates} column's mean and quantile MCSE is checked
#'   against -- e.g. \code{0.1 * epsilon}, where \code{epsilon} is the
#'   approved probability half-width target (a project convention, not
#'   something this function derives). \code{NULL} (default): MCSE is
#'   still reported, just not pass/failed against anything. Extending
#'   sampling addresses a numerically UNSTABLE estimate (MCSE too wide for
#'   the budget); it does NOT make a long but genuinely imprecise
#'   scientific interval (wide because of real data/model uncertainty)
#'   narrower -- these are different problems, and this check only
#'   diagnoses the former.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{status}}{\code{"pass"}, \code{"warning"}, or
#'       \code{"fail"}. \code{"warning"} means only ESS (precision, not
#'       sampler validity) failed to clear \code{ess_min}; any divergence,
#'       treedepth, E-BFMI, or Rhat failure is \code{"fail"}. MCSE-budget
#'       failures (on \code{derived_estimates}) do NOT affect
#'       \code{status} -- they are a separate, reporting-level judgment
#'       the caller makes explicitly, not a sampler-validity verdict.}
#'     \item{\code{reasons}}{character vector of failed checks' names
#'       (empty when \code{status == "pass"})}
#'     \item{\code{checks}}{data frame, one row per check: \code{check},
#'       \code{value}, \code{threshold}, \code{status}}
#'     \item{\code{mcse}}{data frame, one row per \code{variables} entry:
#'       \code{variable}, \code{mcse_mean}, \code{mcse_lower_quantile},
#'       \code{mcse_upper_quantile} (at \code{ci_level}'s endpoints)}
#'     \item{\code{derived_mcse}}{empty (zero rows) unless
#'       \code{derived_estimates} is supplied; one row per
#'       (estimate name, column): \code{estimate_name}, \code{column},
#'       \code{mcse_mean}, \code{mcse_lower_quantile},
#'       \code{mcse_upper_quantile}, \code{budget}, \code{passes_budget}
#'       (\code{NA} when \code{mcse_budget} is \code{NULL})}
#'     \item{\code{plots}}{list of 3 \code{ggplot} objects: \code{trace},
#'       \code{rank}, \code{energy}}
#'   }
#' @export
profile_check_sampler <- function(fit_result,
                                  rhat_max = 1.01,
                                  ess_min = 400L,
                                  ebfmi_min = 0.2,
                                  variables = NULL,
                                  output_dir = NULL,
                                  ci_level = 0.95,
                                  derived_estimates = NULL,
                                  mcse_budget = NULL) {
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("Package 'bayesplot' is required for profile_check_sampler().", call. = FALSE)
  }

  fit <- if (is.list(fit_result) && "fit" %in% names(fit_result)) fit_result$fit else fit_result
  if (!inherits(fit, "CmdStanMCMC")) {
    stop("fit_result must be the list from profile_fit_joint_bernoulli() or a CmdStanMCMC object.")
  }

  summary_tbl <- if (is.list(fit_result) && "summary" %in% names(fit_result)) {
    fit_result$summary
  } else {
    fit$summary()
  }
  diag <- if (is.list(fit_result) && "diagnostics" %in% names(fit_result)) {
    fit_result$diagnostics
  } else {
    fit$diagnostic_summary()
  }

  if (is.null(variables)) {
    variables <- setdiff(summary_tbl$variable, "lp__")
  }
  sub_summary <- summary_tbl[summary_tbl$variable %in% variables, , drop = FALSE]

  checks <- list()
  add_check <- function(name, value, threshold, ok) {
    checks[[length(checks) + 1]] <<- data.frame(
      check = name, value = value, threshold = threshold,
      status = if (ok) "pass" else "fail", stringsAsFactors = FALSE
    )
  }

  n_divergent <- sum(diag$num_divergent)
  add_check("divergent_transitions", n_divergent, 0, n_divergent == 0)

  n_treedepth <- sum(diag$num_max_treedepth)
  add_check("max_treedepth_hits", n_treedepth, 0, n_treedepth == 0)

  min_ebfmi <- min(diag$ebfmi, na.rm = TRUE)
  add_check("min_ebfmi", round(min_ebfmi, 4), ebfmi_min, min_ebfmi >= ebfmi_min)

  max_rhat <- max(sub_summary$rhat, na.rm = TRUE)
  add_check("max_rhat", round(max_rhat, 4), rhat_max, max_rhat <= rhat_max)

  min_ess_bulk <- min(sub_summary$ess_bulk, na.rm = TRUE)
  add_check("min_ess_bulk", round(min_ess_bulk, 1), ess_min, min_ess_bulk >= ess_min)

  min_ess_tail <- min(sub_summary$ess_tail, na.rm = TRUE)
  add_check("min_ess_tail", round(min_ess_tail, 1), ess_min, min_ess_tail >= ess_min)

  checks_tbl <- do.call(rbind, checks)
  reasons <- checks_tbl$check[checks_tbl$status == "fail"]
  status <- if (length(reasons) == 0L) {
    "pass"
  } else if (all(reasons %in% c("min_ess_bulk", "min_ess_tail"))) {
    "warning" # low ESS alone is a precision problem, not a sampler-validity failure
  } else {
    "fail"
  }

  draws <- fit$draws(variables = variables, format = "draws_array")
  plots <- list(
    trace = bayesplot::mcmc_trace(draws) + ggplot2::ggtitle("Trace"),
    rank = bayesplot::mcmc_rank_overlay(draws) + ggplot2::ggtitle("Rank"),
    energy = bayesplot::mcmc_nuts_energy(bayesplot::nuts_params(fit)) + ggplot2::ggtitle("Energy")
  )

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for profile_check_sampler()'s MCSE checks.", call. = FALSE)
  }
  lo_q <- (1 - ci_level) / 2
  hi_q <- 1 - lo_q
  mcse_tbl <- do.call(rbind, lapply(variables, function(v) {
    x <- posterior::subset_draws(draws, variable = v) # keeps the iteration x chain structure
    data.frame(
      variable = v,
      mcse_mean = posterior::mcse_mean(x),
      mcse_lower_quantile = posterior::mcse_quantile(x, probs = lo_q),
      mcse_upper_quantile = posterior::mcse_quantile(x, probs = hi_q),
      stringsAsFactors = FALSE
    )
  }))
  rownames(mcse_tbl) <- NULL

  # ---- derived_estimates: MCSE for quantities not among `variables` ----
  derived_rows <- list()
  if (!is.null(derived_estimates)) {
    for (nm in names(derived_estimates)) {
      mat <- derived_estimates[[nm]]
      for (j in seq_len(ncol(mat))) {
        x <- mat[, j]
        mcse_m <- posterior::mcse_mean(x)
        mcse_lo <- posterior::mcse_quantile(x, probs = lo_q)
        mcse_hi <- posterior::mcse_quantile(x, probs = hi_q)
        worst_mcse <- max(mcse_m, mcse_lo, mcse_hi, na.rm = TRUE)
        derived_rows[[length(derived_rows) + 1]] <- data.frame(
          estimate_name = nm, column = j,
          mcse_mean = mcse_m, mcse_lower_quantile = mcse_lo, mcse_upper_quantile = mcse_hi,
          budget = if (is.null(mcse_budget)) NA_real_ else mcse_budget,
          passes_budget = if (is.null(mcse_budget)) NA else worst_mcse <= mcse_budget,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  derived_mcse_tbl <- if (length(derived_rows) > 0L) {
    out <- do.call(rbind, derived_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(estimate_name = character(0), column = integer(0), mcse_mean = numeric(0),
               mcse_lower_quantile = numeric(0), mcse_upper_quantile = numeric(0),
               budget = numeric(0), passes_budget = logical(0), stringsAsFactors = FALSE)
  }

  result <- list(
    status = status, reasons = reasons, checks = checks_tbl,
    mcse = mcse_tbl, derived_mcse = derived_mcse_tbl, plots = plots
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result[c("status", "reasons", "checks", "mcse", "derived_mcse")], file.path(output_dir, "sampler_check.rds"))
    for (nm in names(plots)) {
      grDevices::pdf(file.path(output_dir, sprintf("diagnostic_%s.pdf", nm)), width = 10, height = 7)
      print(plots[[nm]])
      grDevices::dev.off()
    }
  }

  result
}
