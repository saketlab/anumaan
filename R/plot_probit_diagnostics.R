# plot_probit_diagnostics.R
# Generates a multi-page diagnostic PDF for one pathogen fit from
# fit_bayesian_multivariate_probit().
#
# Pages produced (each skipped gracefully if data unavailable):
#   1. Sampling trace plots      — key params (beta, tau_*, lp__)
#   2. Warmup + sampling traces  — requires CmdStanMCMC CSV files still on disk
#   3. Rank plots                — chain mixing diagnostic
#   4. Beta posterior densities  — 89% / 95% credible intervals
#   5. NUTS energy (E-BFMI)      — geometry health check
#   6. NUTS step size & acceptance rate
#   7. Rhat histogram            — convergence across all parameters
#   8. ESS histogram             — sampling efficiency across all parameters
#
# Called once per pathogen per experiment, after sampling completes.
# Runtime: ~1-2 minutes. Has no effect on sampling speed.

#' Generate diagnostic plots for a fitted Bayesian multivariate probit model
#'
#' @param fit_obj    Named list returned by \code{fit_bayesian_multivariate_probit()}.
#' @param output_dir Directory where the PDF will be saved.
#' @param experiment_id  Character. Experiment identifier (used in filename and titles).
#' @param pathogen   Character. Pathogen name (used in plot titles).
#' @param max_params Integer. Maximum number of parameters shown in trace/rank plots.
#'   Defaults to 20. Only \code{beta}, \code{tau_*}, and \code{lp__} are ever included.
#'
#' @return Invisibly returns the path to the saved PDF, or \code{NULL} if skipped.
#' @export
plot_probit_diagnostics <- function(
    fit_obj,
    output_dir,
    experiment_id,
    pathogen,
    max_params  = 20L
) {

  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    message("[plot_probit_diagnostics] 'bayesplot' not installed — skipping diagnostic plots.")
    message("[plot_probit_diagnostics] Install with: install.packages('bayesplot')")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[plot_probit_diagnostics] 'ggplot2' not installed — skipping diagnostic plots.")
    return(invisible(NULL))
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    message("[plot_probit_diagnostics] 'posterior' not installed — skipping diagnostic plots.")
    return(invisible(NULL))
  }

  draws    <- fit_obj$draws        # draws_array (sampling only, no z_free)
  stan_fit <- fit_obj$fit          # CmdStanMCMC object
  diag     <- fit_obj$diagnostics

  if (is.null(draws)) {
    message("[plot_probit_diagnostics] fit_obj$draws is NULL — skipping.")
    return(invisible(NULL))
  }

  iter_warmup <- if (!is.null(diag) && !is.null(diag$iter_warmup))
                   as.integer(diag$iter_warmup[[1L]]) else 1000L

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path   <- file.path(output_dir, paste0(experiment_id, "_diagnostics.pdf"))
  title_base <- sprintf("%s\n%s", experiment_id, pathogen)

  # -- Select key parameters (exclude z_free / RE arrays which are huge) ------
  all_vars <- dimnames(draws)[[3]]
  priority <- unique(c(
    grep("^lp__$",            all_vars, value = TRUE),
    grep("^beta\\[",          all_vars, value = TRUE),
    grep("^tau_hospital\\[",  all_vars, value = TRUE),
    grep("^tau_patient\\[",   all_vars, value = TRUE),
    grep("^tau_admission\\[", all_vars, value = TRUE)
  ))
  if (length(priority) > max_params) priority <- priority[seq_len(max_params)]

  if (length(priority) == 0L) {
    message("[plot_probit_diagnostics] No plottable parameters found — skipping.")
    return(invisible(NULL))
  }

  message(sprintf("[plot_probit_diagnostics] Writing diagnostic PDF: %s", pdf_path))

  grDevices::pdf(pdf_path, width = 14, height = 8, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Helper: print plot, catch errors silently
  .try_plot <- function(expr, label) {
    tryCatch(
      print(expr),
      error = function(e)
        message(sprintf("  [WARN] %s: %s", label, conditionMessage(e)))
    )
  }

  n_chains    <- dim(draws)[2]
  n_sampling  <- dim(draws)[1]

  # ── 1. Sampling trace plots ──────────────────────────────────────────────────
  .try_plot(
    bayesplot::mcmc_trace(draws, pars = priority) +
      ggplot2::labs(
        title    = paste(title_base, "— Sampling Traces"),
        subtitle = sprintf("%d chains × %d sampling iterations",
                           n_chains, n_sampling)
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
    "sampling trace"
  )

  # ── 2. Warmup + sampling traces (needs CSV files still on disk) ─────────────
  base_vars    <- unique(sub("\\[.*", "", priority))
  warmup_draws <- tryCatch(
    stan_fit$draws(inc_warmup = TRUE, variables = base_vars, format = "draws_array"),
    error = function(e) {
      message(sprintf(
        "  [INFO] Warmup draws unavailable (CSV files may be cleaned up): %s",
        conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(warmup_draws)) {
    wvars <- intersect(priority, dimnames(warmup_draws)[[3]])
    if (length(wvars) > 0L) {
      .try_plot(
        bayesplot::mcmc_trace(warmup_draws, pars = wvars,
                              n_warmup = iter_warmup) +
          ggplot2::labs(
            title    = paste(title_base, "— Warmup + Sampling Traces"),
            subtitle = sprintf(
              "Vertical line at warmup/sampling boundary (iter = %d)", iter_warmup)
          ) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
        "warmup trace"
      )
    }
  }

  # ── 3. Rank plots ────────────────────────────────────────────────────────────
  .try_plot(
    bayesplot::mcmc_rank_overlay(draws, pars = priority) +
      ggplot2::labs(
        title    = paste(title_base, "— Rank Plots"),
        subtitle = "Uniform distribution across ranks indicates good chain mixing"
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
    "rank plot"
  )

  # ── 4. Beta posterior densities ──────────────────────────────────────────────
  beta_vars <- grep("^beta\\[", priority, value = TRUE)
  if (length(beta_vars) > 0L) {
    .try_plot(
      bayesplot::mcmc_areas(draws, pars = beta_vars,
                            prob = 0.89, prob_outer = 0.95) +
        ggplot2::labs(
          title    = paste(title_base, "— Beta Coefficients"),
          subtitle = "Shaded: 89% CI | outer line: 95% CI"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
      "beta area plot"
    )
  }

  # ── 5 & 6. NUTS diagnostics (energy + acceptance) ───────────────────────────
  np <- tryCatch(bayesplot::nuts_params(stan_fit), error = function(e) NULL)
  if (!is.null(np)) {
    .try_plot(
      bayesplot::mcmc_nuts_energy(np) +
        ggplot2::labs(
          title    = paste(title_base, "— NUTS Energy"),
          subtitle = "E-BFMI < 0.3 suggests poor geometry"
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
      "energy plot"
    )
    .try_plot(
      bayesplot::mcmc_nuts_acceptance(np) +
        ggplot2::labs(
          title    = paste(title_base, "— NUTS Step Size & Acceptance Rate"),
          subtitle = sprintf("adapt_delta target: %.2f",
                             if (!is.null(diag$adapt_delta[[1L]])) diag$adapt_delta[[1L]] else 0.8)
        ) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
      "acceptance plot"
    )
  }

  # ── 7. Rhat histogram ────────────────────────────────────────────────────────
  rhat_df <- tryCatch({
    s <- posterior::summarise_draws(draws, "rhat")
    s[!is.na(s$rhat), , drop = FALSE]
  }, error = function(e) NULL)

  if (!is.null(rhat_df) && nrow(rhat_df) > 0L) {
    max_rhat <- round(max(rhat_df$rhat, na.rm = TRUE), 3)
    .try_plot(
      ggplot2::ggplot(rhat_df, ggplot2::aes(x = rhat)) +
        ggplot2::geom_histogram(bins = 40, fill = "#4292C6", colour = "white") +
        ggplot2::geom_vline(xintercept = 1.01, colour = "red",
                            linetype = "dashed", linewidth = 0.8) +
        ggplot2::annotate("text", x = 1.01, y = Inf,
                          label = " Rhat = 1.01", hjust = 0, vjust = 1.5,
                          colour = "red", size = 3.5) +
        ggplot2::labs(
          title    = paste(title_base, "— Rhat Distribution"),
          subtitle = sprintf("max Rhat = %.3f | %d parameters | red line at 1.01",
                             max_rhat, nrow(rhat_df)),
          x = "Rhat", y = "Count"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)),
      "Rhat histogram"
    )
  }

  # ── 8. ESS histogram ─────────────────────────────────────────────────────────
  ess_df <- tryCatch({
    s <- posterior::summarise_draws(draws, "ess_bulk", "ess_tail")
    s[(!is.na(s$ess_bulk)) | (!is.na(s$ess_tail)), , drop = FALSE]
  }, error = function(e) NULL)

  if (!is.null(ess_df) && nrow(ess_df) > 0L) {
    bulk_rows <- data.frame(ess = ess_df$ess_bulk, metric = "ESS bulk",
                            stringsAsFactors = FALSE)
    tail_rows <- data.frame(ess = ess_df$ess_tail, metric = "ESS tail",
                            stringsAsFactors = FALSE)
    ess_long <- rbind(bulk_rows[!is.na(bulk_rows$ess), ],
                      tail_rows[!is.na(tail_rows$ess), ])
    min_bulk <- round(min(ess_df$ess_bulk, na.rm = TRUE))
    min_tail <- round(min(ess_df$ess_tail, na.rm = TRUE))
    .try_plot(
      ggplot2::ggplot(ess_long, ggplot2::aes(x = ess, fill = metric)) +
        ggplot2::geom_histogram(bins = 40, alpha = 0.75, position = "identity") +
        ggplot2::geom_vline(xintercept = 100, colour = "red",
                            linetype = "dashed", linewidth = 0.8) +
        ggplot2::scale_fill_manual(
          values = c("ESS bulk" = "#2166AC", "ESS tail" = "#74ADD1")) +
        ggplot2::labs(
          title    = paste(title_base, "— ESS Distribution"),
          subtitle = sprintf("min ESS bulk = %d | min ESS tail = %d | red line at 100",
                             min_bulk, min_tail),
          x = "Effective Sample Size", y = "Count", fill = NULL
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(plot.title    = ggplot2::element_text(size = 10),
                       legend.position = "top"),
      "ESS histogram"
    )
  }

  message(sprintf("[plot_probit_diagnostics] Done — %s", pdf_path))
  invisible(pdf_path)
}
