# Tests for R/plot_omega_diagnostics.R's degeneracy/summary-table additions
# and R/plot_worst_offender_diagnostics.R.

.omega_draws_2x2 <- function(rho, n_draws = 50L, n_chains = 2L) {
  arr <- array(0,
    dim = c(n_draws, n_chains, 4),
    dimnames = list(NULL, NULL, c("Omega[1,1]", "Omega[2,1]", "Omega[1,2]", "Omega[2,2]"))
  )
  arr[, , "Omega[1,1]"] <- 1
  arr[, , "Omega[2,2]"] <- 1
  arr[, , "Omega[2,1]"] <- rho
  arr[, , "Omega[1,2]"] <- rho
  posterior::as_draws_array(arr)
}

test_that(".omega_degeneracy_stats() computes the analytic smallest eigenvalue for a 2x2 correlation matrix (1 - rho)", {
  fit <- list(residual_structure = "correlated", draws = .omega_draws_2x2(0.6), class_cols = c("A", "B"))
  stats_out <- .omega_degeneracy_stats(fit, c("A", "B"))
  expect_false(is.null(stats_out))
  expect_true(all(abs(stats_out$lambda_min - 0.4) < 1e-6))
})

test_that(".omega_degeneracy_stats()/plot_omega_degeneracy_diagnostic() return NULL for identity fits or missing Omega", {
  identity_fit <- list(residual_structure = "identity", draws = .omega_draws_2x2(0.6), class_cols = c("A", "B"))
  expect_null(.omega_degeneracy_stats(identity_fit, c("A", "B")))
  expect_null(plot_omega_degeneracy_diagnostic(identity_fit, c("A", "B")))

  no_draws_fit <- list(residual_structure = "correlated", draws = NULL, class_cols = c("A", "B"))
  expect_null(.omega_degeneracy_stats(no_draws_fit, c("A", "B")))
})

test_that("plot_omega_degeneracy_diagnostic() flags near-degeneracy at a threshold above 1 - rho, not below it", {
  fit <- list(residual_structure = "correlated", draws = .omega_draws_2x2(0.9), class_cols = c("A", "B"))
  p <- plot_omega_degeneracy_diagnostic(fit, c("A", "B"), degenerate_threshold = 0.5)
  expect_false(is.null(p))
  expect_true(grepl("100\\.0% of draws", p$labels$subtitle))
  p_low <- plot_omega_degeneracy_diagnostic(fit, c("A", "B"), degenerate_threshold = 0.05)
  expect_true(grepl("0\\.0% of draws", p_low$labels$subtitle))
})

test_that("plot_omega_summary_table() returns NULL for empty input and includes ESS tail when present", {
  expect_null(plot_omega_summary_table(NULL, c("A", "B")))
  corr_summary <- data.frame(
    class_1 = "Aminoglycosides", class_2 = "Carbapenems",
    correlation_mean = 0.8, correlation_median = 0.8,
    correlation_lower = 0.7, correlation_upper = 0.9,
    rhat = 1.02, ess_bulk = 300, ess_tail = 500,
    stringsAsFactors = FALSE
  )
  p <- plot_omega_summary_table(corr_summary, c("Aminoglycosides", "Carbapenems"), "test")
  expect_false(is.null(p))
  lines <- ggplot2::ggplot_build(p)$data[[1]]$label
  expect_true(any(grepl("500", lines)))
})

.fake_worst_offender_fit <- (function() {
  vars <- c("beta[1,1]", "beta[2,1]", "Omega[1,2]", "lp__")
  arr <- array(rnorm(200 * length(vars)),
    dim = c(50, 4, length(vars)),
    dimnames = list(NULL, NULL, vars)
  )
  draws <- posterior::as_draws_array(arr)
  monitored <- data.frame(
    variable = vars,
    parameter_group = c("beta", "beta", "Omega", "lp__"),
    rhat = c(1.08, 1.00, 1.02, 1.001),
    ess_bulk = c(60, 900, 150, 950),
    ess_tail = c(150, 1200, 400, 1300),
    stringsAsFactors = FALSE
  )
  X <- matrix(0, nrow = 2, ncol = 2, dimnames = list(NULL, c("(Intercept)", "Age_normalised")))
  list(
    draws = draws,
    X_event = X,
    class_cols = c("Aminoglycosides", "Carbapenems"),
    diagnostics_detail = list(monitored_parameters = monitored)
  )
})()

test_that("plot_probit_worst_offender_diagnostics() relabels the selected variables and returns trace+rank plots", {
  out <- plot_probit_worst_offender_diagnostics(.fake_worst_offender_fit, .fake_worst_offender_fit$draws, "test", n_worst = 2L)
  expect_false(is.null(out))
  expect_true(all(c("trace", "rank") %in% names(out)))
  expect_false(is.null(out$trace))
  expect_false(is.null(out$rank))
})

test_that("plot_probit_worst_offender_diagnostics() returns NULL when monitored_parameters is unavailable", {
  broken_fit <- .fake_worst_offender_fit
  broken_fit$diagnostics_detail <- list()
  expect_null(plot_probit_worst_offender_diagnostics(broken_fit, broken_fit$draws, "test"))
})

# The four ".probit_interpretation_text() ..." tests that lived here exercised
# a function removed along with R/plot_probit_diagnostics.R. Its replacement,
# .probit_classify_fitting_issues() (R/plot_diagnostic_recommendations.R), has
# a different signature and return shape (list(primary=, secondary=,
# not_indicated=, bucket=) instead of a flat "PRIMARY ISSUE:"/"SECONDARY
# ISSUE:" character vector) and does not appear to replicate one piece of the
# old behavior -- naming the specific worst-converging Omega pair (e.g.
# "AMG x FQ") when one pair's Rhat exceeds the threshold. Rather than guess at
# a mechanical port of tests for a feature that may have been deliberately
# dropped, new tests against .probit_classify_fitting_issues()'s actual
# current behavior belong in their own follow-up, written against that
# function directly.
