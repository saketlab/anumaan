# Pathway 2 validation, correlation summaries and diagnostic reports.
# One small correlated fit is shared by every test in this file.

.bvr_cmdstan_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    isTRUE(tryCatch(
      {
        cmdstanr::cmdstan_path()
        TRUE
      },
      error = function(e) FALSE
    ))
}

.bvr_fit <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) {
      return(cache)
    }
    set.seed(7)
    n <- 300L
    shared <- stats::rnorm(n)
    wide <- data.frame(
      center_name = sample(c("H1", "H2", "H3"), n, replace = TRUE), pathogen = "bug",
      Age_normalised = stats::rnorm(n),
      classA = stats::rbinom(n, 1, stats::pnorm(0.2 + shared)),
      classB = stats::rbinom(n, 1, stats::pnorm(-0.3 + shared)),
      classC = stats::rbinom(n, 1, stats::pnorm(-0.6 + 0.5 * shared)),
      event_id = paste0("e", seq_len(n))
    )
    wide$classC[sample(n, 40)] <- NA
    cache <<- suppressWarnings(suppressMessages(fit_bayesian_multivariate_probit(
      event_class_data = wide, class_cols = c("classA", "classB", "classC"),
      fixed_effects = "Age_normalised", random_effects = list(),
      pathogen = "bug", pathogen_col = "pathogen", event_id_col = "event_id",
      profile_group_col = "center_name", residual_structure = "correlated",
      prior_config = list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0),
      sampler_config = list(
        chains = 2L, iter_warmup = 150L, iter_sampling = 150L,
        seed = 7L, parallel_chains = 2L
      ),
      show_messages = FALSE
    )))
    cache
  }
})

skip_without_stan <- function() {
  skip_on_cran()
  skip_if_not(.bvr_cmdstan_available(), "cmdstanr/CmdStan not available")
}

test_that("validate_marginal_calibration reports one coherent row per hospital x class", {
  skip_without_stan()
  m <- suppressMessages(validate_marginal_calibration(.bvr_fit(), n_posterior_draws_for_validation = 100L))
  expect_equal(nrow(m), 9L)
  expect_equal(m$n_resistant + m$n_susceptible, m$n_tested)
  expect_equal(m$observed_resistance, m$n_resistant / m$n_tested, tolerance = 1e-5)
  expect_true(all(m$model_resistance_lower <= m$model_resistance_mean &
    m$model_resistance_mean <= m$model_resistance_upper))
  # no hospital random effect, so per-hospital rates scatter around the pooled
  # fit; the fit should still track them closely on average
  expect_lt(mean(m$absolute_error), 0.1)
})

test_that("validate_pairwise_calibration's four cells are proportions that sum to one", {
  skip_without_stan()
  pw <- suppressMessages(validate_pairwise_calibration(.bvr_fit(), n_posterior_draws_for_validation = 100L))
  expect_equal(nrow(pw), 9L) # 3 class pairs x 3 hospitals
  obs <- pw$observed_RR_prop + pw$observed_RS_prop + pw$observed_SR_prop + pw$observed_SS_prop
  expect_equal(obs, rep(1, nrow(pw)), tolerance = 1e-5)
  mod <- pw$model_pairwise_mean + pw$model_RS_mean + pw$model_SR_mean + pw$model_SS_mean
  expect_equal(mod, rep(1, nrow(pw)), tolerance = 1e-3)
})

test_that("mask_and_validate_ast scores held-out cells reproducibly without touching the RNG", {
  skip_without_stan()
  run <- function() {
    suppressMessages(mask_and_validate_ast(.bvr_fit(),
      fraction_to_mask = 0.1,
      n_posterior_draws_for_validation = 100L, min_tested_after_mask = 10L,
      min_resistant_after_mask = 2L, min_susceptible_after_mask = 2L, seed = 11L
    ))
  }
  set.seed(99)
  before <- .Random.seed
  a <- run()
  expect_identical(.Random.seed, before)
  b <- run()
  expect_identical(a$predictions, b$predictions)

  expect_equal(a$summary$n_masked, nrow(a$predictions))
  expect_equal(a$summary$n_resistant_masked + a$summary$n_susceptible_masked, a$summary$n_masked)
  expect_true(all(a$predictions$predicted_resistance_probability >= 0 &
    a$predictions$predicted_resistance_probability <= 1))
  expect_true(a$summary$auroc >= 0 && a$summary$auroc <= 1)
  expect_null(a$refit_model)
})

test_that("compute_profile_validation_status escalates on miscalibration", {
  skip_without_stan()
  expect_equal(compute_profile_validation_status()$status, "not_evaluated")
  m <- suppressMessages(validate_marginal_calibration(.bvr_fit(), n_posterior_draws_for_validation = 100L))
  bad <- m
  bad$interval_contains_observed <- FALSE
  bad$absolute_error <- 0.5
  st <- compute_profile_validation_status(marginal_tbl = bad)
  expect_match(st$status, "marginal_miscalibration")
  expect_true("warning_marginal_miscalibration" %in% st$reasons)
})

test_that("summarize_fit_correlation_matrix returns valid pairwise correlations", {
  skip_without_stan()
  cm <- summarize_fit_correlation_matrix(.bvr_fit(), "Omega", c("classA", "classB", "classC"))
  expect_equal(nrow(cm), 3L)
  expect_true(all(cm$correlation_mean >= -1 & cm$correlation_mean <= 1))
  expect_true(all(cm$correlation_lower <= cm$correlation_mean & cm$correlation_mean <= cm$correlation_upper))
  # the simulation shares a latent factor across classes
  expect_true(all(cm$correlation_mean > 0))
})

# The "diagnostic and posterior-predictive reports are written to PDF" test
# that lived here exercised plot_probit_diagnostics(),
# simulate_probit_posterior_predictive(), compute_probit_ppc_statistics(),
# and plot_probit_posterior_predictive_checks() -- all removed along with
# R/plot_probit_diagnostics.R and R/probit_posterior_predictive.R.
