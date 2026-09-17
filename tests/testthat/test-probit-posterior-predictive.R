# Tests for the posterior predictive simulation core
# (R/probit_posterior_predictive.R): .probit_predictive_draws_setup() and
# simulate_probit_posterior_predictive(). These build small, hand-constructed
# "fitted_model"-shaped objects with known posterior states rather than
# fitting real cmdstanr models -- the same no-Stan-fit style used by
# test-random-effects-generic.R -- since these tests exercise the generative
# simulation logic, not the sampler. See test-probit-predictive-synthetic-
# recovery.R for tests that DO fit real (tiny) cmdstanr models.

# ---------------------------------------------------------------------------
# Helper: build a hand-constructed fitted_model-shaped object with known,
# fixed posterior states (S states, each an explicit beta/re_effect/L_Omega).
# ---------------------------------------------------------------------------

.ppc_build_test_fit <- function(X_event, re_data, random_effects, class_cols,
                                 beta_states, re_effect_states,
                                 L_omega_states = NULL,
                                 obs_ast, residual_structure = "identity",
                                 profile_group_col = NULL) {
  S <- length(beta_states)
  K <- ncol(X_event)
  D <- length(class_cols)
  re_prep <- prepare_random_effects(re_data, random_effects)
  total_re_levels <- re_prep$total_re_levels
  N_ev <- nrow(X_event)

  args_list <- list()
  for (a in seq_len(K)) for (b in seq_len(D)) {
    key <- sprintf("beta[%d,%d]", a, b)
    args_list[[key]] <- vapply(seq_len(S), function(s) beta_states[[s]][a, b], numeric(1L))
  }
  for (a in seq_len(D)) for (b in seq_len(total_re_levels)) {
    key <- sprintf("re_effect[%d,%d]", a, b)
    args_list[[key]] <- vapply(seq_len(S), function(s) re_effect_states[[s]][a, b], numeric(1L))
  }
  if (!is.null(L_omega_states)) {
    for (a in seq_len(D)) for (b in seq_len(D)) {
      key <- sprintf("L_Omega[%d,%d]", a, b)
      args_list[[key]] <- vapply(seq_len(S), function(s) L_omega_states[[s]][a, b], numeric(1L))
    }
  }
  draws <- do.call(posterior::draws_array, args_list)

  event_meta <- as.data.frame(re_data)
  event_meta$.event_idx <- seq_len(N_ev)
  event_meta$pathogen <- "bug"
  for (d in seq_len(D)) event_meta[[class_cols[d]]] <- obs_ast[, d]

  list(
    draws = draws,
    class_cols = class_cols,
    event_metadata = event_meta,
    upper_re_col = profile_group_col %||% re_prep$group_cols[1],
    profile_group_col = profile_group_col %||% re_prep$group_cols[1],
    pathogen_col = "pathogen",
    random_effects_prep = re_prep,
    residual_structure = residual_structure,
    index_maps = list(class_levels = class_cols),
    X_event = X_event,
    event_re_idx = re_prep$flat_group_index,
    data_long = data.frame(ev_idx = seq_len(N_ev)),
    eligibility_report = NULL,
    prior_config_used = list(beta_sd = 1.5, tau_sd = 1, lkj_eta = 2)
  )
}

# ---------------------------------------------------------------------------
# Part 18A: PPC is NOT DALY conditional completion
# ---------------------------------------------------------------------------

test_that("simulate_probit_posterior_predictive() replicate is unaffected by the observed AST value, given a fixed theta", {
  re_data <- data.frame(hosp = c("H1", "H1", "H2", "H2"))
  X_event <- matrix(1, nrow = 4, ncol = 1)  # intercept only
  beta_states <- list(matrix(c(0.2, -0.1), nrow = 1, ncol = 2))     # K=1, D=2
  re_effect_states <- list(matrix(0, nrow = 2, ncol = 2))           # RE = 0, 2 hospital levels

  obs_ast_1 <- cbind(c(1, NA, 0, 1), c(NA, NA, NA, NA))
  obs_ast_2 <- obs_ast_1
  obs_ast_2[1, 1] <- 0   # flip event 1's observed class-1 outcome from R to S

  fit1 <- .ppc_build_test_fit(X_event, re_data, "hosp", c("classA", "classB"),
                               beta_states, re_effect_states, obs_ast = obs_ast_1)
  fit2 <- .ppc_build_test_fit(X_event, re_data, "hosp", c("classA", "classB"),
                               beta_states, re_effect_states, obs_ast = obs_ast_2)

  ppc1 <- simulate_probit_posterior_predictive(fit1, n_states = 1L, seed = 7L)
  ppc2 <- simulate_probit_posterior_predictive(fit2, n_states = 1L, seed = 7L)

  st1 <- ppc1$generate_state(1L)
  st2 <- ppc2$generate_state(1L)

  # theta (beta/re_effect) and seed are identical between fit1/fit2 -- only
  # the OBSERVED AST value at event 1/class A differs. The generative
  # replicate must be completely unaffected by that change.
  expect_identical(st1$Y_rep_complete, st2$Y_rep_complete)
})

test_that("simulate_probit_posterior_predictive() never calls the DALY conditional-completion helper", {
  # Static regression guard: the generative code path in this file must never
  # reference .gibbs_conditional_profile_probs(), which conditions on
  # observed AST signs (a fundamentally different, DALY-specific estimand).
  src <- paste(deparse(simulate_probit_posterior_predictive), collapse = "\n")
  expect_false(grepl("gibbs_conditional_profile_probs", src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# Part 18B: identity PPC converges to Phi(mu)
# ---------------------------------------------------------------------------

test_that("identity-residual replicate marginal frequency converges to Phi(mu)", {
  n_ev <- 3000L
  re_data <- data.frame(hosp = rep("H1", n_ev))
  X_event <- matrix(1, nrow = n_ev, ncol = 1)
  mu_val <- 0.4
  n_states <- 50L
  beta_states <- rep(list(matrix(mu_val, nrow = 1, ncol = 1)), n_states)
  re_effect_states <- rep(list(matrix(0, nrow = 1, ncol = 1)), n_states)
  obs_ast <- matrix(1, nrow = n_ev, ncol = 1)  # fully observed, irrelevant to generation

  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", "classA",
                              beta_states, re_effect_states, obs_ast = obs_ast)

  ppc <- simulate_probit_posterior_predictive(fit, n_states = n_states, seed = 11L)
  draws_total <- vapply(seq_len(n_states), function(s) mean(ppc$generate_state(s)$Y_rep_complete),
                         numeric(1L))
  observed_freq <- mean(draws_total)
  expect_equal(observed_freq, stats::pnorm(mu_val), tolerance = 0.01)
})

# ---------------------------------------------------------------------------
# Part 18C: correlated PPC reproduces the sign of the residual correlation
# ---------------------------------------------------------------------------

test_that("correlated-residual replicate concordance increases with positive rho and decreases with negative rho", {
  n_ev <- 4000L
  re_data <- data.frame(hosp = rep("H1", n_ev))
  X_event <- matrix(1, nrow = n_ev, ncol = 1)
  beta_states <- list(matrix(0, nrow = 1, ncol = 2))         # mu = 0 for both classes -> Phi(0) = 0.5
  re_effect_states <- list(matrix(0, nrow = 2, ncol = 1))    # single-level hospital block, RE = 0
  obs_ast <- matrix(1, nrow = n_ev, ncol = 2)

  concordance_for_rho <- function(rho, seed) {
    Omega <- matrix(c(1, rho, rho, 1), nrow = 2)
    L <- t(chol(Omega))
    fit <- .ppc_build_test_fit(X_event, re_data, "hosp", c("classA", "classB"),
                                beta_states, re_effect_states,
                                L_omega_states = list(L), obs_ast = obs_ast,
                                residual_structure = "correlated")
    ppc <- simulate_probit_posterior_predictive(fit, n_states = 1L, seed = seed)
    Y <- ppc$generate_state(1L)$Y_rep_complete
    mean(Y[, 1] == Y[, 2])   # RR + SS fraction (concordance)
  }

  conc_neg  <- concordance_for_rho(-0.6, seed = 1L)
  conc_zero <- concordance_for_rho(0.0,  seed = 1L)
  conc_pos  <- concordance_for_rho(0.6,  seed = 1L)

  expect_lt(conc_neg, conc_zero)
  expect_lt(conc_zero, conc_pos)
})

# ---------------------------------------------------------------------------
# Part 18D: observation mask exactness
# ---------------------------------------------------------------------------

test_that("preserve_observation_mask=TRUE reproduces exactly the real tested-cell pattern", {
  re_data <- data.frame(hosp = c("H1", "H1", "H2", "H2", "H2"))
  X_event <- matrix(1, nrow = 5, ncol = 1)
  beta_states <- list(matrix(c(0.1, -0.2, 0.3), nrow = 1, ncol = 3))
  re_effect_states <- list(matrix(0, nrow = 3, ncol = 2))
  obs_ast <- cbind(
    c(1, NA, 0, 1, NA),
    c(NA, NA, 1, NA, 0),
    c(0, 1, NA, NA, NA)
  )

  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", c("classA", "classB", "classC"),
                              beta_states, re_effect_states, obs_ast = obs_ast)

  ppc <- simulate_probit_posterior_predictive(fit, n_states = 1L, seed = 3L,
                                               preserve_observation_mask = TRUE)
  st <- ppc$generate_state(1L)

  expect_identical(is.na(st$Y_rep), is.na(obs_ast))
  expect_false(any(is.na(st$Y_rep_complete)))
})

test_that("preserve_observation_mask=FALSE returns an unmasked streaming replicate", {
  re_data <- data.frame(hosp = c("H1", "H1", "H2", "H2"))
  X_event <- matrix(1, nrow = 4, ncol = 1)
  beta_states <- list(matrix(c(0.1, -0.2), nrow = 1, ncol = 2))
  re_effect_states <- list(matrix(0, nrow = 2, ncol = 2))
  obs_ast <- cbind(c(1, NA, 0, NA), c(NA, 1, NA, 0))
  fit <- .ppc_build_test_fit(
    X_event, re_data, "hosp", c("classA", "classB"),
    beta_states, re_effect_states, obs_ast = obs_ast
  )

  masked <- simulate_probit_posterior_predictive(
    fit, n_states = 1L, seed = 19L, preserve_observation_mask = TRUE,
    return_replicates = TRUE
  )
  complete <- simulate_probit_posterior_predictive(
    fit, n_states = 1L, seed = 19L, preserve_observation_mask = FALSE,
    return_replicates = TRUE
  )

  expect_identical(is.na(masked$generate_state(1L)$Y_rep), is.na(obs_ast))
  expect_false(anyNA(complete$generate_state(1L)$Y_rep))
  expect_identical(complete$generate_state(1L)$Y_rep,
                   complete$generate_state(1L)$Y_rep_complete)
  expect_identical(complete$Y_rep_array, complete$Y_rep_complete_array)
  expect_identical(masked$generate_state(1L)$Y_rep_complete,
                   complete$generate_state(1L)$Y_rep_complete)
})

# ---------------------------------------------------------------------------
# Part 18E: RE contribution correctness for 1/2/3/arbitrary-named blocks
# ---------------------------------------------------------------------------

test_that("reconstructed mu correctly sums RE contributions for 1, 2, 3, and 4 arbitrarily-named blocks", {
  for (R in 1:4) {
    block_cols <- paste0("blk", seq_len(R))
    re_data <- as.data.frame(setNames(
      lapply(seq_len(R), function(r) rep(c("L1", "L2"), each = 2L)[seq_len(4L)]),
      block_cols
    ))
    random_effects <- lapply(seq_len(R), function(r) {
      list(name = paste0("role_", block_cols[r]), group_col = block_cols[r], terms = "intercept")
    })

    X_event <- matrix(0, nrow = 4, ncol = 1)   # beta contribution deliberately zeroed out below
    beta_states <- list(matrix(0, nrow = 1, ncol = 1))

    re_prep_check <- prepare_random_effects(re_data, random_effects)
    total_levels <- re_prep_check$total_re_levels
    set.seed(100 + R)
    re_effect_mat <- matrix(stats::rnorm(1 * total_levels), nrow = 1, ncol = total_levels)
    re_effect_states <- list(re_effect_mat)

    obs_ast <- matrix(1, nrow = 4, ncol = 1)

    fit <- .ppc_build_test_fit(X_event, re_data, random_effects, "classA",
                                beta_states, re_effect_states, obs_ast = obs_ast)

    setup <- anumaan:::.probit_predictive_draws_setup(fit, n_states = 1L, seed = 1L)
    mu_reconstructed <- setup$mu_all_for_draw(1L)

    mu_expected <- re_contribution(re_effect_mat, fit$event_re_idx)
    expect_equal(mu_reconstructed, mu_expected, tolerance = 1e-12,
                 info = sprintf("R = %d blocks", R))
  }
})

test_that("fixed-only posterior and identity prior predictive use eta = X beta", {
  re_data <- data.frame(profile_group = c("H1", "H1", "H2"))
  X_event <- cbind(1, c(-1, 0, 1))
  beta_state <- matrix(c(0.2, -0.1, 0.3, 0.4), nrow = 2, ncol = 2)
  obs_ast <- matrix(c(1, 0, NA, 0, 1, 1), nrow = 3, ncol = 2)
  fit <- .ppc_build_test_fit(
    X_event, re_data, list(), c("classA", "classB"),
    beta_states = list(beta_state), re_effect_states = list(), obs_ast = obs_ast,
    profile_group_col = "profile_group"
  )

  setup <- anumaan:::.probit_predictive_draws_setup(fit, n_states = 1L, seed = 1L)
  expect_equal(setup$mu_all_for_draw(1L), X_event %*% beta_state)

  prior <- simulate_probit_prior_predictive(fit, n_states = 2L, seed = 1L)
  expect_identical(prior$setup$R, 0L)
  expect_true(all(vapply(seq_len(2L), function(s) {
    all(dim(prior$generate_state(s)$mu) == c(3L, 2L))
  }, logical(1L))))
})

# ---------------------------------------------------------------------------
# Part 18F: reproducibility
# ---------------------------------------------------------------------------

test_that("same seed produces identical replicate arrays across independent calls", {
  n_ev <- 40L
  re_data <- data.frame(hosp = rep(c("H1", "H2"), each = n_ev / 2L))
  X_event <- matrix(1, nrow = n_ev, ncol = 1)
  n_states <- 20L
  set.seed(55)
  beta_states <- lapply(seq_len(n_states), function(s) matrix(stats::rnorm(2), nrow = 1, ncol = 2))
  re_effect_states <- lapply(seq_len(n_states), function(s) matrix(stats::rnorm(4), nrow = 2, ncol = 2))
  obs_ast <- matrix(sample(c(0, 1, NA), n_ev * 2, replace = TRUE), nrow = n_ev, ncol = 2)

  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", c("classA", "classB"),
                              beta_states, re_effect_states, obs_ast = obs_ast)

  ppc1 <- simulate_probit_posterior_predictive(fit, n_states = n_states, seed = 42L,
                                                return_replicates = TRUE)
  ppc2 <- simulate_probit_posterior_predictive(fit, n_states = n_states, seed = 42L,
                                                return_replicates = TRUE)

  expect_identical(ppc1$Y_rep_array, ppc2$Y_rep_array)
  expect_identical(ppc1$Y_rep_complete_array, ppc2$Y_rep_complete_array)
})

test_that("return_replicates = TRUE refuses above the safety element-count ceiling", {
  n_ev <- 2000L
  n_states <- 30000L   # n_ev * n_states * D = 6e7 > 5e7 ceiling
  re_data <- data.frame(hosp = rep("H1", n_ev))
  X_event <- matrix(1, nrow = n_ev, ncol = 1)
  beta_states <- rep(list(matrix(0, nrow = 1, ncol = 1)), n_states)
  re_effect_states <- rep(list(matrix(0, nrow = 1, ncol = 1)), n_states)
  obs_ast <- matrix(1, nrow = n_ev, ncol = 1)
  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", "classA",
                              beta_states, re_effect_states, obs_ast = obs_ast)

  expect_error(
    simulate_probit_posterior_predictive(fit, n_states = n_states, return_replicates = TRUE),
    "safety ceiling"
  )
})

# ---------------------------------------------------------------------------
# Part 18G: no y_rep bloat in the Stan generated-quantities blocks
# ---------------------------------------------------------------------------

test_that("neither generic Stan model variant emits a full y_rep generated quantity", {
  identity_code   <- anumaan:::.amr_probit_stan_generic_identity()
  correlated_code <- anumaan:::.amr_probit_stan_generic_correlated()
  expect_false(grepl("y_rep", identity_code, fixed = TRUE))
  expect_false(grepl("y_rep", correlated_code, fixed = TRUE))
})

test_that("fixed-only Stan variants omit random-effect parameters", {
  identity_code <- anumaan:::.amr_probit_stan_fixed_identity()
  correlated_code <- anumaan:::.amr_probit_stan_fixed_correlated()
  for (code in list(identity_code, correlated_code)) {
    expect_false(grepl("z_re|tau_re|L_corr_re|re_effect|R_block", code))
  }
  expect_false(grepl("L_Omega", identity_code, fixed = TRUE))
  expect_true(grepl("L_Omega", correlated_code, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# replicate_random_effects / random_effect_blocks_to_replicate: API-reserved,
# not implemented here -- must error, never silently ignore.
# ---------------------------------------------------------------------------

test_that("replicate_random_effects = TRUE errors and points to simulate_probit_mixed_predictive()", {
  re_data <- data.frame(hosp = c("H1", "H2"))
  X_event <- matrix(1, nrow = 2, ncol = 1)
  beta_states <- list(matrix(0, nrow = 1, ncol = 1))
  re_effect_states <- list(matrix(0, nrow = 1, ncol = 2))
  obs_ast <- matrix(1, nrow = 2, ncol = 1)
  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", "classA",
                              beta_states, re_effect_states, obs_ast = obs_ast)

  expect_error(
    simulate_probit_posterior_predictive(fit, replicate_random_effects = TRUE),
    "simulate_probit_mixed_predictive"
  )
  expect_error(
    simulate_probit_posterior_predictive(fit, random_effect_blocks_to_replicate = "hosp"),
    "simulate_probit_mixed_predictive"
  )
})
