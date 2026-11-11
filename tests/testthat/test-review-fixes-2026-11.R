# test-review-fixes-2026-11.R
#
# Regression tests locking in five defects confirmed during an external
# source-code review (2026-11): (A) validation refits silently dropping
# the original model spec, (B) average precision not invariant to tied
# score order, (C) masking able to leave an event with zero observed
# classes, (D) sample() misbehaving on a length-1 numeric vector, and
# (E) conditional-probability underflow when compatible states score far
# below the event's global max. Each test reproduces the reviewer's own
# worked example or the exact failure mechanism, verified by hand before
# the fix was written.

test_that("D: stratified new-patient sampling is correct even with a length-1 numeric patient stratum (sample() singleton gotcha)", {
  # sample(42, 1) draws from 1:42, NOT 42 -- this reproduces that exact
  # failure mode: patient_id is numeric, and one hospital has exactly one
  # patient, with a numeric id > 1.
  N <- 10L
  patient_id <- c(42, 7, 7, 8, 8, 9, 9, 10, 10, 11) # patient 42 is a SINGLETON
  strata <- c("H_singleton", rep("H_multi", 9))
  Y <- matrix(sample(c(0, 1, NA), N * 2, replace = TRUE), N, 2)
  Y[, 1] <- 1; Y[, 2] <- 0 # keep it simple/deterministic for class coverage
  X <- cbind(1, rnorm(N))
  ids <- data.frame(group = rep("g1", N), stringsAsFactors = FALSE)

  pkg_env <- asNamespace("anumaan")
  original_fn <- get(".profile_validate_event_holdout", envir = pkg_env)
  on.exit({
    unlockBinding(".profile_validate_event_holdout", pkg_env)
    assign(".profile_validate_event_holdout", original_fn, envir = pkg_env)
    lockBinding(".profile_validate_event_holdout", pkg_env)
  })
  mock_fn <- function(Y, X, ids, output_dir, holdout_idx, sampler_config, n_draws, seed, validation_mode, design_spec = list()) {
    list(holdout_idx = holdout_idx, joint_profile_scores = data.frame(),
         n_holdout_events = length(holdout_idx), n_train_events = nrow(as.matrix(Y)) - length(holdout_idx))
  }
  unlockBinding(".profile_validate_event_holdout", pkg_env)
  assign(".profile_validate_event_holdout", mock_fn, envir = pkg_env)
  lockBinding(".profile_validate_event_holdout", pkg_env)

  for (s in 1:15) {
    res <- profile_validate_newpatient(
      Y, X, ids, patient_id, tempfile(), strata = strata,
      holdout_fraction = 0.5, n_draws = 10L, seed = s
    )
    # The singleton stratum (patient 42 alone) must NEVER be held out --
    # it is kept entirely in training by design (length(ps) < 2 guard).
    # If the sample() bug were present, holdout_patient_ids could contain
    # some arbitrary integer in 1:42 that happens to collide with a real
    # id, or simply never include 42 by construction -- the decisive,
    # unambiguous check is that patient 42 is NEVER selected.
    expect_false(42 %in% res$holdout_patient_ids)
  }
})

test_that("D: pooled (non-stratified) sampling never returns an id outside the supplied distinct_patients, even for numeric singleton-like inputs", {
  patient_id <- c(42, 42, 7, 7, 8, 8) # distinct_patients = c(42, 7, 8); numeric, small values
  N <- length(patient_id)
  Y <- matrix(c(1, 0), N, 2, byrow = TRUE)
  X <- cbind(1, rnorm(N))
  ids <- data.frame(group = rep("g1", N), stringsAsFactors = FALSE)

  pkg_env <- asNamespace("anumaan")
  original_fn <- get(".profile_validate_event_holdout", envir = pkg_env)
  on.exit({
    unlockBinding(".profile_validate_event_holdout", pkg_env)
    assign(".profile_validate_event_holdout", original_fn, envir = pkg_env)
    lockBinding(".profile_validate_event_holdout", pkg_env)
  })
  mock_fn <- function(Y, X, ids, output_dir, holdout_idx, sampler_config, n_draws, seed, validation_mode, design_spec = list()) {
    list(holdout_idx = holdout_idx, joint_profile_scores = data.frame(),
         n_holdout_events = length(holdout_idx), n_train_events = nrow(as.matrix(Y)) - length(holdout_idx))
  }
  unlockBinding(".profile_validate_event_holdout", pkg_env)
  assign(".profile_validate_event_holdout", mock_fn, envir = pkg_env)
  lockBinding(".profile_validate_event_holdout", pkg_env)

  for (s in 1:20) {
    res <- profile_validate_newpatient(Y, X, ids, patient_id, tempfile(), holdout_fraction = 0.34, n_draws = 10L, seed = s)
    expect_true(all(res$holdout_patient_ids %in% c(42, 7, 8)))
  }
})

test_that("C: masking never leaves an event with zero observed classes, even when every one of its classes is independently eligible", {
  # Construct an adversarial case: many events with EXACTLY two observed
  # classes, small eligibility floors, and a high fraction_to_mask -- a
  # setting under which independent per-class sampling is likely to pick
  # BOTH of some event's two observed classes unless the final combined
  # mask is explicitly checked.
  set.seed(42)
  N <- 60L
  Y <- matrix(NA_real_, N, 3, dimnames = list(NULL, c("A", "B", "C")))
  for (i in seq_len(N)) {
    classes <- sample(1:3, 2) # exactly 2 observed classes per event
    Y[i, classes] <- sample(c(0, 1), 2, replace = TRUE)
  }
  ids <- data.frame(group = rep("g1", N), stringsAsFactors = FALSE)

  for (seed in 1:10) {
    masked_df <- anumaan:::.profile_select_mask(
      Y, ids, fraction_to_mask = 0.9, max_mask_per_cell = NULL,
      min_tested_after_mask = 0L, min_resistant_after_mask = 0L, min_susceptible_after_mask = 0L,
      seed = seed
    )
    masked_per_event <- table(masked_df$row)
    n_observed <- rowSums(!is.na(Y))
    for (r in as.integer(names(masked_per_event))) {
      expect_lt(masked_per_event[[as.character(r)]], n_observed[r])
    }
  }
})

test_that("B: average precision is invariant to the order of tied prediction scores (reviewer's exact example)", {
  ap <- anumaan:::.average_precision
  expect_equal(ap(c(1, 0), c(0.5, 0.5)), 0.5)
  expect_equal(ap(c(0, 1), c(0.5, 0.5)), 0.5) # reversed truth order, same ties -- must match
})

test_that("B: average precision matches the hand-verified non-tied worked example (no regression)", {
  ap <- anumaan:::.average_precision
  expect_equal(ap(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.2)), 5 / 6, tolerance = 1e-10)
})

# ---- Shared fixture for production-level underflow tests ----------------
# D = 2, intercept-only X, one event. beta[1,1] = 1000 makes class 1 = R
# dominate every score, so any state with class 1 = S scores ~2000 below
# the global max: exp() over ALL states underflows those to exactly 0.
.underflow_fixture <- function(Y) {
  design <- anumaan::profile_build_design(Y = Y, X = matrix(1, nrow(Y), 1))
  m <- matrix(c(1000, 0, 0), nrow = 4, ncol = 3, byrow = TRUE,
              dimnames = list(NULL, c("beta[1,1]", "beta[1,2]", "eta[1]")))
  fit <- list(draws = function(variables = NULL, ...) posterior::as_draws_matrix(m))
  list(design = design, fit_result = list(fit = fit))
}

test_that("E: profile_estimate() keeps finite, normalised conditional probabilities when compatible states underflow globally", {
  skip_if_not_installed("posterior")
  Y <- matrix(c(0, NA), 1, 2, dimnames = list(NULL, c("A", "B"))) # class A observed S
  fx <- .underflow_fixture(Y)
  est <- profile_estimate(fx$fit_result, fx$design, Y, data.frame(group = "g1"),
                          estimand = "observed_preserving_conditional", n_draws = 4L, seed = 1L)
  p <- est$profile_probability[1, ]
  expect_true(all(is.finite(p)))
  expect_equal(sum(p), 1, tolerance = 1e-12)
  compat <- substr(names(p), 1, 1) == "S" # states consistent with A = S
  expect_equal(sum(p[compat]), 1, tolerance = 1e-12)
  expect_true(all(p[!compat] == 0)) # observed result preserved
  expect_equal(unname(p[compat]), c(0.5, 0.5), tolerance = 1e-12) # beta for B = 0, eta = 0

  alloc <- profile_estimate(fx$fit_result, fx$design, Y, data.frame(group = "g1"),
                            estimand = "missing_allocation", n_draws = 4L, seed = 1L)
  expect_true(all(substr(alloc$state_labels[alloc$allocated_state], 1, 1) == "S"))
})

test_that("E: the masked-cell predictor uses the same stable normalisation (no NA draws silently dropped)", {
  skip_if_not_installed("posterior")
  Y <- matrix(c(0, 1), 1, 2, dimnames = list(NULL, c("A", "B"))) # B is masked below
  fx <- .underflow_fixture(Y)
  pred <- anumaan:::.profile_predict_masked_cells(
    fx$fit_result, fx$design, Y, data.frame(row = 1L, class_idx = 2L), n_draws = 4L, seed = 1L
  )
  expect_true(all(is.finite(pred)))
  expect_equal(as.numeric(pred), 0.5, tolerance = 1e-12) # P(B = R | A = S) with beta_B = 0, eta = 0
  expect_equal(attr(pred, "n_draws_used"), 4L)
})

# ---- Mocks for refit entry points ----------------------------------------
.with_refit_mocks <- function(code) {
  pkg_env <- asNamespace("anumaan")
  orig_build <- get("profile_build_design", envir = pkg_env)
  orig_fit <- get("profile_fit_joint_bernoulli", envir = pkg_env)
  captured <- new.env()
  set_binding <- function(name, value) {
    unlockBinding(name, pkg_env); assign(name, value, envir = pkg_env); lockBinding(name, pkg_env)
  }
  on.exit({ set_binding("profile_build_design", orig_build); set_binding("profile_fit_joint_bernoulli", orig_fit) })
  set_binding("profile_build_design", function(Y, X, ...) {
    d <- orig_build(Y, X, ...)
    captured$design <- d
    d
  })
  set_binding("profile_fit_joint_bernoulli", function(design, output_dir, sampler_config = list()) stop("stop-before-stan"))
  code(captured)
}

.refit_fixture <- function() {
  set.seed(1)
  N <- 30L
  Y <- matrix(sample(c(0, 1), N * 3, replace = TRUE), N, 3, dimnames = list(NULL, c("A", "B", "C")))
  list(Y = Y, X = cbind(1, rnorm(N)), ids = data.frame(group = rep(c("g1", "g2"), length.out = N)), N = N)
}

test_that("A: design_spec reaches the built design at all three refit entry points", {
  f <- .refit_fixture()
  spec <- list(max_order = 2L, interaction_sd = 0.9)
  check <- function(captured) {
    expect_equal(captured$design$design$max_order, 2L)
    expect_equal(unique(captured$design$stan_data$eta_sd), 0.9)
    expect_equal(captured$design$stan_data$P, 3L) # 3 pairs, no three-way term at max_order = 2
  }
  .with_refit_mocks(function(captured) {
    expect_error(profile_validate_newpatient(f$Y, f$X, f$ids, seq_len(f$N), tempfile(),
                                             holdout_fraction = 0.2, n_draws = 5L, design_spec = spec),
                 "stop-before-stan")
    check(captured)
  })
  .with_refit_mocks(function(captured) {
    expect_error(profile_validate_temporal(f$Y, f$X, f$ids, event_date = as.Date("2024-01-01") + seq_len(f$N),
                                           cutoff_date = as.Date("2024-01-20"), output_dir = tempfile(),
                                           n_draws = 5L, design_spec = spec),
                 "stop-before-stan")
    check(captured)
  })
  .with_refit_mocks(function(captured) {
    expect_error(profile_validate_heldout(f$Y, f$X, f$ids, tempfile(), fraction_to_mask = 0.2,
                                          min_tested_after_mask = 0L, min_resistant_after_mask = 0L,
                                          min_susceptible_after_mask = 0L, n_draws = 5L, design_spec = spec),
                 "stop-before-stan")
    check(captured)
  })
})

test_that("A: unknown design_spec fields are rejected, and an empty spec is announced", {
  f <- .refit_fixture()
  .with_refit_mocks(function(captured) {
    expect_error(profile_validate_newpatient(f$Y, f$X, f$ids, seq_len(f$N), tempfile(),
                                             holdout_fraction = 0.2, design_spec = list(max_ordr = 2L)),
                 "Unknown design_spec field")
    expect_message(
      expect_error(profile_validate_newpatient(f$Y, f$X, f$ids, seq_len(f$N), tempfile(), holdout_fraction = 0.2),
                   "stop-before-stan"),
      "design_spec is empty"
    )
  })
})

test_that("D: high holdout fractions never empty a stratum; train and test patients never overlap", {
  pkg_env <- asNamespace("anumaan")
  orig <- get(".profile_validate_event_holdout", envir = pkg_env)
  on.exit({ unlockBinding(".profile_validate_event_holdout", pkg_env)
            assign(".profile_validate_event_holdout", orig, envir = pkg_env)
            lockBinding(".profile_validate_event_holdout", pkg_env) })
  unlockBinding(".profile_validate_event_holdout", pkg_env)
  assign(".profile_validate_event_holdout", function(Y, X, ids, output_dir, holdout_idx, ...) list(holdout_idx = holdout_idx),
         envir = pkg_env)
  lockBinding(".profile_validate_event_holdout", pkg_env)

  patient_id <- c(1, 1, 2, 3, 3, 4, 5, 6, 7, 8)           # patient 1 and 3 have repeat events
  strata <- c("h2", "h2", "h2", "h5", "h5", "h5", "h5", "h5", "h5", "h5") # h2 has exactly 2 patients
  N <- length(patient_id)
  Y <- matrix(1, N, 2); X <- cbind(1, seq_len(N)); ids <- data.frame(group = strata)
  for (s in 1:20) {
    res <- profile_validate_newpatient(Y, X, ids, patient_id, tempfile(), strata = strata,
                                       holdout_fraction = 0.9, seed = s, design_spec = list(max_order = 2L))
    held <- res$holdout_patient_ids
    for (h in unique(strata)) {
      expect_gte(length(setdiff(unique(patient_id[strata == h]), held)), 1L) # >= 1 training patient per stratum
    }
    expect_true(all(patient_id[res$holdout_idx] %in% held))                 # every event of a held patient is held
    expect_false(any(patient_id[-res$holdout_idx] %in% held))               # no train/test patient overlap
    expect_true(is.numeric(held))                                           # numeric IDs not coerced to character
  }
  expect_error(profile_validate_newpatient(Y, X, ids, patient_id, tempfile(), holdout_fraction = 1),
               "strictly between 0 and 1")
})

test_that("F: aggregators expose computational_checks_passed_* and no longer emit eligible_for_*", {
  est <- list(estimand = "missing_allocation", state_labels = c("SS", "SR", "RS", "RR"),
              allocated_state = matrix(c(1L, 3L, 1L, 4L), nrow = 2, byrow = TRUE))
  out <- profile_aggregate_counts_with_allocation(est, data.frame(group = c("g", "g")), "p", "c",
                                                  sampler_check = list(status = "pass"),
                                                  min_n_events = 1L, min_n_draws = 1L)
  expect_true(all(c("computational_checks_passed_YLL", "computational_checks_passed_YLD") %in% names(out)))
  expect_false(any(c("eligible_for_YLL", "eligible_for_YLD") %in% names(out)))
  # Exact allocation fixture: one fixed SS event, one partial event -> SS = 1, RS = 0.5, RR = 0.5, SR = 0
  expect_equal(setNames(out$N_mean, out$profile_delta)[c("SS", "SR", "RS", "RR")], c(SS = 1, SR = 0, RS = 0.5, RR = 0.5))
})
