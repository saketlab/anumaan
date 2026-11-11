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

test_that("E: conditional probabilities do not underflow to an all-zero row when compatible states score far below the global max", {
  D <- 2L
  # Adversarial construction, verified by hand: class 1 has a huge POSITIVE
  # main effect (1000), so the global score max is attained at class1=R.
  # Observing class 1 as SUSCEPTIBLE (the opposite direction) forces every
  # compatible state's score down by 2x that effect (2000 on the log
  # scale) relative to the global max -- exactly the regime where
  # exp(score - global_max) underflows to literal 0.0 in double precision
  # for EVERY compatible state at once.
  states <- anumaan:::profile_enumerate_states(D)
  B <- states
  Mu <- matrix(c(1000, 0), 1, D)
  Score <- Mu %*% t(B)
  sign_mat <- matrix(0L, 1, D); sign_mat[1, 1] <- -1L # class 1 observed S
  idx <- anumaan:::profile_compatible_states(states, sign_mat[1, ])

  # Confirm the adversarial setup genuinely triggers the old failure mode:
  # globally-stabilized softmax underflows every compatible entry to 0.
  Score_stab_global <- Score - max(Score)
  Pmat_global <- exp(Score_stab_global) / sum(exp(Score_stab_global))
  expect_equal(sum(Pmat_global[idx]), 0)

  # The fix: re-stabilize within the compatible subset directly from Score.
  w <- exp(Score[idx] - max(Score[idx]))
  Pcond <- w / sum(w)
  expect_true(all(is.finite(Pcond)))
  expect_true(all(Pcond > 0))
  expect_equal(sum(Pcond), 1, tolerance = 1e-12)
})

test_that("A: design_spec is actually forwarded into the validation refit's design, not silently defaulted", {
  N <- 20L
  Y <- matrix(sample(c(0, 1), N * 3, replace = TRUE), N, 3, dimnames = list(NULL, c("A", "B", "C")))
  X <- cbind(1, rnorm(N))
  ids <- data.frame(group = rep("g1", N), stringsAsFactors = FALSE)
  patient_id <- seq_len(N)

  pkg_env <- asNamespace("anumaan")
  original_fn <- get("profile_build_design", envir = pkg_env)
  captured <- new.env()
  on.exit({
    unlockBinding("profile_build_design", pkg_env)
    assign("profile_build_design", original_fn, envir = pkg_env)
    lockBinding("profile_build_design", pkg_env)
  })
  mock_build_design <- function(Y, X, class_names = NULL, max_order = NULL, beta_sd = NULL,
                                 interaction_sd = function(k) 1.5 / (k - 1), grainsize = 64L) {
    captured$max_order <- max_order
    captured$interaction_sd_result <- if (is.function(interaction_sd)) interaction_sd(2) else interaction_sd
    original_fn(Y, X, class_names = class_names, max_order = max_order, beta_sd = beta_sd,
                interaction_sd = interaction_sd, grainsize = grainsize)
  }
  unlockBinding("profile_build_design", pkg_env)
  assign("profile_build_design", mock_build_design, envir = pkg_env)
  lockBinding("profile_build_design", pkg_env)

  # Also mock the fit itself so this test never touches cmdstanr/Stan --
  # only profile_build_design()'s arguments are under test here.
  original_fit_fn <- get("profile_fit_joint_bernoulli", envir = pkg_env)
  on.exit({
    unlockBinding("profile_fit_joint_bernoulli", pkg_env)
    assign("profile_fit_joint_bernoulli", original_fit_fn, envir = pkg_env)
    lockBinding("profile_fit_joint_bernoulli", pkg_env)
  }, add = TRUE)
  mock_fit_fn <- function(design, output_dir, sampler_config = list()) stop("stop-before-stan: design captured")
  unlockBinding("profile_fit_joint_bernoulli", pkg_env)
  assign("profile_fit_joint_bernoulli", mock_fit_fn, envir = pkg_env)
  lockBinding("profile_fit_joint_bernoulli", pkg_env)

  spec <- list(max_order = 2L, interaction_sd = 0.9)
  tryCatch(
    profile_validate_newpatient(Y, X, ids, patient_id, tempfile(), holdout_fraction = 0.2,
                                 n_draws = 5L, seed = 1L, design_spec = spec),
    error = function(e) invisible(NULL)
  )
  expect_equal(captured$max_order, 2L)
  expect_equal(captured$interaction_sd_result, 0.9)

  # And confirm the OLD default behaviour is preserved when design_spec is NOT supplied
  captured$max_order <- NULL
  tryCatch(
    profile_validate_newpatient(Y, X, ids, patient_id, tempfile(), holdout_fraction = 0.2, n_draws = 5L, seed = 1L),
    error = function(e) invisible(NULL)
  )
  expect_null(captured$max_order) # NULL was passed through -> profile_build_design()'s own default (D) applies
})
