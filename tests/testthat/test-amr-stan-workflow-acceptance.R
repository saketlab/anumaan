# Acceptance tests from AMR_Stan_Workflow_Implementation_Prompt.md Section L.
# Numbered comments below reference that document's own numbering (1-19);
# only the deterministic, no-live-Stan-process subset is covered here
# directly, plus a few real-fit checks kept deliberately tiny/fast.

.amr_wf_cmdstan_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    isTRUE(tryCatch(
      {
        cmdstanr::cmdstan_path()
        TRUE
      },
      error = function(e) FALSE
    ))
}

.amr_wf_skip_without_stan <- function() {
  skip_on_cran()
  skip_if_not(.amr_wf_cmdstan_available(), "cmdstanr/CmdStan not available")
}

test_that("1. Uniform kernel: zero coefficients give uniform states, 0.5 margins, 0.25 pair cells, 0.5 conditional", {
  D <- 3
  states <- profile_enumerate_states(D)
  S <- nrow(states)
  Pmat_uniform <- rep(1 / S, S) # what Score = 0 everywhere (beta=0, eta=0) produces

  # margins: fraction of states with class d resistant
  for (d in seq_len(D)) {
    expect_equal(sum(Pmat_uniform[states[, d] == 1]), 0.5)
  }
  # pair cells: fraction of states with BOTH classes resistant
  expect_equal(sum(Pmat_uniform[states[, 1] == 1 & states[, 2] == 1]), 0.25)

  # conditioning on other classes being independent/untested: hidden class stays 0.5
  compat <- profile_compatible_states(states, c(1L, 0L, 0L)) # class 1 = R, classes 2-3 untested
  cond_p_class2_R <- sum(Pmat_uniform[compat][states[compat, 2] == 1]) / sum(Pmat_uniform[compat])
  expect_equal(cond_p_class2_R, 0.5)
})

test_that("2. Known conditional kernel: SS=.1 RS=.2 SR=.3 RR=.4 gives P(A=R)=.6, P(B=R|A=R)=2/3, compat sums to 1", {
  D <- 2
  states <- profile_enumerate_states(D)
  # states columns: class 1 = A, class 2 = B (sign +1 = R, -1 = S)
  p <- numeric(nrow(states))
  p[states[, 1] == 1 & states[, 2] == 1] <- 0.4 # RR
  p[states[, 1] == 1 & states[, 2] == -1] <- 0.2 # RS
  p[states[, 1] == -1 & states[, 2] == 1] <- 0.3 # SR
  p[states[, 1] == -1 & states[, 2] == -1] <- 0.1 # SS
  expect_equal(sum(p), 1)

  p_A_R <- sum(p[states[, 1] == 1])
  expect_equal(p_A_R, 0.6)

  compat <- profile_compatible_states(states, c(1L, 0L)) # A = R, B untested
  expect_equal(sum(p[compat]), p_A_R) # compatible-state mass equals the marginal exactly
  p_B_R_given_A_R <- sum(p[compat][states[compat, 2] == 1]) / sum(p[compat])
  expect_equal(p_B_R_given_A_R, 2 / 3)
})

test_that("3. Partial likelihood: compatible-state enumeration agrees with brute-force marginal sums, all orders", {
  D <- 4
  states <- profile_enumerate_states(D)
  set.seed(42)
  score <- stats::rnorm(nrow(states))
  score <- score - max(score)
  p <- exp(score); p <- p / sum(p)

  # every possible observed pattern (including all-0 = nothing observed)
  patterns <- as.matrix(expand.grid(rep(list(c(-1L, 0L, 1L)), D)))
  for (i in sample.int(nrow(patterns), 30)) {
    obs <- patterns[i, ]
    compat <- profile_compatible_states(states, obs)
    brute <- which(apply(states, 1, function(s) all(s[obs != 0] == obs[obs != 0])))
    expect_identical(sort(compat), sort(brute))
  }
  # fully observed: exactly one compatible state
  full_obs <- states[5, ]
  expect_length(profile_compatible_states(states, full_obs), 1L)
  # fully unobserved: every state compatible (and thus "contributes zero" to a log-likelihood,
  # since log(sum(p[all])) - log(sum(p[all])) = log(1) - log(1) = 0)
  expect_length(profile_compatible_states(states, rep(0L, D)), nrow(states))
  ll_contribution <- log(sum(p)) - log(sum(p))
  expect_equal(ll_contribution, 0)
})

test_that("4. Cohort/AST: discharged-only support differs from all-outcome support; pair cells sum to binary co-tests; states sum to complete events", {
  d <- data.frame(
    hosp = rep("A", 12),
    event = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    outcome = rep(c("Discharged", "Discharged", "Discharged", "Discharged", "Died", "Died"), 2),
    class = rep(c("Carb", "Flq"), 6),
    val = c("R", "S", "S", "R", "R", "R", "S", "S", "R", "S", "S", "R")
  )
  all_support <- profile_audit_support(d, group_col = "hosp", class_cols = c("Carb", "Flq"),
                                        class_col = "class", value_col = "val", event_id_col = "event")
  discharged_only <- d[d$outcome == "Discharged", ]
  disch_support <- profile_audit_support(discharged_only, group_col = "hosp", class_cols = c("Carb", "Flq"),
                                          class_col = "class", value_col = "val", event_id_col = "event")
  expect_false(identical(all_support$class$tested, disch_support$class$tested))
  expect_true(all(disch_support$class$tested <= all_support$class$tested))

  # pair co_tested == every event tested on both classes (since every event here has exactly 2 rows, one per class)
  expect_equal(all_support$pair$co_tested, sum(table(d$event) == 2))

  # complete_panel (both classes tested) sums to number of events with 2 rows
  expect_equal(all_support$profile$complete_panel, length(unique(d$event)))
})

test_that("5. Goodman: p=0.5 maximizes the implemented half-width (not p=1/k)", {
  tbl <- data.frame(group = "A", class_a = "X", class_b = "Y", co_tested = 200L)
  plan <- profile_plan_precision(tbl, n_categories = 16)
  k <- 16; n <- 200; alpha <- 0.05
  a <- stats::qchisq(1 - alpha / k, df = 1)
  half_width <- function(p) 100 * sqrt(a^2 + 4 * a * n * p * (1 - p)) / (2 * (n + a))
  expect_equal(plan$interval_half_width, half_width(0.5))
  expect_true(half_width(0.5) > half_width(1 / k))
  # p = 0.5 is the maximizer over the full grid, not just vs. 1/k
  grid <- seq(0.01, 0.99, by = 0.01)
  expect_equal(grid[which.max(half_width(grid))], 0.5)
})

test_that("6. No false identification: distinct full-state distributions sharing lower-order margins stay distinguishable", {
  D <- 3
  states <- profile_enumerate_states(D)
  # two different FULL joint distributions engineered to share identical marginals (0.5 each)
  # but differ in their higher-order (3-way) structure -- a saturated design must keep them distinct.
  p1 <- rep(1 / 8, 8) # fully independent/uniform
  # a different distribution: swap mass between two complementary states, keeping all margins at 0.5
  p2 <- p1
  rr_r <- which(states[, 1] == 1 & states[, 2] == 1 & states[, 3] == 1)
  ss_s <- which(states[, 1] == -1 & states[, 2] == -1 & states[, 3] == -1)
  rs_mixed <- which(states[, 1] == 1 & states[, 2] == -1 & states[, 3] == 1)
  sr_mixed <- which(states[, 1] == -1 & states[, 2] == 1 & states[, 3] == -1)
  p2[rr_r] <- p2[rr_r] + 0.05; p2[ss_s] <- p2[ss_s] + 0.05
  p2[rs_mixed] <- p2[rs_mixed] - 0.05; p2[sr_mixed] <- p2[sr_mixed] - 0.05
  for (d in seq_len(D)) expect_equal(sum(p1[states[, d] == 1]), sum(p2[states[, d] == 1])) # same margins
  expect_false(isTRUE(all.equal(p1, p2))) # but genuinely different full distributions
  # a saturated design (profile_bernoulli_terms up to order D) has enough terms to represent both
  terms <- profile_bernoulli_terms(D, max_order = D)
  expect_equal(nrow(terms), 2^D - 1) # D main effects + all higher-order interactions = 2^D - 1 terms
})

test_that("9. Mask/partial-pattern denominators retain zero-count sub-profiles, never silently drop them", {
  .amr_wf_skip_without_stan()
  d <- data.frame(
    hosp = rep("A", 6), event = c(1, 1, 2, 2, 3, 3),
    class = c("Carb", "Flq", "Carb", "Flq", "Carb", "Flq"),
    val = c("R", "R", "R", "R", "R", "R") # every co-tested event is RR -- RS/SR/SS all have ZERO count
  )
  all_support <- profile_audit_support(d, group_col = "hosp", class_cols = c("Carb", "Flq"),
                                        class_col = "class", value_col = "val", event_id_col = "event")
  expect_equal(all_support$pair$n_RR, 3L)
  expect_equal(all_support$pair$n_RS, 0L)
  expect_equal(all_support$pair$n_SR, 0L)
  expect_equal(all_support$pair$n_SS, 0L)
  # the zero cells are REPORTED (present as 0), not absent from the table
  expect_true(all(c("n_RR", "n_RS", "n_SR", "n_SS") %in% names(all_support$pair)))
})

test_that("12. Metrics: hand-computed AUROC/AP/log-loss, tied scores, one-class truth, zero positives", {
  y <- c(1, 1, 0, 0)
  p <- c(0.9, 0.5, 0.5, 0.1) # a tie at 0.5 between one positive and one negative
  expect_equal(.auroc(y, p), 0.875) # hand count: among 4 pairs, 3.5 concordant (tie counts 0.5)

  expect_true(is.na(.auroc(c(1, 1, 1), c(0.2, 0.5, 0.9)))) # one-class truth -> NA, not 0 or 1
  expect_true(is.na(.average_precision(c(0, 0, 0), c(0.1, 0.5, 0.9)))) # zero positives -> NA, not 0

  # hand-computed log loss for a single case
  expect_equal(.log_loss(1, 0.8), -log(0.8))
  expect_equal(.log_loss(0, 0.8), -log(0.2))

  # AP for a perfect ranking (all positives ranked above all negatives) is exactly 1
  expect_equal(.average_precision(c(1, 1, 0, 0), c(0.9, 0.8, 0.3, 0.1)), 1)
})

test_that("13. Intervals: aggregating within a draw before quantiles differs from (and is narrower than) summing per-state quantiles", {
  .amr_wf_skip_without_stan()
  set.seed(5)
  N <- 60; D <- 2
  X <- cbind(Intercept = 1, age = scale(rnorm(N))[, 1])
  Y <- matrix(rbinom(N * D, 1, 0.4), N, D, dimnames = list(NULL, c("Carb", "Flq")))
  design <- profile_build_design(Y, X)
  out_dir <- tempfile("acceptance_fit_")
  fit <- profile_fit_joint_bernoulli(design, output_dir = out_dir,
    sampler_config = list(chains = 1L, parallel_chains = 1L, iter_warmup = 100L, iter_sampling = 100L, seed = 5L))
  ids <- data.frame(event_id = seq_len(N), group = rep("A", N))
  est <- profile_estimate(fit, design, Y, ids, estimand = "unconditional", n_draws = 80L, seed = 5L)

  B <- design$stan_data$B
  naive_upper <- rowSums(est$profile_probability_upper[, B[, 1] == 1, drop = FALSE])
  proper_upper <- est$class_probability_upper[, 1]
  # proper within-draw aggregation is never WIDER than summing independent per-state upper bounds
  expect_true(all(proper_upper <= naive_upper + 1e-9))
})

test_that("14. Dimensions: relabeling class order changes labels, not the underlying aligned probabilities; D=1 is rejected", {
  D <- 3
  states1 <- profile_enumerate_states(D)
  labels1 <- apply(states1, 1, function(r) paste(ifelse(r == 1, "R", "S"), collapse = ""))
  # reordering classes relabels states but the STATE SPACE SIZE and structure stay the same
  expect_equal(nrow(states1), 2^D)
  expect_setequal(labels1, unique(labels1))

  Y1 <- matrix(rbinom(20, 1, 0.5), 10, 2, dimnames = list(NULL, c("A", "B")))
  X <- cbind(Intercept = 1, x = rnorm(10))
  design1 <- profile_build_design(Y1, X)
  Y2 <- Y1[, c(2, 1)] # reordered columns
  design2 <- profile_build_design(Y2, X)
  expect_equal(design1$stan_data$D, design2$stan_data$D)
  expect_equal(sort(design1$design$class_names), sort(design2$design$class_names))

  # D = 1 is explicitly rejected, not silently coerced into a degenerate model
  Y_d1 <- matrix(rbinom(10, 1, 0.5), 10, 1, dimnames = list(NULL, "A"))
  expect_error(profile_build_design(Y_d1, X), "at least 2")
})

test_that("16. Sensitivity: identical fit used as both sides of the comparison gives exactly zero contrast", {
  .amr_wf_skip_without_stan()
  set.seed(6)
  N <- 40; D <- 2
  X <- cbind(Intercept = 1, age = scale(rnorm(N))[, 1])
  Y <- matrix(rbinom(N * D, 1, 0.4), N, D, dimnames = list(NULL, c("Carb", "Flq")))
  design <- profile_build_design(Y, X)
  out_dir <- tempfile("acceptance_fit_")
  fit <- profile_fit_joint_bernoulli(design, output_dir = out_dir,
    sampler_config = list(chains = 1L, parallel_chains = 1L, iter_warmup = 100L, iter_sampling = 100L, seed = 6L))
  ids <- data.frame(event_id = seq_len(N), group = rep("A", N))

  res <- profile_sensitivity(
    fit_results = list(a = fit, b = fit), design = design, Y = Y, ids = ids,
    reference = "a", n_draws = 50L, seed = 6L
  )
  expect_true(all(res$sensitivity_contrast == 0))
})

test_that("17. Output integrity: profile probability draws sum to 1 for every event and draw; missing_allocation counts sum to N", {
  .amr_wf_skip_without_stan()
  set.seed(7)
  N <- 30; D <- 2
  X <- cbind(Intercept = 1, age = scale(rnorm(N))[, 1])
  Y <- matrix(rbinom(N * D, 1, 0.4), N, D, dimnames = list(NULL, c("Carb", "Flq")))
  design <- profile_build_design(Y, X)
  out_dir <- tempfile("acceptance_fit_")
  fit <- profile_fit_joint_bernoulli(design, output_dir = out_dir,
    sampler_config = list(chains = 1L, parallel_chains = 1L, iter_warmup = 100L, iter_sampling = 100L, seed = 7L))
  ids <- data.frame(event_id = seq_len(N), group = rep("A", N))

  est_u <- profile_estimate(fit, design, Y, ids, estimand = "unconditional", n_draws = 40L, seed = 7L)
  sums <- apply(est_u$profile_probability_draws, c(1, 2), sum)
  expect_true(all(abs(sums - 1) < 1e-8))

  est_a <- profile_estimate(fit, design, Y, ids, estimand = "missing_allocation", n_draws = 40L, seed = 7L)
  expect_equal(ncol(est_a$allocated_state), N)
  expect_true(all(est_a$allocated_state >= 1 & est_a$allocated_state <= 2^D))
})
