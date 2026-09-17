# Tests for compute_probit_ppc_statistics(), plot_probit_posterior_predictive_checks(),
# and compute_posterior_predictive_status() (R/probit_posterior_predictive.R,
# second half). Two tiers:
#   1. Fast, hand-constructed-fit unit tests (no cmdstanr required) covering
#      output schema, the status/threshold logic, and admission/patient role
#      mapping in isolation.
#   2. Part 19's five synthetic-recovery scenarios: real (small/fast)
#      cmdstanr fits where the TRUE generative process is known, proving
#      each PPC statistic family actually detects the specific model
#      deficiency it is designed to detect (not just that the pipeline runs).

# ---------------------------------------------------------------------------
# Tier 1: fast unit tests (reuses the .ppc_build_test_fit() helper defined in
# test-probit-posterior-predictive.R, sourced into the same environment when
# the full suite runs; redefined here for standalone-file execution).
# ---------------------------------------------------------------------------

if (!exists(".ppc_build_test_fit")) {
  .ppc_build_test_fit <- function(X_event, re_data, random_effects, class_cols,
                                   beta_states, re_effect_states,
                                   L_omega_states = NULL,
                                   obs_ast, residual_structure = "identity") {
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
      upper_re_col = re_prep$group_cols[1],
      pathogen_col = "pathogen",
      random_effects_prep = re_prep,
      residual_structure = residual_structure,
      index_maps = list(class_levels = class_cols),
      X_event = X_event,
      event_re_idx = re_prep$flat_group_index,
      data_long = data.frame(ev_idx = seq_len(N_ev)),
      eligibility_report = NULL
    )
  }
}

test_that("compute_probit_ppc_statistics() returns the exact Part 7 output schema", {
  re_data <- data.frame(hosp = rep(c("H1", "H2"), each = 30))
  X_event <- matrix(1, nrow = 60, ncol = 1)
  S <- 20L
  set.seed(1)
  beta_states <- lapply(seq_len(S), function(s) matrix(rnorm(2, sd = 0.3), nrow = 1, ncol = 2))
  re_effect_states <- lapply(seq_len(S), function(s) matrix(rnorm(4, sd = 0.2), nrow = 2, ncol = 2))
  obs_ast <- cbind(rbinom(60, 1, 0.4), rbinom(60, 1, 0.5))
  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", c("classA", "classB"),
                              beta_states, re_effect_states, obs_ast = obs_ast)

  stats_tbl <- compute_probit_ppc_statistics(fit, n_states = S, seed = 1,
                                              statistics = c("marginal", "resistant_count"))
  expect_true(is.data.frame(stats_tbl))
  expect_identical(names(stats_tbl), c(
    "statistic_name", "stratum", "observed_value", "replicated_mean", "replicated_sd",
    "replicated_q025", "replicated_q50", "replicated_q975",
    "ppc_tail_probability", "ppc_two_sided", "n_replications", "support_status"
  ))
  expect_false(any(names(stats_tbl) == "p_value"))
  expect_true(all(stats_tbl$ppc_two_sided[!is.na(stats_tbl$ppc_two_sided)] <= 1))
})

test_that("finite-replicate PPC tails cannot be zero", {
  low <- anumaan:::.ppc_summarize_statistic(
    "synthetic", "all", -10, rep(0, 9), "supported", 0.95
  )
  high <- anumaan:::.ppc_summarize_statistic(
    "synthetic", "all", 10, rep(0, 9), "supported", 0.95
  )

  expect_equal(low$ppc_tail_probability, 1)
  expect_equal(low$ppc_two_sided, 0.2)
  expect_equal(high$ppc_tail_probability, 0.1)
  expect_equal(high$ppc_two_sided, 0.2)
  expect_gt(low$ppc_two_sided, 0)
  expect_gt(high$ppc_tail_probability, 0)
})

test_that("compute_posterior_predictive_status() family-level logic", {
  base <- tibble::tibble(
    statistic_name = c("marginal_resistance", "pairwise_RR", "profile_shannon_entropy",
                        "hospital_heterogeneity_sd", "cluster_within_group_same_class_agreement"),
    stratum = "s", observed_value = 0.5,
    replicated_mean = 0.5, replicated_sd = 0.1,
    replicated_q025 = 0.3, replicated_q50 = 0.5, replicated_q975 = 0.7,
    ppc_tail_probability = 0.5, ppc_two_sided = 1,
    n_replications = 100L, support_status = "supported"
  )

  clean <- compute_posterior_predictive_status(base)
  expect_equal(clean$status, "pass")

  severe <- base
  severe$ppc_tail_probability[severe$statistic_name == "pairwise_RR"] <- 0.001
  res <- compute_posterior_predictive_status(severe)
  expect_equal(res$status, "fail_major_ppc_misfit")
  expect_true(any(grepl("pairwise", res$reasons)))

  extreme_frac <- dplyr::bind_rows(
    base[base$statistic_name == "profile_shannon_entropy", ],
    base[base$statistic_name == "profile_shannon_entropy", ],
    base[base$statistic_name == "profile_shannon_entropy", ]
  )
  extreme_frac$stratum <- paste0("s", seq_len(3))
  extreme_frac$ppc_tail_probability <- c(0.5, 0.01, 0.01)
  res2 <- compute_posterior_predictive_status(extreme_frac)
  expect_equal(res2$status, "warning_profile_ppc")

  res3 <- compute_posterior_predictive_status(tibble::tibble())
  expect_equal(res3$status, "insufficient_ppc_support")

  th_res <- compute_posterior_predictive_status(severe, thresholds = list(tail_severe = 0.0001))
  expect_false(identical(th_res$status, "fail_major_ppc_misfit"))
})

test_that("admission-role clustering statistic detects strong within-admission agreement when correctly mapped", {
  set.seed(11)
  n_adm <- 60L
  ev_per_adm <- 3L
  adm_effect <- rep(rnorm(n_adm, sd = 2), each = ev_per_adm)   # strong per-admission latent tendency
  N_ev <- n_adm * ev_per_adm
  re_data <- data.frame(
    hosp = rep("H1", N_ev),
    adm = rep(paste0("A", seq_len(n_adm)), each = ev_per_adm)
  )
  X_event <- matrix(1, nrow = N_ev, ncol = 1)
  mu_true <- adm_effect   # class A driven entirely by admission-level tendency
  obs_ast <- cbind(classA = rbinom(N_ev, 1, stats::pnorm(mu_true)))

  # Fit declares ONLY "hosp" as a random-effect block -- admission is present
  # as a raw column but never modelled.
  S <- 1L
  beta_states <- list(matrix(0, nrow = 1, ncol = 1))
  re_prep_check <- prepare_random_effects(re_data, "hosp")
  re_effect_states <- list(matrix(0, nrow = 1, ncol = re_prep_check$total_re_levels))

  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", "classA",
                              beta_states, re_effect_states, obs_ast = obs_ast)

  stats_tbl <- compute_probit_ppc_statistics(
    fit, n_states = 1L, seed = 1, statistics = "admission_clustering",
    random_effect_roles = c(adm = "admission")
  )
  agree_row <- stats_tbl[stats_tbl$statistic_name == "cluster_within_group_same_class_agreement", ]
  expect_equal(agree_row$support_status, "supported")
  expect_gt(agree_row$observed_value, 0.7)   # strong true within-admission agreement
})

test_that("clustering statistics report role_not_mapped / insufficient_support cleanly, never as failures", {
  re_data <- data.frame(hosp = c("H1", "H2"))
  X_event <- matrix(1, nrow = 2, ncol = 1)
  beta_states <- list(matrix(0, nrow = 1, ncol = 1))
  re_effect_states <- list(matrix(0, nrow = 1, ncol = 2))
  obs_ast <- matrix(1, nrow = 2, ncol = 1)
  fit <- .ppc_build_test_fit(X_event, re_data, "hosp", "classA",
                              beta_states, re_effect_states, obs_ast = obs_ast)

  no_role <- compute_probit_ppc_statistics(fit, n_states = 1L, statistics = "admission_clustering")
  expect_true(all(no_role$support_status == "role_not_mapped"))
  expect_true(all(is.na(no_role$observed_value)))

  bad_col <- compute_probit_ppc_statistics(fit, n_states = 1L, statistics = "admission_clustering",
                                            random_effect_roles = c(nonexistent_col = "admission"))
  expect_true(all(bad_col$support_status == "role_not_mapped"))
})

# ---------------------------------------------------------------------------
# Tier 2: Part 19 synthetic recovery scenarios (real, small/fast cmdstanr fits)
# ---------------------------------------------------------------------------

.ppcrec_cmdstan_available <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({ cmdstanr::cmdstan_path(); TRUE }, error = function(e) FALSE))
}

.ppcrec_fit <- function(wide, class_cols, random_effects, residual_structure = "identity",
                         fixed_effects = "Age_normalised", prior_config = list(),
                         chains = 2L, iter = 200L, seed = 1L) {
  suppressWarnings(fit_bayesian_multivariate_probit(
    event_class_data   = wide,
    class_cols         = class_cols,
    fixed_effects       = fixed_effects,
    random_effects      = random_effects,
    pathogen            = "bug",
    pathogen_col        = "pathogen",
    event_id_col        = "event_id",
    residual_structure  = residual_structure,
    prior_config        = utils::modifyList(list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0), prior_config),
    sampler_config      = list(chains = chains, iter_warmup = iter, iter_sampling = iter,
                               seed = seed, parallel_chains = chains, adapt_delta = 0.9),
    show_messages       = FALSE
  ))
}

test_that("Scenario 1: identity-generated + identity-fit -- PPC broadly calibrated", {
  skip_if_not(.ppcrec_cmdstan_available(), "cmdstanr/CmdStan not available")
  set.seed(101)
  n_hosp <- 4L; n_per_hosp <- 50L
  hosp_eff <- rnorm(n_hosp, sd = 0.3)
  wide <- do.call(rbind, lapply(seq_len(n_hosp), function(h) {
    n <- n_per_hosp
    age <- rnorm(n)
    mu1 <- 0.2 + 0.3 * age + hosp_eff[h]
    mu2 <- -0.1 + 0.2 * age + hosp_eff[h]
    data.frame(
      center_name = paste0("H", h), pathogen = "bug", Age_normalised = age,
      classA = rbinom(n, 1, stats::pnorm(mu1)), classB = rbinom(n, 1, stats::pnorm(mu2)),
      event_id = paste0("H", h, "_", seq_len(n))
    )
  }))

  fit <- .ppcrec_fit(wide, c("classA", "classB"), random_effects = c("center_name"),
                      residual_structure = "identity")
  stats_tbl <- compute_probit_ppc_statistics(
    fit, n_states = 150L, seed = 5L,
    statistics = c("marginal", "resistant_count", "pairwise", "complete_profile", "hospital_heterogeneity"),
    min_hospital_support = 20L, min_complete_profile_events = 20L
  )
  status <- compute_posterior_predictive_status(stats_tbl)
  expect_false(identical(status$status, "fail_major_ppc_misfit"))
})

test_that("Scenario 2: strong correlated residual generated + identity-fit -- pairwise/profile PPC should fail", {
  skip_if_not(.ppcrec_cmdstan_available(), "cmdstanr/CmdStan not available")
  set.seed(102)
  n <- 220L
  age <- rnorm(n)
  mu <- cbind(0.0 + 0.2 * age, 0.0 - 0.1 * age)
  rho <- 0.75
  Omega <- matrix(c(1, rho, rho, 1), 2, 2)
  L <- t(chol(Omega))
  Z <- mu + t(L %*% matrix(rnorm(2 * n), 2, n))
  Y <- (Z > 0) * 1L
  wide <- data.frame(
    center_name = "H1", pathogen = "bug", Age_normalised = age,
    classA = Y[, 1], classB = Y[, 2], event_id = paste0("E", seq_len(n))
  )

  fit <- .ppcrec_fit(wide, c("classA", "classB"), random_effects = c("center_name"),
                      residual_structure = "identity")
  stats_tbl <- compute_probit_ppc_statistics(
    fit, n_states = 150L, seed = 5L,
    statistics = c("marginal", "pairwise", "complete_profile"),
    min_complete_profile_events = 20L
  )
  status <- compute_posterior_predictive_status(stats_tbl)
  pairwise_bad <- !is.null(status$family_status$pairwise) && status$family_status$pairwise$status != "ok"
  profile_bad  <- !is.null(status$family_status$profile)  && status$family_status$profile$status != "ok"
  expect_true(pairwise_bad || profile_bad)
})

test_that("Scenario 3: strong hospital heterogeneity + heavily shrunk hospital RE -- hospital heterogeneity PPC should fail", {
  skip_if_not(.ppcrec_cmdstan_available(), "cmdstanr/CmdStan not available")
  set.seed(103)
  n_per_hosp <- 50L
  hosp_p <- c(0.08, 0.90, 0.10, 0.88)   # extreme, well-separated hospital rates
  wide <- do.call(rbind, lapply(seq_along(hosp_p), function(h) {
    n <- n_per_hosp
    data.frame(
      center_name = paste0("H", h), pathogen = "bug", Age_normalised = rnorm(n),
      classA = rbinom(n, 1, hosp_p[h]),
      classB = rbinom(n, 1, 0.3),   # unrelated nuisance class -- probit requires D >= 2
      event_id = paste0("H", h, "_", seq_len(n))
    )
  }))

  # tau_sd forced tiny: the declared hospital RE is present in the model but
  # its prior strongly shrinks estimated hospital effects toward zero,
  # emulating "the model does not effectively use hospital-level variation"
  # (Stage-1 architecture requires >=1 declared RE block, so a literal
  # omission of the hospital block is not expressible via this API).
  fit <- .ppcrec_fit(wide, c("classA", "classB"), random_effects = c("center_name"),
                      residual_structure = "identity", prior_config = list(tau_sd = 0.02))
  stats_tbl <- compute_probit_ppc_statistics(
    fit, n_states = 150L, seed = 5L,
    statistics = c("marginal", "hospital_heterogeneity"),
    min_hospital_support = 20L
  )
  status <- compute_posterior_predictive_status(stats_tbl)
  expect_true(!is.null(status$family_status$hospital_heterogeneity) &&
              status$family_status$hospital_heterogeneity$status != "ok")
})

test_that("Scenario 4: repeated-admission clustering, admission not declared as RE -- admission clustering PPC should fail", {
  skip_if_not(.ppcrec_cmdstan_available(), "cmdstanr/CmdStan not available")
  set.seed(104)
  n_adm <- 70L; ev_per_adm <- 3L
  adm_id <- rep(paste0("A", seq_len(n_adm)), each = ev_per_adm)
  adm_eff <- rep(rnorm(n_adm, sd = 1.8), each = ev_per_adm)
  n <- n_adm * ev_per_adm
  age <- rnorm(n)
  mu <- 0.1 * age + adm_eff
  wide <- data.frame(
    center_name = "H1", pathogen = "bug", Age_normalised = age,
    classA = rbinom(n, 1, stats::pnorm(mu)),
    classB = rbinom(n, 1, 0.3),   # unrelated nuisance class -- probit requires D >= 2
    event_id = paste0("E", seq_len(n)),
    admission_id = adm_id
  )

  fit <- .ppcrec_fit(wide, c("classA", "classB"), random_effects = c("center_name"),
                      residual_structure = "identity")
  stats_tbl <- compute_probit_ppc_statistics(
    fit, n_states = 150L, seed = 5L, statistics = "admission_clustering",
    random_effect_roles = c(admission_id = "admission")
  )
  status <- compute_posterior_predictive_status(stats_tbl)
  agree_row <- stats_tbl[stats_tbl$statistic_name == "cluster_within_group_same_class_agreement", ]
  expect_equal(agree_row$support_status, "supported")
  expect_true(!is.null(status$family_status$cluster) && status$family_status$cluster$status != "ok")
})

test_that("Scenario 5: correctly specified correlated + hospital/admission model -- core PPC broadly reproduced", {
  skip_if_not(.ppcrec_cmdstan_available(), "cmdstanr/CmdStan not available")
  set.seed(105)
  n_hosp <- 3L; n_adm_per_hosp <- 20L; ev_per_adm <- 2L
  rows <- list()
  for (h in seq_len(n_hosp)) {
    hosp_eff <- rnorm(1, sd = 0.3)
    for (a in seq_len(n_adm_per_hosp)) {
      adm_eff <- rnorm(1, sd = 0.4)
      age <- rnorm(ev_per_adm)
      mu1 <- 0.1 + 0.2 * age + hosp_eff + adm_eff
      mu2 <- -0.1 + 0.15 * age + hosp_eff + adm_eff
      rho <- 0.4
      Omega <- matrix(c(1, rho, rho, 1), 2, 2); L <- t(chol(Omega))
      Z <- cbind(mu1, mu2) + t(L %*% matrix(rnorm(2 * ev_per_adm), 2, ev_per_adm))
      Y <- (Z > 0) * 1L
      rows[[length(rows) + 1L]] <- data.frame(
        center_name = paste0("H", h), pathogen = "bug", Age_normalised = age,
        classA = Y[, 1], classB = Y[, 2],
        event_id = paste0("H", h, "_A", a, "_", seq_len(ev_per_adm)),
        admission_id = paste0("H", h, "_A", a)
      )
    }
  }
  wide <- do.call(rbind, rows)

  fit <- .ppcrec_fit(wide, c("classA", "classB"),
                      random_effects = list(
                        list(name = "hospital", group_col = "center_name", terms = "intercept"),
                        list(name = "admission", group_col = "admission_id", terms = "intercept")
                      ), residual_structure = "correlated")
  stats_tbl <- compute_probit_ppc_statistics(
    fit, n_states = 150L, seed = 5L,
    statistics = c("marginal", "resistant_count", "pairwise", "complete_profile",
                    "hospital_heterogeneity", "admission_clustering"),
    random_effect_roles = c(admission_id = "admission"),
    min_hospital_support = 15L, min_complete_profile_events = 15L
  )
  status <- compute_posterior_predictive_status(stats_tbl)
  expect_false(identical(status$status, "fail_major_ppc_misfit"))
})
