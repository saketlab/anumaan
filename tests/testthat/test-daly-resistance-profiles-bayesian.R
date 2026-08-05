# Tests for the Pathway 2 (Bayesian multivariate probit) resistance-profile
# pipeline: fit_bayesian_multivariate_probit(), compute_event_profile_probabilities(),
# and aggregate_profiles_for_daly(). These fit real (tiny) cmdstanr models, so
# they are skipped when cmdstanr / CmdStan are not available.

.bprofile_cmdstan_available <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE))
}

.bprofile_make_wide_data <- function(seed = 42, n_hosp = 3, n_ev_per_hosp = 80,
                                      class_cols = c("classA", "classB", "classC"),
                                      sparse_class_at_last_hosp = TRUE,
                                      all_died_at_last_hosp = FALSE) {
  set.seed(seed)
  make_hosp_data <- function(h, n) {
    base_p <- c(0.3, 0.6, 0.45)[seq_along(class_cols)] + (h - 2) * 0.05
    m <- sapply(base_p, function(p) rbinom(n, 1, p))
    colnames(m) <- class_cols
    df <- as.data.frame(m)
    df$center_name <- paste0("H", h)
    df$pathogen <- "bug"
    df$Age_normalised <- rnorm(n)
    df$gender <- sample(c("M", "F"), n, replace = TRUE)
    df$final_outcome <- if (all_died_at_last_hosp && h == n_hosp) {
      rep("Died", n)
    } else {
      sample(c("Discharged", "Died"), n, replace = TRUE, prob = c(0.85, 0.15))
    }
    df$event_id <- paste0("H", h, "_", seq_len(n))
    df
  }
  wide <- do.call(rbind, lapply(seq_len(n_hosp), make_hosp_data, n = n_ev_per_hosp))

  if (sparse_class_at_last_hosp && "classC" %in% class_cols) {
    idx <- which(wide$center_name == paste0("H", n_hosp))
    keep_tested <- sample(idx, 5)
    wide$classC[setdiff(idx, keep_tested)] <- NA
  }
  # ~30% missing on classA/classB so imputation is exercised
  for (cc in intersect(c("classA", "classB"), class_cols)) {
    na_idx <- sample(seq_len(nrow(wide)), floor(0.3 * nrow(wide)))
    wide[[cc]][na_idx] <- NA
  }
  wide
}

.bprofile_fit <- function(wide, class_cols, residual_structure = "identity",
                          chains = 2L, iter = 250L, seed = 1L) {
  suppressWarnings(fit_bayesian_multivariate_probit(
    event_class_data   = wide,
    class_cols         = class_cols,
    fixed_effects      = c("Age_normalised", "gender"),
    random_effects     = c("center_name"),
    pathogen           = "bug",
    pathogen_col       = "pathogen",
    event_id_col       = "event_id",
    outcome_col        = "final_outcome",
    residual_structure = residual_structure,
    prior_config       = list(beta_sd = 1.5, tau_sd = 1.0, lkj_eta = 2.0),
    sampler_config     = list(chains = chains, iter_warmup = iter, iter_sampling = iter,
                              seed = seed, parallel_chains = chains, adapt_delta = 0.9),
    show_messages      = FALSE
  ))
}

test_that("identity residual: fully observed event gets a degenerate profile", {
  skip_if_not(.bprofile_cmdstan_available(), "cmdstanr/CmdStan not available")

  wide <- .bprofile_make_wide_data()
  class_cols <- c("classA", "classB", "classC")
  fit <- .bprofile_fit(wide, class_cols, "identity")

  profs <- compute_event_profile_probabilities(
    fitted_model = fit, n_posterior_draws_for_profiles = 150,
    outcome_col = "final_outcome", seed = 1
  )

  ev_meta <- fit$event_metadata
  fully_obs_idx <- which(ev_meta$center_name != "H3" &
                          rowSums(is.na(ev_meta[, class_cols])) == 0)
  skip_if(length(fully_obs_idx) == 0L, "no fully-observed event in this draw of synthetic data")
  test_ev <- ev_meta$.event_idx[fully_obs_idx[1]]
  obs_row <- ev_meta[ev_meta$.event_idx == test_ev, class_cols]
  expected_label <- paste(ifelse(unlist(obs_row) == 1, "R", "S"), collapse = "")

  rows <- profs$event_profiles[profs$event_profiles$event_idx == test_ev, ]
  match_row <- rows[rows$profile_delta == expected_label, ]
  other_rows <- rows[rows$profile_delta != expected_label, ]

  expect_equal(match_row$profile_probability, 1, tolerance = 1e-9)
  expect_true(all(other_rows$profile_probability < 1e-9))
  expect_false(match_row$fully_model_imputed)
  expect_equal(match_row$n_classes_missing, 0L)
})

test_that("identity residual: partially observed event's profile distribution sums to 1 and respects observed cells", {
  skip_if_not(.bprofile_cmdstan_available(), "cmdstanr/CmdStan not available")

  wide <- .bprofile_make_wide_data()
  class_cols <- c("classA", "classB", "classC")
  fit <- .bprofile_fit(wide, class_cols, "identity")

  profs <- compute_event_profile_probabilities(
    fitted_model = fit, n_posterior_draws_for_profiles = 150,
    outcome_col = "final_outcome", seed = 1
  )

  ev_meta <- fit$event_metadata
  partial_idx <- which(ev_meta$center_name != "H3" &
                        rowSums(is.na(ev_meta[, class_cols])) > 0 &
                        rowSums(is.na(ev_meta[, class_cols])) < length(class_cols))
  skip_if(length(partial_idx) == 0L, "no partially-observed event in this draw of synthetic data")
  test_ev <- ev_meta$.event_idx[partial_idx[1]]
  obs_row <- ev_meta[ev_meta$.event_idx == test_ev, class_cols]

  rows <- profs$event_profiles[profs$event_profiles$event_idx == test_ev, ]
  expect_equal(sum(rows$profile_probability), 1, tolerance = 1e-6)
  expect_true(rows$n_classes_missing[1] > 0L)
  expect_false(rows$fully_model_imputed[1])

  known_cols <- names(obs_row)[!is.na(unlist(obs_row))]
  for (kc in known_cols) {
    observed_letter <- if (obs_row[[kc]] == 1) "R" else "S"
    pos <- match(kc, strsplit(rows$profile_class_set[1], "\\|")[[1]])
    mismatched <- rows[substr(rows$profile_delta, pos, pos) != observed_letter, ]
    expect_true(all(mismatched$profile_probability < 1e-9))
  }
})

test_that("identity residual: panel excludes marginally-ineligible class but not on pairwise grounds", {
  skip_if_not(.bprofile_cmdstan_available(), "cmdstanr/CmdStan not available")

  wide <- .bprofile_make_wide_data()
  class_cols <- c("classA", "classB", "classC")
  fit <- .bprofile_fit(wide, class_cols, "identity")

  elig <- fit$eligibility_report$marginal
  h3_classC <- elig[elig$center_name == "H3" & elig$antibiotic_class == "classC", ]
  expect_false(h3_classC$eligible)

  profs <- compute_event_profile_probabilities(
    fitted_model = fit, n_posterior_draws_for_profiles = 100,
    outcome_col = "final_outcome", seed = 1
  )

  h3_sets <- unique(profs$event_profiles$profile_class_set[profs$event_profiles$center_name == "H3"])
  expect_false(any(grepl("classC", h3_sets)))

  h3_rows <- profs$event_profiles[profs$event_profiles$center_name == "H3", ]
  expect_true(all(h3_rows$panel_eligibility_method == "marginal_only"))
  expect_true(all(grepl("classC", h3_rows$classes_excluded)))
})

test_that("aggregate_profiles_for_daly: zero-nonfatal cohort yields NA (not NaN) and an explicit exclusion reason", {
  skip_if_not(.bprofile_cmdstan_available(), "cmdstanr/CmdStan not available")

  wide <- .bprofile_make_wide_data(n_hosp = 2, class_cols = c("classA", "classB"),
                                    sparse_class_at_last_hosp = FALSE,
                                    all_died_at_last_hosp = TRUE)
  class_cols <- c("classA", "classB")
  fit <- .bprofile_fit(wide, class_cols, "identity")

  profs <- compute_event_profile_probabilities(
    fitted_model = fit, n_posterior_draws_for_profiles = 100,
    outcome_col = "final_outcome", seed = 1
  )
  daly_tbl <- aggregate_profiles_for_daly(
    profile_output = profs, hospital_col = fit$upper_re_col,
    pathogen_col = "pathogen", estimand = fit$estimand
  )

  h2_rows <- daly_tbl[daly_tbl$center_name == "H2", ]
  expect_true(all(h2_rows$n_events_nonfatal == 0L))
  expect_true(all(is.na(h2_rows$R_NF_mean)))
  expect_true(all(!is.nan(h2_rows$R_NF_mean[is.na(h2_rows$R_NF_mean)])))  # NA, not NaN
  expect_true(all(h2_rows$exclusion_reason_YLD == "no_nonfatal_events"))
  expect_true(all(!h2_rows$eligible_for_YLD))
  # H2's known-outcome cohort is fully populated (every event is "Died", a
  # known outcome) so YLL eligibility should be unaffected.
  expect_true(all(h2_rows$eligible_for_YLL))
})

test_that("correlated residual: profile output uses conditional Gibbs imputation and is DALY-eligible", {
  skip_if_not(.bprofile_cmdstan_available(), "cmdstanr/CmdStan not available")

  # Commit 6 replaced the old unconditional-simulation correlated path (which
  # was hard-flagged ineligible for DALY use) with conditional Gibbs
  # imputation on the truncated multivariate normal, conditioning on observed
  # AST signs. Correlated fits are therefore no longer categorically
  # DALY-ineligible -- eligibility now follows the same
  # sampler/events/draws-based rules as identity fits.
  wide <- .bprofile_make_wide_data(class_cols = c("classA", "classB", "classC"),
                                    sparse_class_at_last_hosp = FALSE)
  class_cols <- c("classA", "classB", "classC")
  fit <- .bprofile_fit(wide, class_cols, "correlated")

  profs <- compute_event_profile_probabilities(
    fitted_model = fit, n_posterior_draws_for_profiles = 150,
    outcome_col = "final_outcome", seed = 1
  )
  expect_true(all(profs$event_profiles$profile_generation_method ==
                  "conditional_gibbs_correlated"))

  # Profile probabilities must still sum to 1 per event under the Gibbs path,
  # exactly as under the identity/analytic path.
  by_event <- tapply(profs$event_profiles$profile_probability,
                      profs$event_profiles$event_idx, sum)
  expect_true(all(abs(by_event - 1) < 1e-6))

  daly_tbl <- aggregate_profiles_for_daly(
    profile_output = profs, hospital_col = fit$upper_re_col,
    pathogen_col = "pathogen", estimand = fit$estimand
  )
  expect_true(all(daly_tbl$eligible_for_profile_inference))
  expect_true(any(daly_tbl$eligible_for_YLL))
  expect_true(any(daly_tbl$eligible_for_YLD))
  expect_true(all(is.na(daly_tbl$exclusion_reason_YLL[daly_tbl$eligible_for_YLL])))
})
