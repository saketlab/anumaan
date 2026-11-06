# Tests for daly_fit_los_distn_by_profile() (Section 5.3, profile-specific median LOS)
# in R/daly_rr_and_los.R. Profiles are SUPPLIED (as from the profile_* pipeline);
# the function must not derive them.

# n patients per profile; LOS drawn from a log-normal with the given medians.
.profile_fixture <- function(medians = c(SS = 10, RS = 15, RR = 25),
                             n = c(SS = 60, RS = 40, RR = 30), seed = 5) {
  set.seed(seed)
  prof <- rep(names(medians), times = n[names(medians)])
  los <- round(stats::rlnorm(length(prof), log(medians[prof]), 0.4))
  los[los < 1] <- 1
  ids <- sprintf("P%03d", seq_along(prof))
  adm <- as.Date("2025-01-01")
  data <- data.frame(
    PatientInformation_id = ids,
    center_name = "H1",
    organism_name = "Klebsiella pneumoniae",
    date_of_admission = adm,
    final_outcome_date = adm + los,
    final_outcome = "Discharged",
    stringsAsFactors = FALSE
  )
  profiles <- data.frame(
    PatientInformation_id = ids,
    profile_delta = prof,
    stringsAsFactors = FALSE
  )
  list(data = data, profiles = profiles, los = stats::setNames(los, ids), prof = prof)
}

.fit_prof <- function(...) suppressWarnings(suppressMessages(daly_fit_los_distn_by_profile(...)))

test_that("one row per profile with the best-AIC median and RR vs the all-S reference", {
  fx <- .profile_fixture()
  out <- .fit_prof(fx$data, fx$profiles)

  expect_setequal(out$profile_delta, c("SS", "RS", "RR"))
  expect_true(all(out$status == "fitted"))
  expect_equal(out$profile_delta[out$is_reference], "SS")
  expect_equal(out$n_patients[out$profile_delta == "RR"], 30L)

  # median for each profile = theoretical median of its lowest-AIC fit
  dists <- c("weibull", "lnorm", "gamma")
  for (pr in c("SS", "RS", "RR")) {
    x <- fx$los[fx$prof == pr]
    fits <- Filter(Negate(is.null), stats::setNames(lapply(dists, safe_fit, x = x), dists))
    best <- names(which.min(sapply(fits, `[[`, "aic")))
    row <- out[out$profile_delta == pr, ]
    expect_equal(row$best_dist, best)
    expect_equal(row$median_LOS, round(.compute_median_from_fit(fits[[best]], best), 3L))
  }

  ref <- out$median_LOS[out$profile_delta == "SS"]
  expect_true(all(out$median_LOS_ref == ref))
  expect_equal(out$RR_LOS, round(out$median_LOS / ref, 4L))
  expect_equal(out$RR_LOS[out$profile_delta == "SS"], 1)
  expect_gt(out$RR_LOS[out$profile_delta == "RR"], out$RR_LOS[out$profile_delta == "RS"])
})

test_that("profiles below min_n are reported as insufficient_n with NA estimates", {
  fx <- .profile_fixture(n = c(SS = 60, RS = 40, RR = 5))
  out <- .fit_prof(fx$data, fx$profiles, min_n = 10)
  rr <- out[out$profile_delta == "RR", ]
  expect_equal(rr$status, "insufficient_n")
  expect_true(is.na(rr$median_LOS) && is.na(rr$RR_LOS))
  expect_equal(rr$n_los, 5L)
})

test_that("a profile whose distribution fit fails is reported as fit_failed", {
  fx <- .profile_fixture()
  rs_ids <- fx$profiles$PatientInformation_id[fx$profiles$profile_delta == "RS"]
  fx$data$final_outcome_date[fx$data$PatientInformation_id %in% rs_ids] <-
    fx$data$date_of_admission[fx$data$PatientInformation_id %in% rs_ids] + 12 # identical LOS
  out <- .fit_prof(fx$data, fx$profiles)
  expect_equal(out$status[out$profile_delta == "RS"], "fit_failed")
  expect_true(is.na(out$RR_LOS[out$profile_delta == "RS"]))
})

test_that("profiles are taken as supplied, patients without a profile are excluded", {
  fx <- .profile_fixture()
  drop <- fx$profiles$PatientInformation_id[1:5] # 5 SS patients lose their profile
  out <- .fit_prof(fx$data, fx$profiles[!fx$profiles$PatientInformation_id %in% drop, ])
  expect_equal(out$n_patients[out$profile_delta == "SS"], 55L)
  expect_message(
    daly_fit_los_distn_by_profile(fx$data, fx$profiles[!fx$profiles$PatientInformation_id %in% drop, ]),
    "5 patient\\(s\\) with valid LOS have no profile"
  )
})

test_that("patients with conflicting profiles are excluded", {
  fx <- .profile_fixture()
  extra <- data.frame(PatientInformation_id = fx$profiles$PatientInformation_id[1:3],
                      profile_delta = "RR", stringsAsFactors = FALSE)
  out <- .fit_prof(fx$data, rbind(fx$profiles, extra))
  expect_equal(out$n_patients[out$profile_delta == "SS"], 57L)
})

test_that("custom column names, explicit reference and pathogen-specific profiles work", {
  fx <- .profile_fixture()
  prof <- data.frame(
    event_id = fx$profiles$PatientInformation_id,
    pathogen = "Klebsiella pneumoniae",
    state = fx$profiles$profile_delta,
    stringsAsFactors = FALSE
  )
  out <- .fit_prof(fx$data, prof,
    profile_id_col = "event_id", profile_col = "state",
    profile_pathogen_col = "pathogen", ref_profile = "RS"
  )
  expect_equal(out$profile_delta[out$is_reference], "RS")
  expect_equal(out$RR_LOS[out$profile_delta == "RS"], 1)
  expect_equal(unique(out$pathogen), "Klebsiella pneumoniae")
})

test_that("no all-S profile gives NA RR_LOS but still fits medians", {
  fx <- .profile_fixture(medians = c(RS = 15, RR = 25), n = c(RS = 40, RR = 30))
  out <- .fit_prof(fx$data, fx$profiles)
  expect_true(all(is.na(out$RR_LOS)))
  expect_true(all(!out$is_reference))
  expect_true(all(is.finite(out$median_LOS)))
})

test_that("only discharged patients with valid LOS are used", {
  fx <- .profile_fixture()
  ss <- which(fx$profiles$profile_delta == "SS")
  fx$data$final_outcome[ss[1:4]] <- "Death"
  fx$data$final_outcome_date[ss[5:6]] <- NA
  fx$data$final_outcome_date[ss[7]] <- fx$data$date_of_admission[ss[7]] - 1
  out <- .fit_prof(fx$data, fx$profiles)
  expect_equal(out$n_los[out$profile_delta == "SS"], 60L - 7L)
})

test_that("missing profile columns give a clear error", {
  fx <- .profile_fixture()
  expect_error(
    daly_fit_los_distn_by_profile(fx$data, fx$profiles[, "PatientInformation_id", drop = FALSE]),
    "profile_delta"
  )
})
