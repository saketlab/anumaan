# Tests for the distribution-based median LOS workflow in
# daly_fit_los_rr_distribution() (R/daly_rr_and_los.R):
#   - .compute_median_from_fit(): theoretical medians
#       Weibull    lambda * log(2)^(1/k)
#       Log-normal exp(meanlog)
#       Gamma      qgamma(0.5, shape, rate)
#   - AIC-based selection of the fitted distribution
#   - RR_LOS = median_LOS_R / median_LOS_S
#   - min_n, failed fits, missing/invalid LOS
#   - propagation of median_LOS_R/S into daly_assign_rr_to_profiles()

.fake_fit <- function(...) list(estimate = c(...))

test_that(".compute_median_from_fit() returns the theoretical Weibull median", {
  med <- .compute_median_from_fit(.fake_fit(shape = 1.7, scale = 12), "weibull")
  expect_equal(med, 12 * log(2)^(1 / 1.7))
  expect_equal(med, stats::qweibull(0.5, shape = 1.7, scale = 12))
})

test_that(".compute_median_from_fit() returns the theoretical log-normal median", {
  med <- .compute_median_from_fit(.fake_fit(meanlog = 2.4, sdlog = 0.6), "lnorm")
  expect_equal(med, exp(2.4))
  expect_equal(med, stats::qlnorm(0.5, meanlog = 2.4, sdlog = 0.6))
})

test_that(".compute_median_from_fit() returns the gamma median via qgamma", {
  med <- .compute_median_from_fit(.fake_fit(shape = 3, rate = 0.25), "gamma")
  expect_equal(med, stats::qgamma(0.5, shape = 3, rate = 0.25))
  # the gamma median is below its mean (shape / rate) for a right-skewed fit
  expect_lt(med, 3 / 0.25)
})

test_that(".compute_median_from_fit() returns NA for invalid parameters, unknown distributions and errors", {
  expect_true(is.na(.compute_median_from_fit(.fake_fit(shape = -1, scale = 12), "weibull")))
  expect_true(is.na(.compute_median_from_fit(.fake_fit(shape = 2, scale = 0), "weibull")))
  expect_true(is.na(.compute_median_from_fit(.fake_fit(meanlog = 2, sdlog = 0), "lnorm")))
  expect_true(is.na(.compute_median_from_fit(.fake_fit(meanlog = NA, sdlog = 1), "lnorm")))
  expect_true(is.na(.compute_median_from_fit(.fake_fit(shape = 2, rate = -0.1), "gamma")))
  expect_true(is.na(.compute_median_from_fit(.fake_fit(shape = 2, rate = 0.5), "exp")))
  expect_true(is.na(.compute_median_from_fit(NULL, "gamma")))
})

# Long-format fixture: one class "A", resistance drives LOS.
.median_fixture <- function(n = 200, seed = 3, los_r = NULL, los_s = NULL) {
  set.seed(seed)
  res <- rep(c("R", "S"), length.out = n)
  los <- ifelse(res == "R",
    round(stats::rlnorm(n, log(20), 0.5)),
    round(stats::rlnorm(n, log(10), 0.5))
  )
  los[los < 1] <- 1
  if (!is.null(los_r)) los[res == "R"] <- los_r
  if (!is.null(los_s)) los[res == "S"] <- los_s
  adm <- as.Date("2025-01-01")
  data.frame(
    PatientInformation_id = sprintf("P%03d", seq_len(n)),
    center_name = "H1",
    organism_name = "Klebsiella pneumoniae",
    infection_type = "CAI",
    antibiotic_class = "A",
    antibiotic_name = "drugA",
    antibiotic_value = res,
    date_of_admission = adm,
    final_outcome_date = adm + los,
    final_outcome = "Discharged",
    stringsAsFactors = FALSE
  )
}

.fit_dist <- function(d, ...) {
  suppressWarnings(suppressMessages(daly_fit_los_rr_distribution(d, ...)))
}

test_that("daly_fit_los_rr_distribution() returns median columns and RR = ratio of medians", {
  d <- .median_fixture()
  out <- .fit_dist(d, min_n = 10)
  expect_equal(nrow(out), 1L)
  expect_true(all(c("median_LOS_R", "median_LOS_S", "median_LOS_R_years", "median_LOS_S_years") %in% names(out)))
  expect_false(any(grepl("^mean_LOS", names(out))))

  # medians match the gamma fit of each group (default distribution)
  los <- as.numeric(d$final_outcome_date - d$date_of_admission)
  fit_r <- safe_fit(los[d$antibiotic_value == "R"], "gamma")
  fit_s <- safe_fit(los[d$antibiotic_value == "S"], "gamma")
  med_r <- stats::qgamma(0.5, fit_r$estimate[["shape"]], fit_r$estimate[["rate"]])
  med_s <- stats::qgamma(0.5, fit_s$estimate[["shape"]], fit_s$estimate[["rate"]])
  expect_equal(out$median_LOS_R, round(med_r, 3L))
  expect_equal(out$median_LOS_S, round(med_s, 3L))
  expect_equal(out$RR_LOS, round(med_r / med_s, 4L))
  expect_equal(out$median_LOS_R_years, round(med_r / 365, 6L))
  expect_equal(out$n_R + out$n_S, nrow(d))
  expect_gt(out$RR_LOS, 1)
})

test_that("daly_fit_los_rr_distribution() selects the lowest-AIC distribution and uses its median", {
  d <- .median_fixture(n = 400, seed = 9)
  dists <- c("weibull", "lnorm", "gamma")
  out <- .fit_dist(d, min_n = 10, distributions = dists)
  expect_equal(nrow(out), 1L)

  los <- as.numeric(d$final_outcome_date - d$date_of_admission)
  los_r <- los[d$antibiotic_value == "R"]
  fits <- Filter(Negate(is.null), setNames(lapply(dists, safe_fit, x = los_r), dists))
  best <- names(which.min(sapply(fits, `[[`, "aic")))
  expect_equal(out$best_dist_R, best)
  expect_equal(out$median_LOS_R, round(.compute_median_from_fit(fits[[best]], best), 3L))
})

test_that("daly_fit_los_rr_distribution() skips a class when a group is below min_n", {
  d <- .median_fixture(n = 30)
  out <- .fit_dist(d, min_n = 20) # 15 R and 15 S
  expect_equal(nrow(out), 0L)
  expect_true(all(c("median_LOS_R", "median_LOS_S") %in% names(out)))
})

test_that("daly_fit_los_rr_distribution() skips a class when the distribution fit fails", {
  # all resistant patients have identical LOS -> safe_fit() returns NULL
  d <- .median_fixture(n = 60, los_r = 12)
  out <- .fit_dist(d, min_n = 10)
  expect_equal(nrow(out), 0L)
})

test_that("daly_fit_los_rr_distribution() ignores missing and invalid LOS", {
  d <- .median_fixture(n = 100)
  d$final_outcome_date[1:3] <- NA # missing discharge date
  d$date_of_admission[4:5] <- NA # missing admission date
  d$final_outcome_date[6:7] <- d$date_of_admission[6:7] - 2 # negative LOS
  d$final_outcome_date[8] <- d$date_of_admission[8] # zero LOS
  out <- .fit_dist(d, min_n = 10)
  expect_equal(nrow(out), 1L)
  expect_equal(out$n_R + out$n_S, nrow(d) - 8L)
  expect_true(is.finite(out$median_LOS_R) && is.finite(out$median_LOS_S))
})

test_that("median_LOS_R/S propagate into daly_assign_rr_to_profiles()", {
  d <- .median_fixture()
  rr <- .fit_dist(d, min_n = 10)
  rr$pathogen <- "Test_Pathogen"
  profiles_output <- list(Test_Pathogen = list(
    profiles = data.frame(
      profile = c("S", "A"), A = c(0L, 1L), probability = c(0.6, 0.4),
      stringsAsFactors = FALSE
    ),
    classes = "A"
  ))
  out <- daly_assign_rr_to_profiles(
    profiles_output, rr,
    los_r_col = "median_LOS_R", los_s_col = "median_LOS_S", los_unit = "days"
  )$Test_Pathogen
  expect_equal(out$LOSR_years[out$profile == "A"], round(rr$median_LOS_R / 365, 6L))
  expect_equal(out$LOSS_years[out$profile == "A"], round(rr$median_LOS_S / 365, 6L))
  expect_true(is.na(out$LOSR_years[out$profile == "S"]))
  expect_equal(out$RR_LOS_profile[out$profile == "A"], rr$RR_LOS)
})
