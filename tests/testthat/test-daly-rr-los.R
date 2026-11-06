# Relative risk (mortality, LOS) and LOS distribution helpers

# Synthetic cohort where resistance multiplies LOS by 1.6 and raises death
# probability from 0.12 to 0.35 (true OR ~= 3.9).
rr_cohort <- function(n = 300, seed = 11) {
  set.seed(seed)
  res <- sample(c("R", "S"), n, replace = TRUE, prob = c(0.4, 0.6))
  adm <- as.Date("2024-01-01") + sample(0:200, n, replace = TRUE)
  cult <- adm + sample(0:5, n, replace = TRUE)
  los <- round(rgamma(n, shape = 3, rate = 0.3) * ifelse(res == "R", 1.6, 1)) + 1
  died <- stats::rbinom(n, 1, ifelse(res == "R", 0.35, 0.12)) == 1
  data.frame(
    PatientInformation_id = sprintf("P%03d", seq_len(n)),
    center_name = sample(c("H1", "H2"), n, replace = TRUE),
    organism_name = "Klebsiella pneumoniae",
    syndrome = "BSI",
    infection_type = ifelse(cult - adm > 2, "HAI", "CAI"),
    type_of_infection = ifelse(cult - adm > 2, "HAI", "CAI"),
    antibiotic_class = "Carbapenems",
    antibiotic_name = "Meropenem",
    antibiotic_value = res,
    unit_type = sample(c("ICU", "Ward"), n, replace = TRUE),
    date_of_admission = adm,
    date_of_first_positive_culture = cult,
    final_outcome_date = cult + los,
    final_outcome = ifelse(died, "Death", "Discharged"),
    Age = sample(18:85, n, replace = TRUE),
    Gender = sample(c("M", "F"), n, replace = TRUE)
  )
}

test_that("daly_fit_los_rr recovers the simulated LOS ratio per hospital", {
  out <- suppressWarnings(suppressMessages(daly_fit_los_rr(rr_cohort(), min_n = 10)))
  expect_setequal(out$hospital, c("H1", "H2"))
  expect_true(all(out$RR_LOS > 1.2 & out$RR_LOS < 2.1))
  expect_true(all(out$CI_lower < out$RR_LOS & out$RR_LOS < out$CI_upper))
  expect_equal(out$RR_LOS, exp(out$beta_log_rr))
  expect_true(all(out$mean_los_resistant > out$mean_los_susceptible))
})

test_that("daly_fit_mortality_rr recovers a raised odds of death for resistance", {
  out <- suppressWarnings(suppressMessages(
    daly_fit_mortality_rr(rr_cohort(), min_n = 10, min_deaths = 5)
  ))
  expect_setequal(out$hospital, c("H1", "H2"))
  expect_true(all(out$OR_death > 1.5))
  expect_equal(out$OR_death, exp(out$beta_logit))
  expect_true(all(out$OR_CI_lower > 1))
  # counterfactual RR is bounded by the OR when baseline risk > 0
  expect_true(all(out$RR_death > 1 & out$RR_death < out$OR_death))
})

test_that("daly_derive_hai_cai_for_mortality applies the 48 h admission-to-culture rule", {
  d <- rr_cohort()[1:4, ]
  d$type_of_infection <- NA
  d$date_of_first_positive_culture <- d$date_of_admission + c(0, 1, 3, 5)
  out <- suppressMessages(daly_derive_hai_cai_for_mortality(d))
  expect_equal(out$infection_type_derived, c("CAI", "CAI", "HAI", "HAI"))
})

test_that("daly_compute_patient_los and get_los_by_resistance split discharged LOS by phenotype", {
  los <- data.frame(
    PatientInformation_id = c("A", "A", "B", "C"),
    center_name = "H1",
    date_of_admission = as.Date(c("2024-01-01", "2024-01-01", "2024-01-05", "2024-01-10")),
    final_outcome_date = as.Date(c("2024-01-11", "2024-01-11", "2024-01-08", "2024-12-31")),
    final_outcome = c("Discharged", "Discharged", "Death", "Discharged")
  )
  prepped <- suppressMessages(daly_compute_patient_los(los))
  # deduplicated, deaths dropped, LOS above max_los (200) dropped
  expect_equal(prepped$PatientInformation_id, "A")
  expect_equal(prepped$LOS_days, 10)

  d <- rr_cohort()
  abx <- transform(d, organism_clean = organism_name, sus = antibiotic_value)
  los_d <- suppressMessages(daly_compute_patient_los(d))
  split <- suppressMessages(get_los_by_resistance(
    abx, los_d,
    organism = "Klebsiella pneumoniae"
  ))
  expect_named(split, c("R", "S"))
  expect_gt(mean(split$R), mean(split$S))

  # centre filter works because daly_compute_patient_los() returns center_name
  split_h1 <- suppressMessages(get_los_by_resistance(
    abx, los_d,
    organism = "Klebsiella pneumoniae", center = "H1"
  ))
  expect_named(split_h1, c("R", "S"))
  expect_lte(length(split_h1$R) + length(split_h1$S), length(split$R) + length(split$S))
})

test_that("LOS distribution fitting, summary, AIC comparison and plot work together", {
  set.seed(1)
  x <- rgamma(200, shape = 2, rate = 0.3)
  fits <- suppressMessages(fit_distributions(x))
  expect_named(fits, c("weibull", "lnorm", "gamma"))

  s <- suppressMessages(summarise_distribution(fits$gamma, "gamma"))
  expect_equal(s$Mean_LOS, 2 / 0.3, tolerance = 0.15)

  aic <- suppressMessages(compare_distribution_aic(x, fits))
  expect_equal(unname(unlist(aic)), unname(vapply(fits, `[[`, numeric(1), "aic")))

  expect_s3_class(plot_los_distributions(x, "LOS", fits = fits), "ggplot")
})

test_that("daly_load_rr_reference reads the bundled GBD RR table", {
  rr <- suppressMessages(daly_load_rr_reference())
  expect_named(rr, c("pathogen", "drug", "sample_size", "mean_rr", "lower_bound", "upper_bound"))
  expect_gt(nrow(rr), 50L)
  expect_true(all(rr$lower_bound <= rr$mean_rr & rr$mean_rr <= rr$upper_bound, na.rm = TRUE))
})

test_that("daly_add_rr_mappings maps organisms and classes onto RR categories", {
  d <- data.frame(
    organism_normalized = c("Klebsiella pneumoniae", "Escherichia coli"),
    antibiotic_class = c("Carbapenems", "Third-generation-cephalosporins")
  )
  out <- suppressMessages(daly_add_rr_mappings(d))
  expect_equal(out$rr_pathogen, tolower(d$organism_normalized)) # lower-case internal key
  expect_equal(out$rr_drug, d$antibiotic_class)
})

test_that("daly_filter_profiles_to_rr_classes keeps reference + matchable profiles and renormalises", {
  prof <- list(K = data.frame(
    profile = c("SS", "RS", "SR", "RR"), A = c(0, 1, 0, 1), B = c(0, 0, 1, 1),
    probability = c(0.4, 0.2, 0.3, 0.1)
  ))
  rr <- data.frame(pathogen = "K", antibiotic_class = "A")
  out <- suppressMessages(daly_filter_profiles_to_rr_classes(prof, rr))
  kp <- out$K
  expect_s3_class(kp, "data.frame")
  expect_setequal(kp$profile, c("SS", "RS", "RR"))
  expect_equal(sum(kp$probability), 1)
  expect_equal(kp$probability[kp$profile == "SS"], 0.4 / 0.7)

  expect_warning(suppressMessages(daly_filter_profiles_to_rr_classes(
    prof, data.frame(pathogen = "K", antibiotic_class = "Z")
  )), "no profile classes")
})

test_that(".encode_comorbidity_mortality handles numeric, binary, ordinal and absent inputs", {
  enc <- function(x) {
    suppressMessages(.encode_comorbidity_mortality(
      data.frame(PatientInformation_id = seq_along(x), comorbidities = x)
    ))
  }
  num <- enc(c(0, 2, 5))
  expect_equal(attr(num, "comorbidity_encoding"), "numeric")
  expect_equal(num$comorbidity_encoded, c(0, 2, 5))

  bin <- enc(c("present", "none", "unknown"))
  expect_equal(attr(bin, "comorbidity_encoding"), "binary")
  expect_equal(as.numeric(bin$comorbidity_encoded), c(1, 0, NA))

  ord <- enc(c("none", "mild", "severe"))
  expect_equal(attr(ord, "comorbidity_encoding"), "ordinal")

  absent <- suppressMessages(.encode_comorbidity_mortality(data.frame(PatientInformation_id = 1)))
  expect_equal(attr(absent, "comorbidity_encoding"), "absent")
})

test_that("daly_fit_mortality_rr accepts a comorbidity covariate", {
  d <- rr_cohort()
  d$comorbidities <- sample(0:6, nrow(d), replace = TRUE)
  out <- suppressWarnings(suppressMessages(daly_fit_mortality_rr(d,
    comorbidity_col = "comorbidities",
    min_n = 10, min_deaths = 5
  )))
  expect_true(all(out$OR_death > 1))
})
