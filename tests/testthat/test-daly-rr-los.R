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
  out <- suppressWarnings(suppressMessages(daly_fit_los_rr(
    rr_cohort(), min_n = 10, hospital_effect = "stratified", ci_method = "delta"
  )))$rr
  expect_setequal(out$hospital, c("H1", "H2"))
  expect_true(all(out$RR_LOS > 1.2 & out$RR_LOS < 2.1))
  expect_true(all(out$CI_lower < out$RR_LOS & out$RR_LOS < out$CI_upper))
  expect_equal(out$RR_LOS, exp(out$beta_log_rr))
  expect_true(all(out$mean_los_resistant > out$mean_los_susceptible))
})

test_that("daly_fit_los_rr hospital fixed effect shares RR and gives per-hospital LOSR/LOSS", {
  fit <- suppressWarnings(suppressMessages(daly_fit_los_rr(
    rr_cohort(n = 400), min_n = 10, ci_method = "bootstrap", n_boot = 50,
    seed = 3, return_draws = TRUE
  )))
  rr <- fit$rr
  h <- fit$hospital
  expect_equal(nrow(rr), 1L)
  expect_true(is.na(rr$hospital))
  expect_true(rr$RR_LOS > 1.2 & rr$RR_LOS < 2.1)
  expect_true(rr$CI_lower < rr$RR_LOS & rr$RR_LOS < rr$CI_upper)
  expect_setequal(h$hospital, c("H1", "H2"))
  # LOSR / LOSS reproduces the shared ratio in every hospital
  expect_equal(h$LOSR / h$LOSS, rep(rr$RR_LOS, 2), tolerance = 1e-8)
  expect_equal(rr$LOSR / rr$LOSS, rr$RR_LOS, tolerance = 1e-8)
  expect_true(all(h$LOSR_lower < h$LOSR & h$LOSR < h$LOSR_upper))
  expect_true(all(h$excess_lower < h$excess_LOS & h$excess_LOS < h$excess_upper))
  expect_true(is.finite(rr$RxHospital_p))
  d <- fit$draws[["Klebsiella pneumoniae | Carbapenems"]]
  expect_length(d$RR_LOS, rr$n_boot_ok)
  expect_equal(dim(d$LOSR), c(rr$n_boot_ok, 2L))
  expect_equal(colnames(d$LOSS), c("H1", "H2"))
})

test_that("daly_fit_los_rr: delta and bootstrap give similar RR intervals", {
  args <- list(rr_cohort(n = 400), min_n = 10)
  a <- suppressWarnings(suppressMessages(do.call(daly_fit_los_rr, c(args, ci_method = "delta"))))$rr
  b <- suppressWarnings(suppressMessages(do.call(daly_fit_los_rr, c(args, ci_method = "bootstrap", n_boot = 200, seed = 9))))$rr
  expect_equal(a$RR_LOS, b$RR_LOS)
  expect_equal(a$CI_lower, b$CI_lower, tolerance = 0.05)
  expect_equal(a$CI_upper, b$CI_upper, tolerance = 0.05)
  n <- suppressWarnings(suppressMessages(do.call(daly_fit_los_rr, c(args, ci_method = "none"))))$rr
  expect_true(is.na(n$CI_lower))
})

test_that("daly_fit_los_rr uses only the pathogen's own results and its earliest culture", {
  base <- rr_cohort(n = 200, seed = 5)
  # every admission also has an E. coli culture that is carbapenem-resistant,
  # taken later and labelled with the opposite infection type
  other <- base
  other$organism_name <- "Escherichia coli"
  other$antibiotic_value <- "R"
  other$date_of_first_positive_culture <- other$date_of_first_positive_culture + 3
  other$infection_type <- ifelse(base$infection_type == "HAI", "CAI", "HAI")
  both <- rbind(base, other)
  one <- suppressWarnings(suppressMessages(daly_fit_los_rr(
    base, min_n = 10, ci_method = "none", organism_name = "Klebsiella pneumoniae"
  )))$rr
  two <- suppressWarnings(suppressMessages(daly_fit_los_rr(
    both, min_n = 10, ci_method = "none", organism_name = "Klebsiella pneumoniae"
  )))$rr
  expect_equal(two$n_resistant, one$n_resistant)
  expect_equal(two$RR_LOS, one$RR_LOS)
})

test_that("daly_fit_los_rr classes accepts a per-pathogen panel and reports skipped classes", {
  d <- rr_cohort()
  fit <- suppressWarnings(suppressMessages(daly_fit_los_rr(
    d, min_n = 10, ci_method = "none",
    classes = list("klebsiella pneumoniae" = c("Carbapenems", "Colistin"))
  )))
  expect_equal(fit$rr$antibiotic_class, "Carbapenems")
  expect_equal(fit$skipped$antibiotic_class, "Colistin")
})

test_that("daly_fit_los_rr deprecated arguments warn", {
  expect_warning(
    suppressMessages(daly_fit_los_rr(rr_cohort(), patient_id_col = "PatientInformation_id",
                                     min_n = 10, ci_method = "none")),
    "deprecated"
  )
})

test_that(".los_admission_covariates takes HAI from the earliest culture of the pathogen", {
  d <- data.frame(
    id = c("A", "A", "B"), org = "K",
    cul = as.Date(c("2024-01-09", "2024-01-02", "2024-01-05")),
    inf = c("HAI", "CAI", "Unknown"), age = c(40, 40, 60), sex = c("M", "M", "F")
  )
  cv <- .los_admission_covariates(d, "id", "org", "cul", "inf", "HAI", "CAI",
                                  "age", "sex", NULL, NULL, NULL, NULL)
  expect_equal(cv$HAI[cv$.id == "A"], 0L) # earliest culture (2 Jan) is CAI
  expect_true(is.na(cv$HAI[cv$.id == "B"])) # "Unknown" -> missing
})

test_that(".los_class_status applies the worst phenotype and the colistin intermediate rule", {
  d <- data.frame(
    id = c("A", "A", "B", "C", "D"), org = "K",
    cls = c("Carbapenems", "Carbapenems", "Carbapenems", "Colistin", "Colistin"),
    drug = c("meropenem", "imipenem", "meropenem", "colistin", "colistin"),
    val = c("S", "R", "I", "I", "S")
  )
  st <- .los_class_status(d, "id", "org", "cls", "drug", "val", intermediate_as_resistant = TRUE)
  st <- st[order(st$.id), ]
  expect_equal(st$R, c(1L, 1L, 0L, 0L)) # A any R; B I -> R; C colistin I -> S; D S
  st2 <- .los_class_status(d, "id", "org", "cls", "drug", "val", intermediate_as_resistant = FALSE)
  expect_false(any(st2$.id %in% c("B", "C")))
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

test_that("daly_fit_mortality_rr builds resistance from each pathogen's own results", {
  base <- rr_cohort(n = 300, seed = 21)
  other <- base
  other$organism_name <- "Escherichia coli"
  other$antibiotic_value <- "R"
  one <- suppressWarnings(suppressMessages(daly_fit_mortality_rr(
    base, min_n = 10, min_deaths = 5, organism_name = "Klebsiella pneumoniae", facility_col = NULL
  )))
  both <- suppressWarnings(suppressMessages(daly_fit_mortality_rr(
    rbind(base, other), min_n = 10, min_deaths = 5, facility_col = NULL
  )))
  kp <- both[both$pathogen == "Klebsiella pneumoniae", ]
  expect_equal(kp$n_resistant, one$n_resistant)
  expect_equal(kp$OR_death, one$OR_death)
})
