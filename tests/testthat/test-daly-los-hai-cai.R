# Tests for LOS calculation in daly_compute_patient_los() and its callers
# daly_fit_los_rr() / daly_fit_los_rr_distribution() (R/daly_rr_and_los.R).
#
# Methodology: LOS = final outcome (discharge) date - admission date for every
# patient. LOS does NOT depend on the culture date or on HAI/CAI status.
# Only discharged patients are kept; LOS must be > 0 and <= max_los.
# HAI/CAI is still read directly from infection_type_col (not re-derived from
# dates) and used as a covariate in daly_fit_los_rr().

.synthetic_los_data <- function(infection_type_value) {
  set.seed(42)
  n <- 30
  admission <- as.Date("2026-01-01")
  culture <- admission + 7 # fixed 7-day gap, every patient
  # Jittered discharge offsets, independent of R/S assignment below, so each
  # group has non-degenerate within-group variance for the gamma fit.
  discharge <- culture + sample(15:25, n, replace = TRUE)
  data.frame(
    PatientInformation_id = as.character(seq_len(n)),
    center_name = "Hosp1",
    organism_name = "Test_Pathogen",
    syndrome = "BSI",
    infection_type = infection_type_value, # pre-computed, already final
    antibiotic_class = "A",
    antibiotic_name = "drugA",
    antibiotic_value = sample(c("R", "S"), n, replace = TRUE),
    date_of_admission = admission,
    date_of_first_positive_culture = culture,
    final_outcome_date = discharge,
    final_outcome = "Discharged",
    stringsAsFactors = FALSE
  )
}

# Minimal patient rows for direct daly_compute_patient_los() tests
.los_rows <- function(id, adm, out, outcome = "Discharged",
                      infection_type = "CAI", culture = NA) {
  data.frame(
    PatientInformation_id = id,
    center_name = "Hosp1",
    date_of_admission = as.Date(adm),
    final_outcome_date = as.Date(out),
    final_outcome = outcome,
    infection_type = infection_type,
    date_of_first_positive_culture = as.Date(culture),
    stringsAsFactors = FALSE
  )
}

.los <- function(d, ...) suppressMessages(daly_compute_patient_los(d, ...))

test_that("LOS = final outcome date - admission date", {
  out <- .los(.los_rows("P1", "2026-01-01", "2026-01-06"))
  expect_equal(out$LOS_days, 5)
  expect_named(out, c("PatientInformation_id", "center_name", "LOS_days"))
})

test_that("HAI, CAI and unknown infection types give the same LOS", {
  d <- rbind(
    .los_rows("HAI1", "2026-01-01", "2026-01-10", infection_type = "HAI", culture = "2026-01-04"),
    .los_rows("CAI1", "2026-01-01", "2026-01-10", infection_type = "CAI", culture = "2026-01-04"),
    .los_rows("UNK1", "2026-01-01", "2026-01-10", infection_type = "Unknown", culture = "2026-01-04")
  )
  out <- .los(d)
  expect_equal(out$LOS_days[match(c("HAI1", "CAI1", "UNK1"), out$PatientInformation_id)], c(9, 9, 9))
})

test_that("culture date does not change LOS", {
  a <- .los(.los_rows("P1", "2026-01-01", "2026-01-10", infection_type = "HAI", culture = "2026-01-02"))
  b <- .los(.los_rows("P1", "2026-01-01", "2026-01-10", infection_type = "HAI", culture = "2026-01-08"))
  c0 <- .los(.los_rows("P1", "2026-01-01", "2026-01-10", infection_type = "HAI", culture = NA))
  expect_equal(a$LOS_days, 9)
  expect_equal(b$LOS_days, 9)
  expect_equal(c0$LOS_days, 9)
})

test_that("missing admission or final outcome date gives NA LOS and is excluded", {
  d <- rbind(
    .los_rows("NOADM", NA, "2026-01-10"),
    .los_rows("NOOUT", "2026-01-01", NA),
    .los_rows("OK", "2026-01-01", "2026-01-10")
  )
  out <- .los(d)
  expect_equal(out$PatientInformation_id, "OK")
})

test_that("zero and negative LOS are excluded", {
  d <- rbind(
    .los_rows("ZERO", "2026-01-05", "2026-01-05"),
    .los_rows("NEG", "2026-01-10", "2026-01-05"),
    .los_rows("OK", "2026-01-01", "2026-01-03")
  )
  out <- .los(d)
  expect_equal(out$PatientInformation_id, "OK")
  expect_equal(out$LOS_days, 2)
})

test_that("LOS equal to max_los is kept and LOS above max_los is dropped", {
  d <- rbind(
    .los_rows("L200", "2026-01-01", as.character(as.Date("2026-01-01") + 200)),
    .los_rows("L201", "2026-01-01", as.character(as.Date("2026-01-01") + 201))
  )
  out <- .los(d, max_los = 200)
  expect_equal(out$PatientInformation_id, "L200")
  expect_equal(out$LOS_days, 200)
})

test_that("only discharged patients are kept and multi-row patients are not duplicated", {
  d <- rbind(
    .los_rows("DEAD", "2026-01-01", "2026-01-10", outcome = "Death"),
    .los_rows("D1", "2026-01-01", "2026-01-10", infection_type = "HAI", culture = "2026-01-03"),
    .los_rows("D1", "2026-01-01", "2026-01-10", infection_type = "CAI", culture = "2026-01-06")
  )
  out <- .los(d)
  expect_equal(out$PatientInformation_id, "D1")
  expect_equal(nrow(out), 1L)
  expect_equal(out$LOS_days, 9)
})

test_that("daly_fit_los_rr_distribution(): HAI and CAI with identical dates give identical LOS", {
  hai_data <- .synthetic_los_data("HAI")
  cai_data <- .synthetic_los_data("CAI")

  common_args <- list(
    patient_id_col = "PatientInformation_id",
    facility_col = "center_name",
    organism_col = "organism_name",
    syndrome_col = "syndrome",
    infection_type_col = "infection_type",
    antibiotic_class_col = "antibiotic_class",
    antibiotic_name_col = "antibiotic_name",
    antibiotic_value_col = "antibiotic_value",
    date_admission_col = "date_of_admission",
    date_discharge_col = "final_outcome_date",
    date_culture_col = "date_of_first_positive_culture",
    final_outcome_col = "final_outcome",
    min_n = 5
  )

  out_hai <- suppressMessages(do.call(daly_fit_los_rr_distribution, c(list(data = hai_data), common_args)))
  out_cai <- suppressMessages(do.call(daly_fit_los_rr_distribution, c(list(data = cai_data), common_args)))

  expect_equal(nrow(out_hai), 1L)
  expect_equal(nrow(out_cai), 1L)
  expect_equal(out_cai$median_LOS_R, out_hai$median_LOS_R)
  expect_equal(out_cai$median_LOS_S, out_hai$median_LOS_S)
  # every patient counted once, in exactly one of R / S
  expect_equal(out_hai$n_R + out_hai$n_S, nrow(hai_data))
  expect_equal(out_hai$n_R, sum(hai_data$antibiotic_value == "R"))
})

test_that("daly_fit_los_rr() works with the new LOS output and keeps HAI as a covariate", {
  set.seed(7)
  n <- 120
  adm <- as.Date("2026-01-01") + sample(0:30, n, replace = TRUE)
  res <- sample(c("R", "S"), n, replace = TRUE)
  inf <- sample(c("HAI", "CAI"), n, replace = TRUE)
  los <- round(rgamma(n, shape = 3, rate = 0.3) * ifelse(res == "R", 1.5, 1)) + 1
  d <- data.frame(
    PatientInformation_id = sprintf("P%03d", seq_len(n)),
    center_name = "H1",
    organism_name = "Klebsiella pneumoniae",
    syndrome = "BSI",
    infection_type = inf,
    antibiotic_class = "Carbapenems",
    antibiotic_name = "Meropenem",
    antibiotic_value = res,
    date_of_admission = adm,
    final_outcome_date = adm + los,
    final_outcome = "Discharged",
    Age = sample(18:85, n, replace = TRUE),
    Gender = sample(c("M", "F"), n, replace = TRUE)
  )
  # no culture-date column: no longer required
  out <- suppressWarnings(suppressMessages(daly_fit_los_rr(d, min_n = 10)))
  expect_equal(nrow(out), 1L)
  expect_equal(out$n_patients, n)
  expect_equal(out$n_resistant, sum(res == "R"))
  expect_equal(out$n_susceptible, sum(res == "S"))
  expect_equal(out$mean_los_resistant, mean(los[res == "R"]))
  expect_equal(out$mean_los_susceptible, mean(los[res == "S"]))
})

test_that("daly_fit_los_rr_distribution() defaults infection_type_col to 'infection_type'", {
  expect_equal(
    formals(daly_fit_los_rr_distribution)$infection_type_col,
    "infection_type"
  )
})

test_that("daly_fit_los_rr() defaults infection_type_col to 'infection_type'", {
  expect_equal(
    formals(daly_fit_los_rr)$infection_type_col,
    "infection_type"
  )
})

test_that("daly_fit_los_rr_distribution() errors clearly when infection_type_col is missing", {
  bad_data <- .synthetic_los_data("HAI")
  bad_data$infection_type <- NULL
  expect_error(
    daly_fit_los_rr_distribution(
      bad_data,
      patient_id_col = "PatientInformation_id",
      facility_col = "center_name",
      organism_col = "organism_name",
      syndrome_col = "syndrome",
      infection_type_col = "infection_type",
      antibiotic_class_col = "antibiotic_class",
      antibiotic_name_col = "antibiotic_name",
      antibiotic_value_col = "antibiotic_value",
      date_admission_col = "date_of_admission",
      date_discharge_col = "final_outcome_date",
      date_culture_col = "date_of_first_positive_culture",
      final_outcome_col = "final_outcome"
    ),
    "infection_type"
  )
})
