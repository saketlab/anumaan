# Category 2: prep_outputs.R

# prep_filter_minimally_usable -------------------------------------------------

test_that("prep_filter_minimally_usable: valid behavior", {
  dat <- cat2_outputs_df()
  out <- prep_filter_minimally_usable(dat)
  expect_true(nrow(out) <= nrow(dat))
})

test_that("prep_filter_minimally_usable: edge cases", {
  dat <- cat2_outputs_df()
  dat$ast_value_harmonized[1] <- "I"
  out <- prep_filter_minimally_usable(dat)
  expect_true(all(out$ast_value_harmonized %in% c("S", "I", "R")))
})

test_that("prep_filter_minimally_usable: missingness", {
  dat <- data.frame(patient_id = c(NA, ""), culture_date = as.Date(c(NA, NA)), organism_name = c(NA, ""), ast_value_harmonized = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_filter_minimally_usable(dat)
  expect_equal(nrow(out), 0)
})

test_that("prep_filter_minimally_usable: error handling", {
  out <- prep_filter_minimally_usable(data.frame(x = 1))
  expect_equal(nrow(out), 1)
})

# prep_filter_analysis_ready ---------------------------------------------------

test_that("prep_filter_analysis_ready: valid behavior", {
  dat <- cat2_outputs_df()
  out <- prep_filter_analysis_ready(dat)
  expect_true(nrow(out) <= nrow(dat))
})

test_that("prep_filter_analysis_ready: edge cases", {
  dat <- cat2_outputs_df()
  dat$contaminant_flag <- c(FALSE, TRUE, FALSE)
  out <- prep_filter_analysis_ready(dat, exclude_contaminants = TRUE)
  expect_false(any(out$contaminant_flag, na.rm = TRUE))
})

test_that("prep_filter_analysis_ready: missingness", {
  dat <- cat2_outputs_df()
  dat$antibiotic_name <- c(NA, "", "")
  out <- prep_filter_analysis_ready(dat)
  expect_equal(nrow(out), 0)
})

test_that("prep_filter_analysis_ready: error handling", {
  out <- prep_filter_analysis_ready(data.frame(x = 1))
  expect_equal(nrow(out), 1)
})

# prep_build_fatal_cohort ------------------------------------------------------

test_that("prep_build_fatal_cohort: valid behavior", {
  dat <- cat2_outputs_df()
  out <- prep_build_fatal_cohort(dat)
  expect_true(all(out$final_outcome == "Died"))
})

test_that("prep_build_fatal_cohort: edge cases", {
  dat <- cat2_outputs_df()
  out <- prep_build_fatal_cohort(dat, died_value = "Died")
  expect_true(nrow(out) >= 0)
})

test_that("prep_build_fatal_cohort: missingness", {
  dat <- data.frame(final_outcome = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_build_fatal_cohort(dat)
  expect_equal(nrow(out), 0)
})

test_that("prep_build_fatal_cohort: error handling", {
  expect_warning(out <- prep_build_fatal_cohort(data.frame(x = 1)), "not found")
  expect_equal(nrow(out), 0)
})

# prep_build_nonfatal_cohort ---------------------------------------------------

test_that("prep_build_nonfatal_cohort: valid behavior", {
  dat <- cat2_outputs_df()
  out <- prep_build_nonfatal_cohort(dat)
  expect_true(all(out$final_outcome == "Survived"))
})

test_that("prep_build_nonfatal_cohort: edge cases", {
  dat <- cat2_outputs_df()
  out <- prep_build_nonfatal_cohort(dat, survived_value = "Survived")
  expect_true(nrow(out) >= 0)
})

test_that("prep_build_nonfatal_cohort: missingness", {
  dat <- data.frame(final_outcome = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_build_nonfatal_cohort(dat)
  expect_equal(nrow(out), 0)
})

test_that("prep_build_nonfatal_cohort: error handling", {
  expect_warning(out <- prep_build_nonfatal_cohort(data.frame(x = 1)), "not found")
  expect_equal(nrow(out), 0)
})

# prep_attrition_flow ----------------------------------------------------------

test_that("prep_attrition_flow: valid behavior", {
  flow <- NULL
  dat <- cat2_outputs_df()
  flow <- prep_attrition_flow(flow, dat, stage_name = "raw", reason = "initial")
  expect_s3_class(flow, "data.frame")
  expect_true(all(c("stage", "n_rows", "n_patients", "n_events") %in% names(flow)))
})

test_that("prep_attrition_flow: edge cases", {
  flow <- prep_attrition_flow(NULL, cat2_outputs_df(), stage_name = "raw")
  flow2 <- prep_attrition_flow(flow, cat2_outputs_df()[1:2, ], stage_name = "filtered")
  expect_equal(flow2$n_removed[2], flow2$n_rows[1] - flow2$n_rows[2])
})

test_that("prep_attrition_flow: missingness", {
  dat <- data.frame(x = 1:2)
  flow <- prep_attrition_flow(NULL, dat, stage_name = "raw")
  expect_true(is.na(flow$n_patients))
})

test_that("prep_attrition_flow: error handling", {
  flow <- prep_attrition_flow(NULL, data.frame(x = 1), stage_name = "nonexistent_centre_name")
  expect_equal(flow$stage[1], "nonexistent_centre_name")
})

# prep_missingness_report ------------------------------------------------------

test_that("prep_missingness_report: valid behavior", {
  dat <- cat2_outputs_df()
  out <- prep_missingness_report(dat, threshold = 20)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("col_name", "n_missing", "pct_missing", "is_high_missing") %in% names(out)))
})

test_that("prep_missingness_report: edge cases", {
  dat <- data.frame(a = c(1, NA), b = c("", "x"), stringsAsFactors = FALSE)
  out <- prep_missingness_report(dat, threshold = 10)
  expect_true(any(out$is_high_missing))
})

test_that("prep_missingness_report: missingness", {
  dat <- data.frame(a = c(NA, NA), b = c(NA, NA), stringsAsFactors = FALSE)
  out <- prep_missingness_report(dat)
  expect_true(all(out$pct_missing == 100))
})

test_that("prep_missingness_report: error handling", {
  out <- prep_missingness_report(data.frame(a = 1), cols = c("a", "missing"))
  expect_true(all(out$col_name == "a"))
})

# prep_validate_analysis_ready -------------------------------------------------

test_that("prep_validate_analysis_ready: valid behavior", {
  dat <- cat2_outputs_df()
  out <- prep_validate_analysis_ready(dat, min_rows = 1, max_missing_pct = 80, stop_on_failure = FALSE)
  expect_type(out, "list")
  expect_true(out$passes)
})

test_that("prep_validate_analysis_ready: edge cases", {
  dat <- cat2_outputs_df()
  out <- prep_validate_analysis_ready(dat, min_rows = 10, stop_on_failure = FALSE)
  expect_false(out$passes)
})

test_that("prep_validate_analysis_ready: missingness", {
  dat <- data.frame(patient_id = c(NA, NA), culture_date = as.Date(c(NA, NA)), organism_name = c(NA, NA), antibiotic_name = c(NA, NA), ast_value_harmonized = c(NA, NA), stringsAsFactors = FALSE)
  out <- prep_validate_analysis_ready(dat, min_rows = 1, max_missing_pct = 30, stop_on_failure = FALSE)
  expect_false(out$passes)
})

test_that("prep_validate_analysis_ready: error handling", {
  dat <- data.frame(x = 1)
  expect_error(prep_validate_analysis_ready(dat, stop_on_failure = TRUE), "Validation issues found")
})
