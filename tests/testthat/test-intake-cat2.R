# Category 2: prep_intake.R

# prep_check_columns -----------------------------------------------------------

test_that("prep_check_columns: valid behavior", {
  dat <- cat2_intake_df()
  rep <- prep_check_columns(
    dat,
    required = c("patient_id", "date_of_culture", "final_outcome_date"),
    stop_on_missing = TRUE
  )
  expect_s3_class(rep, "data.frame")
  expect_true(all(c("column", "required", "present", "type_ok") %in% names(rep)))
})

test_that("prep_check_columns: edge cases", {
  dat <- cat2_intake_df()
  rep <- prep_check_columns(dat, expected_types = c(patient_id = "character"), stop_on_missing = TRUE)
  expect_true(any(rep$column == "patient_id"))
  expect_true(any(rep$column == "previous_hospitalisation"))
})

test_that("prep_check_columns: missingness", {
  dat <- cat2_intake_df()
  dat$final_outcome_date <- NULL
  expect_warning(
    prep_check_columns(dat, required = c("final_outcome_date"), stop_on_missing = FALSE),
    "Missing required"
  )
})

test_that("prep_check_columns: error handling", {
  expect_error(prep_check_columns(1:3), "must be a data frame")
})

# prep_check_keys --------------------------------------------------------------

test_that("prep_check_keys: valid behavior", {
  dat <- cat2_intake_df()
  rep <- prep_check_keys(dat, key_col = "patient_id")
  expect_s3_class(rep, "data.frame")
  expect_equal(rep$key_col, "patient_id")
  expect_equal(rep$n_missing_both_dates, 0)
})

test_that("prep_check_keys: edge cases", {
  dat <- data.frame(
    k = c("a", "a", "b", "c"),
    date_of_admission = c("2025-01-01", NA, "2025-01-03", NA),
    date_of_culture = c("2025-01-02", NA, "2025-01-04", NA),
    stringsAsFactors = FALSE
  )
  expect_warning(
    rep <- prep_check_keys(dat, key_col = "k"),
    "both 'date_of_admission' and 'date_of_culture' missing"
  )
  expect_equal(rep$n_duplicated, 1)
  expect_equal(rep$n_missing_both_dates, 2)
})

test_that("prep_check_keys: missingness", {
  dat <- data.frame(
    k = c("", "NA", NA, "x"),
    date_of_admission = c(NA, "", "2025-01-03", NA),
    date_of_culture = c(NA, "", NA, "2025-01-05"),
    stringsAsFactors = FALSE
  )
  expect_warning(
    rep <- prep_check_keys(dat, key_col = "k"),
    "both 'date_of_admission' and 'date_of_culture' missing"
  )
  expect_true(rep$n_missing >= 3)
  expect_true(rep$n_missing_both_dates >= 2)
})

test_that("prep_check_keys: error handling", {
  expect_warning(
    out <- prep_check_keys(data.frame(x = 1), key_col = "k"),
    "Key column 'k' not found"
  )
  expect_null(out)
})

# prep_validate_table ----------------------------------------------------------

test_that("prep_validate_table: valid behavior", {
  dat <- cat2_intake_df()
  out <- prep_validate_table(
    dat,
    required_cols = c("patient_id", "date_of_culture", "final_outcome_date"),
    key_col = "patient_id"
  )
  expect_type(out, "list")
  expect_true(all(c("data", "col_report", "key_report") %in% names(out)))
  expect_equal(out$key_report$n_missing_both_dates, 0)
})

test_that("prep_validate_table: edge cases", {
  dat <- cat2_intake_date_edge_cases()
  expect_warning(
    out <- prep_validate_table(
      dat,
      required_cols = c("patient_id", "date_of_admission", "date_of_culture", "final_outcome_date"),
      date_cols = c("date_of_admission", "date_of_culture", "final_outcome_date")
    ),
    "could not be parsed as dates"
  )
  expect_s3_class(out$data$date_of_admission, "Date")
  expect_s3_class(out$data$date_of_culture, "Date")
  expect_equal(out$data$date_of_admission[1], as.Date("2020-01-01"))
  expect_equal(out$data$date_of_admission[2], as.Date("2020-01-13"))
  expect_true(is.na(out$data$date_of_admission[3]))
  expect_true(is.na(out$data$date_of_admission[4]))
  expect_equal(out$data$date_of_culture[1], as.Date("2020-01-02"))
  expect_equal(out$data$date_of_culture[2], as.Date("2020-01-14"))
})

test_that("prep_validate_table: missingness", {
  dat <- data.frame(patient_id = c(NA, "p2"), date_of_culture = c(NA, "2020-01-02"), stringsAsFactors = FALSE)
  expect_warning(
    out <- prep_validate_table(dat, required_cols = c("patient_id", "date_of_culture"), key_col = "patient_id", stop_on_missing = FALSE),
    "Key 'patient_id'"
  )
  expect_s3_class(out$col_report, "data.frame")
})

test_that("prep_validate_table: error handling", {
  dat <- cat2_intake_df()
  expect_error(prep_validate_table(dat, required_cols = c("missing_col"), stop_on_missing = TRUE), "Missing required")
})

# validate_required_fields -----------------------------------------------------

test_that("validate_required_fields: valid behavior", {
  dat <- cat2_intake_df()
  out <- validate_required_fields(
    dat,
    required_cols = c("patient_id", "date_of_culture", "final_outcome_date"),
    stop_on_failure = FALSE
  )
  expect_true(out$valid)
})

test_that("validate_required_fields: edge cases", {
  dat <- data.frame(a = c(1, NA, 3), b = c(1, 2, 3))
  out <- validate_required_fields(dat, required_cols = c("a", "b"), min_completeness = 0.5, stop_on_failure = FALSE)
  expect_true(out$valid)
})

test_that("validate_required_fields: missingness", {
  dat <- data.frame(a = c(NA, NA, 1))
  out <- validate_required_fields(dat, required_cols = c("a"), min_completeness = 0.8, stop_on_failure = FALSE)
  expect_false(out$valid)
})

test_that("validate_required_fields: error handling", {
  dat <- data.frame(a = 1)
  expect_error(validate_required_fields(dat, required_cols = c("a", "b"), stop_on_failure = TRUE), "validation failed")
})

# validate_data_quality --------------------------------------------------------

test_that("validate_data_quality: valid behavior", {
  dat <- data.frame(patient_id = c("p1", "p2"), organism_normalized = c("e coli", "k pn"))
  out <- validate_data_quality(dat, min_rows = 2, stop_on_failure = FALSE)
  expect_true(out$passes_quality)
})

test_that("validate_data_quality: edge cases", {
  dat <- data.frame(patient_id = c("p1"), organism_normalized = c("e coli"))
  out <- validate_data_quality(dat, min_rows = 2, stop_on_failure = FALSE)
  expect_false(out$passes_quality)
})

test_that("validate_data_quality: missingness", {
  dat <- data.frame(patient_id = c("p1", NA), organism_normalized = c("e coli", NA))
  out <- validate_data_quality(dat, max_missing_pct = 0, stop_on_failure = FALSE)
  expect_false(out$passes_quality)
})

test_that("validate_data_quality: error handling", {
  dat <- data.frame(patient_id = c("p1"), organism_normalized = c("e coli"))
  expect_error(validate_data_quality(dat, min_rows = 2, stop_on_failure = TRUE), "failed")
})

# prep_log_source --------------------------------------------------------------

test_that("prep_log_source: valid behavior", {
  dat <- cat2_intake_df()
  out <- prep_log_source(dat, study_type = "generic", centre_name = "KGMU")
  prov <- attr(out, ".provenance")
  expect_type(prov, "list")
  expect_equal(prov$centre_name, "KGMU")
})

test_that("prep_log_source: edge cases", {
  dat <- cat2_intake_df()
  out <- prep_log_source(dat, study_type = "aiims_icu_bsi", centre_name = "AIIMS", sheet_name = "Sheet1")
  prov <- attr(out, ".provenance")
  expect_equal(prov$study_type, "aiims_icu_bsi")
  expect_equal(prov$sheet_name, "Sheet1")
})

test_that("prep_log_source: missingness", {
  dat <- cat2_intake_df()
  out <- prep_log_source(dat)
  prov <- attr(out, ".provenance")
  expect_true(is.na(prov$centre_name))
})

test_that("prep_log_source: error handling", {
  dat <- cat2_intake_df()
  out <- prep_log_source(dat, centre_name = "NONEXISTENT_CENTRE")
  prov <- attr(out, ".provenance")
  expect_equal(prov$centre_name, "NONEXISTENT_CENTRE")
})

# prep_inventory_columns -------------------------------------------------------

test_that("prep_inventory_columns: valid behavior", {
  dat <- cat2_intake_df()
  rep <- prep_inventory_columns(dat)
  expect_s3_class(rep, "tbl_df")
  expect_true(all(c("col_name", "class", "n_distinct", "pct_missing") %in% names(rep)))
})

test_that("prep_inventory_columns: edge cases", {
  dat <- data.frame(a = c(1, 1, 1), b = c("x", "", NA), stringsAsFactors = FALSE)
  rep <- prep_inventory_columns(dat)
  expect_true(any(rep$col_name == "a"))
})

test_that("prep_inventory_columns: missingness", {
  dat <- data.frame(a = c(NA, NA), b = c(NA, "x"), stringsAsFactors = FALSE)
  rep <- prep_inventory_columns(dat)
  expect_true(any(rep$pct_missing == 100))
})

test_that("prep_inventory_columns: error handling", {
  expect_error(prep_inventory_columns(NULL), NA)
})

# prep_detect_schema_drift -----------------------------------------------------

test_that("prep_detect_schema_drift: valid behavior", {
  d1 <- data.frame(a = 1, b = 2)
  d2 <- data.frame(a = 3, c = 4)
  out <- prep_detect_schema_drift(list(c1 = d1, c2 = d2))
  expect_s3_class(out, "tbl_df")
  expect_true(all(c("column", "present_in", "missing_from", "is_universal") %in% names(out)))
})

test_that("prep_detect_schema_drift: edge cases", {
  d1 <- data.frame(a = 1, b = 2)
  d2 <- data.frame(a = 3, b = 4)
  out <- prep_detect_schema_drift(list(c1 = d1, c2 = d2), reference_centre = "c1")
  expect_true(all(out$is_universal))
})

test_that("prep_detect_schema_drift: missingness", {
  d1 <- data.frame(a = 1)
  d2 <- data.frame()
  out <- prep_detect_schema_drift(list(c1 = d1, c2 = d2))
  expect_true(any(!out$is_universal))
})

test_that("prep_detect_schema_drift: error handling", {
  expect_error(prep_detect_schema_drift(list()), "non-empty")
})
