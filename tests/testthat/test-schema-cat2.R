# Category 2: prep_schema.R

# prep_standardize_column_names ------------------------------------------------

test_that("prep_standardize_column_names: valid behavior", {
  dat <- cat2_schema_df()
  mapping <- list(patient_id = c("PID"), culture_date = c("CultureDate"), organism_name = c("Organism"))
  out <- prep_standardize_column_names(dat, mapping = mapping, fuzzy_match = FALSE)
  expect_type(out, "list")
  expect_true(all(c("patient_id", "culture_date", "organism_name") %in% names(out$data)))
})

test_that("prep_standardize_column_names: edge cases", {
  dat <- cat2_schema_df()
  mapping <- list(patient_id = c("PID"), culture_date = c("culture_date_typo"))
  out <- prep_standardize_column_names(dat, mapping = mapping, fuzzy_match = TRUE, fuzzy_threshold = 0.4)
  expect_true("patient_id" %in% names(out$data))
})

test_that("prep_standardize_column_names: missingness", {
  dat <- data.frame(PID = c("p1", NA), stringsAsFactors = FALSE)
  out <- prep_standardize_column_names(dat, mapping = list(patient_id = c("PID")), fuzzy_match = FALSE)
  expect_true("patient_id" %in% names(out$data))
})

test_that("prep_standardize_column_names: error handling", {
  dat <- cat2_schema_df()
  out <- prep_standardize_column_names(dat, mapping = list(x = c("does_not_exist")), fuzzy_match = FALSE)
  expect_true("x" %in% names(out$mapping_log) || length(out$mapping_log) == 0)
})

# prep_build_column_map --------------------------------------------------------

test_that("prep_build_column_map: valid behavior", {
  dat <- cat2_schema_df()
  m <- prep_build_column_map(dat, column_map = c(patient_id = "PID", culture_date = "CultureDate"))
  expect_equal(unname(m), c("PID", "CultureDate"))
})

test_that("prep_build_column_map: edge cases", {
  dat <- cat2_schema_df()
  m <- prep_build_column_map(dat, column_map = c(patient_id = "PID"), custom_map = c(organism_name = "Organism"))
  expect_true(all(c("patient_id", "organism_name") %in% names(m)))
})

test_that("prep_build_column_map: missingness", {
  dat <- cat2_schema_df()
  m <- prep_build_column_map(dat, column_map = c(patient_id = "PID", bad = "missing"))
  expect_false("bad" %in% names(m))
})

test_that("prep_build_column_map: error handling", {
  dat <- cat2_schema_df()
  expect_error(prep_build_column_map(dat, column_map = list(a = "PID")), "named character vector")
})

# prep_apply_column_map --------------------------------------------------------

test_that("prep_apply_column_map: valid behavior", {
  dat <- cat2_schema_df()
  out <- prep_apply_column_map(dat, c(patient_id = "PID", culture_date = "CultureDate"))
  expect_true(all(c("patient_id", "culture_date") %in% names(out)))
})

test_that("prep_apply_column_map: edge cases", {
  dat <- cat2_schema_df()
  out <- prep_apply_column_map(dat, c(patient_id = "PID", sample_type = "Sample_Type"))
  expect_true("sample_type" %in% names(out))
})

test_that("prep_apply_column_map: missingness", {
  dat <- cat2_schema_df()
  out <- prep_apply_column_map(dat, c(patient_id = "MISSING"))
  expect_false("patient_id" %in% names(out))
})

test_that("prep_apply_column_map: error handling", {
  dat <- cat2_schema_df()
  out <- prep_apply_column_map(dat, character())
  expect_identical(out, dat)
})

# prep_assert_standard_names ---------------------------------------------------

test_that("prep_assert_standard_names: valid behavior", {
  dat <- data.frame(patient_id = "p1", culture_date = as.Date("2020-01-01"))
  expect_invisible(prep_assert_standard_names(dat, c("patient_id", "culture_date"), strict = TRUE))
})

test_that("prep_assert_standard_names: edge cases", {
  dat <- data.frame(patient_id = "p1")
  expect_warning(prep_assert_standard_names(dat, c("patient_id", "culture_date"), strict = FALSE), "missing")
})

test_that("prep_assert_standard_names: missingness", {
  dat <- data.frame(patient_id = c("p1", NA))
  expect_invisible(prep_assert_standard_names(dat, c("patient_id"), strict = TRUE))
})

test_that("prep_assert_standard_names: error handling", {
  dat <- data.frame(patient_id = "p1")
  expect_error(prep_assert_standard_names(dat, c("patient_id", "culture_date"), strict = TRUE), "missing")
})

# detect_preprocessing_capabilities --------------------------------------------

test_that("detect_preprocessing_capabilities: valid behavior", {
  dat <- data.frame(patient_id = "p1", culture_date = as.Date("2020-01-01"), organism_name = "e coli", sample_type = "blood", antibiotic_name = "amikacin", antibiotic_value = "R", stringsAsFactors = FALSE)
  caps <- detect_preprocessing_capabilities(dat)
  expect_type(caps, "logical")
  expect_true("build_outputs" %in% names(caps))
})

test_that("detect_preprocessing_capabilities: edge cases", {
  dat <- data.frame(patient_id = "p1", dob = as.Date("2000-01-01"), culture_date = as.Date("2020-01-01"))
  caps <- detect_preprocessing_capabilities(dat)
  expect_true(caps[["derive_age"]])
})

test_that("detect_preprocessing_capabilities: missingness", {
  caps <- detect_preprocessing_capabilities(data.frame(x = 1))
  expect_false(any(caps))
})

test_that("detect_preprocessing_capabilities: error handling", {
  expect_error(detect_preprocessing_capabilities(NULL), NA)
})

# prep_report_capabilities -----------------------------------------------------

test_that("prep_report_capabilities: valid behavior", {
  caps <- c(parse_dates = TRUE, derive_hai = FALSE)
  out <- prep_report_capabilities(caps)
  expect_identical(out, caps)
})

test_that("prep_report_capabilities: edge cases", {
  dat <- data.frame(patient_id = "p1", culture_date = as.Date("2020-01-01"), organism_name = "e coli")
  out <- prep_report_capabilities(dat)
  expect_type(out, "logical")
})

test_that("prep_report_capabilities: missingness", {
  out <- prep_report_capabilities(data.frame())
  expect_type(out, "logical")
})

test_that("prep_report_capabilities: error handling", {
  expect_error(prep_report_capabilities(1), NA)
})

# prep_clean_optional_columns --------------------------------------------------

test_that("prep_clean_optional_columns: valid behavior", {
  dat <- data.frame(infection_type = c("community acquired", "hospital acquired"), unit_type = c("general ward", "icu"), stringsAsFactors = FALSE)
  out <- prep_clean_optional_columns(dat)
  expect_equal(out$infection_type, c("CAI", "HAI"))
  expect_equal(out$unit_type, c("Ward", "ICU"))
})

test_that("prep_clean_optional_columns: edge cases", {
  dat <- data.frame(comorbidities = c("DM", "HTN"), previous_history = c("x", "y"), stringsAsFactors = FALSE)
  out <- prep_clean_optional_columns(dat)
  expect_equal(out$comorbidities, c("DM", "HTN"))
})

test_that("prep_clean_optional_columns: missingness", {
  dat <- data.frame(infection_type = c("", NA), stringsAsFactors = FALSE)
  out <- prep_clean_optional_columns(dat)
  expect_true(all(is.na(out$infection_type)))
})

test_that("prep_clean_optional_columns: error handling", {
  out <- prep_clean_optional_columns(data.frame(x = 1))
  expect_identical(out, data.frame(x = 1))
})
