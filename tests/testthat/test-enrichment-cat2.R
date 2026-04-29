# Category 2: prep_derivation.R enrichment functions

# prep_fill_age ----------------------------------------------------------------

test_that("prep_fill_age: valid behavior", {
  dat <- cat2_enrich_df()
  out <- prep_fill_age(dat)
  expect_true("Age" %in% names(out))
  expect_false(is.na(out$Age[1]))
})

test_that("prep_fill_age: edge cases", {
  dat <- cat2_enrich_df()
  dat$Age <- c(20, 35)
  out <- prep_fill_age(dat, overwrite = FALSE)
  expect_equal(out$Age, c(20, 35))
})

test_that("prep_fill_age: missingness", {
  dat <- data.frame(Age = c(NA, NA), stringsAsFactors = FALSE)
  out <- prep_fill_age(dat, age_col = "Age", dob_col = "DOB", date_col = "date_of_culture")
  expect_true(all(is.na(out$Age)))
})

test_that("prep_fill_age: error handling", {
  dat <- data.frame(x = 1)
  out <- prep_fill_age(dat, age_col = "Age")
  expect_true("Age" %in% names(out))
})

# prep_infer_department --------------------------------------------------------

test_that("prep_infer_department: valid behavior", {
  dat <- cat2_enrich_df()
  out <- prep_infer_department(dat)
  expect_true("hospital_department" %in% names(out))
  expect_true("department_method" %in% names(out))
})

test_that("prep_infer_department: edge cases", {
  dat <- data.frame(Age = c(5, 40), specimen_type = c("blood", "peritoneal abscess"), diagnosis_1 = c("", ""), stringsAsFactors = FALSE)
  out <- prep_infer_department(dat)
  expect_equal(out$hospital_department[1], "Pediatrics")
})

test_that("prep_infer_department: missingness", {
  dat <- data.frame(hospital_department = c(NA, NA), stringsAsFactors = FALSE)
  out <- prep_infer_department(dat)
  expect_true(all(is.na(out$hospital_department)))
})

test_that("prep_infer_department: error handling", {
  dat <- data.frame(x = 1)
  out <- prep_infer_department(dat)
  expect_true("hospital_department" %in% names(out))
})

# prep_derive_los_from_dates ---------------------------------------------------

test_that("prep_derive_los_from_dates: valid behavior", {
  dat <- data.frame(admission_date = as.Date(c("2020-01-01", "2020-01-05")), outcome_date = as.Date(c("2020-01-03", "2020-01-05")), stringsAsFactors = FALSE)
  out <- prep_derive_los_from_dates(dat)
  expect_equal(out$los_days, c(2, 0))
})

test_that("prep_derive_los_from_dates: edge cases", {
  dat <- cat2_enrichment_los_date_edge_cases()
  expect_warning(
    dat <- prep_coerce_dates(dat, cols = c("admission_date", "outcome_date"), table_label = "cat2_los_edge_cases"),
    "could not be parsed as dates"
  )
  expect_warning(out <- prep_derive_los_from_dates(dat), "negative LOS")
  expect_equal(out$los_days[1], 2)
  expect_equal(out$los_days[2], 4)
  expect_true(is.na(out$los_days[3]))
  expect_equal(out$los_days[4], -1)
  expect_equal(out$los_days[5], 0)
})

test_that("prep_derive_los_from_dates: missingness", {
  dat <- data.frame(unit_admission_date = as.Date(c("2020-01-01", NA)), unit_duration_days = c(3, NA), stringsAsFactors = FALSE)
  out <- prep_derive_los_from_dates(dat)
  expect_equal(out$los_days[1], 3)
  expect_true(is.na(out$los_days[2]))
})

test_that("prep_derive_los_from_dates: error handling", {
  dat <- data.frame(x = 1)
  out <- prep_derive_los_from_dates(dat)
  expect_true("los_days" %in% names(out))
})

# prep_derive_dob_from_components ----------------------------------------------

test_that("prep_derive_dob_from_components: valid behavior", {
  dat <- data.frame(dob_year = 2000, dob_month = 1, dob_day = 15, stringsAsFactors = FALSE)
  out <- prep_derive_dob_from_components(dat)
  expect_equal(out$dob, as.Date("2000-01-15"))
})

test_that("prep_derive_dob_from_components: edge cases", {
  dat <- data.frame(age_years = 20, admission_date = as.Date("2020-01-01"), stringsAsFactors = FALSE)
  out <- prep_derive_dob_from_components(dat)
  expect_s3_class(out$dob, "Date")
})

test_that("prep_derive_dob_from_components: missingness", {
  dat <- data.frame(dob_year = NA, dob_month = NA, dob_day = NA, age_years = NA, admission_date = as.Date(NA), stringsAsFactors = FALSE)
  out <- prep_derive_dob_from_components(dat)
  expect_true(all(is.na(out$dob)))
})

test_that("prep_derive_dob_from_components: error handling", {
  dat <- data.frame(x = 1)
  out <- prep_derive_dob_from_components(dat)
  expect_true("dob" %in% names(out))
})
