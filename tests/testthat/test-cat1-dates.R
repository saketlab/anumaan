# Category 1: prep_dates.R

# prep_parse_date_column -------------------------------------------------------

test_that("prep_parse_date_column: valid behavior", {
  x <- c("43831", "1577836800000", "2020-12-31")
  out <- prep_parse_date_column(x, col_name = "d", table_label = "t")
  expect_s3_class(out, "Date")
  expect_equal(out, as.Date(c("2020-01-01", "2020-01-01", "2020-12-31")))
})

test_that("prep_parse_date_column: edge cases", {
  out <- prep_parse_date_column(c("20201301", "2020/01/05", "01-31-2020"), col_name = "d", table_label = "t")
  # Numeric-like YYYYDDMM currently remains unparsed (NA) in this implementation.
  expect_true(is.na(out[1]))
  expect_equal(out[2], as.Date("2020-01-05"))
  expect_equal(out[3], as.Date("2020-01-31"))
})

test_that("prep_parse_date_column: missingness", {
  out <- prep_parse_date_column(c(NA, "", "NULL", "N/A"), col_name = "d", table_label = "t")
  expect_true(all(is.na(out)))
})

test_that("prep_parse_date_column: error handling", {
  expect_warning(
    out <- prep_parse_date_column(c("bad-date", "still_bad"), col_name = "d", table_label = "t"),
    "could not be parsed"
  )
  expect_true(all(is.na(out)))
})

# prep_coerce_dates ------------------------------------------------------------

test_that("prep_coerce_dates: valid behavior", {
  dat <- data.frame(date_of_culture = c("2020-01-01", "2020-01-02"), event_date = c("43831", "43832"), stringsAsFactors = FALSE)
  out <- prep_coerce_dates(dat)
  expect_s3_class(out$date_of_culture, "Date")
  expect_s3_class(out$event_date, "Date")
})

test_that("prep_coerce_dates: edge cases", {
  dat <- data.frame(collected_on = c("2020-01-01", "2020-01-02"), stringsAsFactors = FALSE)
  out <- prep_coerce_dates(dat, cols = "collected_on")
  expect_s3_class(out$collected_on, "Date")
})

test_that("prep_coerce_dates: missingness", {
  dat <- data.frame(date_of_culture = c(NA, "", "N/A"), stringsAsFactors = FALSE)
  out <- prep_coerce_dates(dat)
  expect_true(all(is.na(out$date_of_culture)))
})

test_that("prep_coerce_dates: error handling", {
  dat <- data.frame(name = c("a", "b"), stringsAsFactors = FALSE)
  out <- prep_coerce_dates(dat)
  expect_identical(out, dat)
})

# prep_validate_date_logic -----------------------------------------------------

test_that("prep_validate_date_logic: valid behavior", {
  dat <- cat1_min_dates()
  expect_no_warning(prep_validate_date_logic(dat))
})

test_that("prep_validate_date_logic: edge cases", {
  dat <- data.frame(
    admission_date = as.Date("2020-01-03"),
    culture_date = as.Date("2020-01-03"),
    outcome_date = as.Date("2020-01-03"),
    dob = as.Date("2020-01-03"),
    age = 0,
    final_outcome = "Survived",
    stringsAsFactors = FALSE
  )
  expect_no_warning(prep_validate_date_logic(dat))
})

test_that("prep_validate_date_logic: missingness", {
  dat <- data.frame(
    admission_date = as.Date(c("2020-01-01", NA)),
    culture_date = as.Date(c("2020-01-02", NA)),
    outcome_date = as.Date(c(NA, NA)),
    dob = as.Date(c(NA, NA)),
    age = c(NA, NA),
    final_outcome = c("Died", NA),
    stringsAsFactors = FALSE
  )
  expect_warning(prep_validate_date_logic(dat), "no outcome date")
})

test_that("prep_validate_date_logic: error handling", {
  bad <- data.frame(
    admission_date = as.Date("2020-01-03"),
    culture_date = as.Date("2020-01-02"),
    outcome_date = as.Date("2020-01-01"),
    dob = as.Date("2021-01-01"),
    age = 130,
    final_outcome = "Died",
    stringsAsFactors = FALSE
  )
  expect_warning(prep_validate_date_logic(bad), "admission_date > culture_date")
})
