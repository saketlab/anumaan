# Category 1: prep_derivation.R

# get_age_bins -----------------------------------------------------------------

test_that("get_age_bins: valid behavior", {
  gbd <- get_age_bins("GBD_standard")
  expect_true(is.character(gbd))
  expect_true("<1" %in% gbd)
})

test_that("get_age_bins: edge cases", {
  ped <- get_age_bins("pediatric")
  ger <- get_age_bins("geriatric")
  neo <- get_age_bins("neonatal")
  expect_true("18+" %in% ped)
  expect_true("90+" %in% ger)
  expect_true("<0.02" %in% neo)
  expect_true("0.02-0.08" %in% neo)
  expect_true("18+" %in% neo)
})

test_that("get_age_bins: missingness", {
  expect_error(get_age_bins(NA_character_), "Unknown age bin type")
})

test_that("get_age_bins: error handling", {
  expect_error(get_age_bins("wrong_mode_name"), "Unknown age bin type")
})

# prep_assign_age_bins ---------------------------------------------------------

test_that("prep_assign_age_bins: valid behavior", {
  dat <- cat1_age_years_edge()[c(4, 9, 12), c("age"), drop = FALSE]
  out <- prep_assign_age_bins(dat, age_col = "age")
  expect_true("Age_bin" %in% names(out))
  expect_equal(as.character(out$Age_bin), c("<1", "10-15", "85+"))
})

test_that("prep_assign_age_bins: edge cases - age_unit conversion", {
  # Age in months: 6 months -> 0.5 years -> "<1" bin
  dat <- data.frame(Age = c(6, 18, 36), stringsAsFactors = FALSE)
  out <- prep_assign_age_bins(dat, age_unit = "months")
  expect_equal(as.character(out$Age_bin), c("<1", "1-5", "1-5"))

  # Age in days: 10 days -> 0.027 years -> "0.02-0.08" bin (neonatal preset)
  dat2 <- data.frame(Age = c(3, 15, 45, 200, 400), stringsAsFactors = FALSE)
  out2 <- prep_assign_age_bins(dat2, age_unit = "days", bins = "neonatal")
  expect_equal(as.character(out2$Age_bin), c("<0.02", "0.02-0.08", "0.08-0.25", "0.25-1", "1-5"))
})

test_that("prep_assign_age_bins: edge cases - compound columns", {
  # Reuse fixture with year+month+day components.
  dat <- cat1_age_from_components()[1:3, c("age_years", "age_months", "age_days"), drop = FALSE]
  out <- prep_assign_age_bins(
    dat,
    age_col        = "age_years",
    bins           = "neonatal",
    age_months_col = "age_months",
    age_days_col   = "age_days"
  )
  expect_equal(as.character(out$Age_bin[1]), "<0.02")      # 3 days
  expect_equal(as.character(out$Age_bin[2]), "0.02-0.08")  # 15 days
  expect_equal(as.character(out$Age_bin[3]), "0.08-0.25")  # 1 month
})

test_that("prep_assign_age_bins: missingness", {
  # NA in primary age stays NA; NA in component columns treated as 0.
  dat <- cat1_age_from_components()[c(13, 24), c("age_years", "age_months", "age_days"), drop = FALSE]
  out <- prep_assign_age_bins(
    dat,
    age_col = "age_years",
    age_months_col = "age_months",
    age_days_col = "age_days",
    bins = "neonatal"
  )
  expect_false(is.na(out$Age_bin[1]))
  expect_true(is.na(out$Age_bin[2]))
})

test_that("prep_assign_age_bins: age-years edge values from shared fixture", {
  dat <- cat1_age_years_edge()[c(1, 2, 4, 16, 21), c("age"), drop = FALSE]
  dat$year <- c(0, 0, dat$age[3], dat$age[4], dat$age[5])
  dat$months <- c(0, 3, 0, 0, 0)
  dat$age_days <- c(10, 0, 0, 0, 0)

  out <- prep_assign_age_bins(
    dat,
    age_col = "age",
    bins = "neonatal",
    negative_age_strategy = "fallback",
    fallback_years_col = "year",
    fallback_months_col = "months",
    fallback_days_col = "age_days"
  )
  expect_equal(as.character(out$Age_bin[1]), "0.02-0.08")  # recovered from 0y 0m 10d
  expect_equal(as.character(out$Age_bin[2]), "0.25-1")     # recovered from 0y 3m 0d
  expect_equal(as.character(out$Age_bin[3]), "0.08-0.25") # 0.2 years
  expect_equal(as.character(out$Age_bin[4]), "0.02-0.08") # 0.04 years
  expect_equal(as.character(out$Age_bin[5]), "0.02-0.08") # 0.02 years
})

test_that("prep_assign_age_bins: negative ages can be forced to NA", {
  dat <- data.frame(age = c(-1, -0.2, 0.1), stringsAsFactors = FALSE)
  out <- prep_assign_age_bins(
    dat,
    age_col = "age",
    bins = "neonatal",
    negative_age_strategy = "na"
  )
  expect_true(is.na(out$Age_bin[1]))
  expect_true(is.na(out$Age_bin[2]))
  expect_false(is.na(out$Age_bin[3]))
})

test_that("prep_fill_age: derives age from DOB-only scenario fixture", {
  dat <- cat1_age_from_dob_only()
  expect_warning(
    dat <- prep_coerce_dates(dat, cols = "dob", table_label = "cat1_age_from_dob_only"),
    "appear encrypted/undecodable"
  )
  expect_warning(
    dat <- prep_coerce_dates(dat, cols = "culture_date", table_label = "cat1_age_from_dob_only"),
    "appear encrypted/undecodable"
  )
  out <- prep_fill_age(dat, age_col = "age", dob_col = "dob", date_col = "culture_date")
  expect_true(any(!is.na(out$age)))
  expect_true(any(is.na(out$age))) # invalid/encoded dates remain missing
  expect_true("age_method" %in% names(out))
})

test_that("prep_assign_age_bins: error handling", {
  expect_error(prep_assign_age_bins(data.frame(x = 1), age_col = "Age"), "not found")
  expect_error(prep_assign_age_bins(data.frame(Age = 10), bins = "invalid_mode"), "Unknown age bin type")
  expect_error(
    prep_assign_age_bins(data.frame(Age = 1), age_months_col = "missing_col"),
    "not found"
  )
  expect_error(
    prep_assign_age_bins(data.frame(Age = 1), age_days_col = "missing_col"),
    "not found"
  )
})

# find_extdata_file ------------------------------------------------------------

test_that("find_extdata_file: valid behavior", {
  p <- anumaan:::find_extdata_file("organisms.csv")
  expect_true(is.character(p))
  expect_true(nchar(p) > 0)
  expect_true(file.exists(p))
})

test_that("find_extdata_file: edge cases", {
  p <- anumaan:::find_extdata_file("common_commensals.csv")
  expect_true(file.exists(p))
})

test_that("find_extdata_file: missingness", {
  p <- anumaan:::find_extdata_file("")
  expect_true(is.character(p))
})

test_that("find_extdata_file: error handling", {
  p <- anumaan:::find_extdata_file("nonexistent_centre_name_zz.csv")
  expect_equal(p, "")
})
