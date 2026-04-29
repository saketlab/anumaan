# Category 1: prep_types.R

# prep_standardize_sex ---------------------------------------------------------

test_that("prep_standardize_sex: valid behavior", {
  dat <- data.frame(gender = c("M", "female", "man", "WOMAN"), stringsAsFactors = FALSE)
  out <- suppressWarnings(prep_standardize_sex(dat, col = "gender"))
  expect_equal(out$gender, c("M", "F", "M", "F"))
})

test_that("prep_standardize_sex: edge cases", {
  dat <- data.frame(gender = c(" M ", "f", "FEMALE", "MALE", "unknown"), stringsAsFactors = FALSE)
  expect_warning(out <- prep_standardize_sex(dat, col = "gender"), "could not be standardized")
  expect_equal(out$gender, c(NA_character_, "F", "F", "M", NA_character_))
})

test_that("prep_standardize_sex: missingness", {
  dat <- data.frame(gender = c(NA, "", "  "), stringsAsFactors = FALSE)
  expect_warning(out <- prep_standardize_sex(dat, col = "gender"), "could not be standardized")
  expect_true(all(is.na(out$gender)))
})

test_that("prep_standardize_sex: error handling", {
  expect_error(prep_standardize_sex(data.frame(sex = "M"), col = "gender"), "not found")
})

# prep_standardize_outcome -----------------------------------------------------

test_that("prep_standardize_outcome: valid behavior", {
  dat <- data.frame(final_outcome = c("alive", "expired"), stringsAsFactors = FALSE)
  expect_warning(out <- prep_standardize_outcome(dat), "prep_standardize_final_outcome")
  expect_equal(out$outcome_std, c("Survived", "Died"))
})

test_that("prep_standardize_outcome: edge cases", {
  dat <- data.frame(final_outcome = c("Survived/Discharged", "Died/Expired"), stringsAsFactors = FALSE)
  expect_warning(out <- prep_standardize_outcome(dat), "prep_standardize_final_outcome")
  expect_equal(out$outcome_std, c("Survived", "Died"))
})

test_that("prep_standardize_outcome: missingness", {
  dat <- data.frame(final_outcome = c(NA, "", "LAMA"), stringsAsFactors = FALSE)
  expect_warning(out <- prep_standardize_outcome(dat), "prep_standardize_final_outcome")
  expect_true(all(is.na(out$outcome_std)))
})

test_that("prep_standardize_outcome: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_standardize_outcome(dat, col = "final_outcome"), "not found")
  expect_true(all(is.na(out$outcome_std)))
})

# prep_standardize_final_outcome -----------------------------------------------

test_that("prep_standardize_final_outcome: valid behavior", {
  dat <- cat1_min_outcome()
  out <- prep_standardize_final_outcome(dat)
  expect_equal(out$outcome_std[1:2], c("Survived", "Died"))
})

test_that("prep_standardize_final_outcome: edge cases", {
  dat <- data.frame(final_outcome = c(" Discharged Alive ", "DEAD", "Left Against Medical Advice"), stringsAsFactors = FALSE)
  out <- prep_standardize_final_outcome(dat)
  expect_equal(out$outcome_std, c("Survived", "Died", NA_character_))
})

test_that("prep_standardize_final_outcome: missingness", {
  dat <- data.frame(final_outcome = c(NA, "", "NA"), stringsAsFactors = FALSE)
  out <- prep_standardize_final_outcome(dat)
  expect_true(all(is.na(out$outcome_std)))
})

test_that("prep_standardize_final_outcome: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_standardize_final_outcome(dat, col = "final_outcome"), "not found")
  expect_true(all(is.na(out$outcome_std)))
})

# prep_standardize_infection_type ----------------------------------------------

test_that("prep_standardize_infection_type: valid behavior", {
  out <- prep_standardize_infection_type(cat1_min_infection())
  expect_equal(out$infection_type[1:2], c("HAI", "CAI"))
})

test_that("prep_standardize_infection_type: edge cases", {
  dat <- data.frame(infection_type = c("Hospital-Acquired", "Healthcare Associated", "Community-Acquired"), stringsAsFactors = FALSE)
  out <- prep_standardize_infection_type(dat)
  expect_equal(out$infection_type, c("HAI", "HAI", "CAI"))
})

test_that("prep_standardize_infection_type: missingness", {
  dat <- data.frame(infection_type = c(NA, "", "  "), stringsAsFactors = FALSE)
  out <- prep_standardize_infection_type(dat)
  expect_true(all(is.na(out$infection_type)))
})

test_that("prep_standardize_infection_type: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_standardize_infection_type(dat, col = "infection_type"), "not found")
  expect_false("infection_type" %in% names(out))
})
