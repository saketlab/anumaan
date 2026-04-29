# Category 1: prep_standardize_antibiotics.R (target helper subset)

# prep_classify_aware ----------------------------------------------------------

test_that("prep_classify_aware: valid behavior", {
  dat <- data.frame(antibiotic_normalized = c("Amikacin", "Ampicillin"), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = c("Amikacin", "Ampicillin"), Category = c("Access", "Watch"), stringsAsFactors = FALSE)
  out <- prep_classify_aware(dat, who_table = who)
  expect_equal(out$aware_category, c("Access", "Watch"))
})

test_that("prep_classify_aware: edge cases", {
  dat <- data.frame(antibiotic_normalized = c("Amikacin", "Amikacin", "UnknownDrug"), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = "Amikacin", Category = "Access", stringsAsFactors = FALSE)
  out <- prep_classify_aware(dat, who_table = who)
  expect_equal(out$aware_category[1], "Access")
  expect_equal(out$aware_category[2], "Access")
  expect_true(is.na(out$aware_category[3]))
})

test_that("prep_classify_aware: missingness", {
  dat <- data.frame(antibiotic_normalized = c(NA, ""), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = "Amikacin", Category = "Access", stringsAsFactors = FALSE)
  out <- prep_classify_aware(dat, who_table = who)
  expect_true(all(is.na(out$aware_category)))
})

test_that("prep_classify_aware: error handling", {
  dat <- data.frame(x = "Amikacin")
  expect_error(prep_classify_aware(dat, antibiotic_col = "antibiotic_normalized"), "not found")
  dat2 <- data.frame(antibiotic_normalized = "Amikacin")
  expect_warning(out <- prep_classify_aware(dat2, who_table = NULL), "not provided")
  expect_equal(out$aware_source, "needs_who_table")
})

# prep_decode_antibiotic_code --------------------------------------------------

test_that("prep_decode_antibiotic_code: valid behavior", {
  map <- tempfile(fileext = ".csv")
  readr::write_csv(data.frame(code = c("AMK", "AMP"), antibiotic_name_full = c("Amikacin", "Ampicillin")), map)
  dat <- data.frame(antibiotic_name_raw = c("amk", "AMP"), stringsAsFactors = FALSE)
  out <- prep_decode_antibiotic_code(dat, map_path = map)
  expect_equal(out$antibiotic_name_std, c("Amikacin", "Ampicillin"))
})

test_that("prep_decode_antibiotic_code: edge cases", {
  map <- tempfile(fileext = ".csv")
  readr::write_csv(data.frame(code = "AMK", antibiotic_name_full = "Amikacin"), map)
  dat <- data.frame(antibiotic_name_raw = c("AMK", "AMK", "XYZ"), stringsAsFactors = FALSE)
  out <- prep_decode_antibiotic_code(dat, map_path = map)
  expect_equal(out$antibiotic_name_std[1:2], c("Amikacin", "Amikacin"))
  expect_equal(out$antibiotic_name_std[3], "XYZ")
})

test_that("prep_decode_antibiotic_code: missingness", {
  map <- tempfile(fileext = ".csv")
  readr::write_csv(data.frame(code = "AMK", antibiotic_name_full = "Amikacin"), map)
  dat <- data.frame(antibiotic_name_raw = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_decode_antibiotic_code(dat, map_path = map)
  expect_true(is.na(out$antibiotic_name_std[1]))
  expect_equal(out$antibiotic_name_std[2], "")
})

test_that("prep_decode_antibiotic_code: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_decode_antibiotic_code(dat, code_col = "antibiotic_name_raw"), "not found")
  expect_true("antibiotic_name_std" %in% names(out))

  bad_map <- tempfile(fileext = ".csv")
  readr::write_csv(data.frame(code = "AMK"), bad_map)
  dat2 <- data.frame(antibiotic_name_raw = "AMK")
  expect_error(prep_decode_antibiotic_code(dat2, map_path = bad_map), "must have columns")
})
