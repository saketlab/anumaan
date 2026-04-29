# Category 1: prep_standardize_specimens.R

# prep_standardize_specimens ---------------------------------------------------

test_that("prep_standardize_specimens: valid behavior", {
  dat <- data.frame(specimen_type = c("Blood culture", "midstream urine c/s"), stringsAsFactors = FALSE)
  out <- prep_standardize_specimens(dat)
  expect_equal(out$specimen_normalized, c("Blood", "Urine"))
  expect_true(all(c("sample_category", "sterile_classification") %in% names(out)))
})

test_that("prep_standardize_specimens: edge cases", {
  dat <- data.frame(
    specimen_type = c("blood culture", "blood culture", "Urine culture / sensitivity", "unknown specimen"),
    stringsAsFactors = FALSE
  )
  out <- prep_standardize_specimens(dat)
  expect_equal(out$specimen_normalized[1], out$specimen_normalized[2])
  expect_equal(out$specimen_normalized[3], "Urine")
  expect_true(is.na(out$specimen_normalized[4]))
})

test_that("prep_standardize_specimens: edge cases - typo spellings", {
  dat <- data.frame(
    specimen_type = c("blood culure", "urne c/s", "csff"),
    stringsAsFactors = FALSE
  )
  out <- prep_standardize_specimens(dat)
  expect_equal(out$specimen_normalized[1], "Blood")
  expect_equal(out$specimen_normalized[2], "Urine")
  expect_true(grepl("^CSF", out$specimen_normalized[3]))
})

test_that("prep_standardize_specimens: missingness", {
  dat <- data.frame(specimen_type = c(NA, "", "   "), stringsAsFactors = FALSE)
  out <- prep_standardize_specimens(dat)
  expect_true(all(is.na(out$specimen_normalized)))
})

test_that("prep_standardize_specimens: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_standardize_specimens(dat, specimen_col = "specimen_type"), "not found")
  expect_false("specimen_normalized" %in% names(out))
})
