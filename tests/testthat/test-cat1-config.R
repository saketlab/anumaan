# Category 1: config.R

# amr_config -----------------------------------------------------------------

test_that("amr_config: valid behavior", {
  cfg <- amr_config(hai_cutoff = 3, strict_validation = TRUE)
  expect_s3_class(cfg, "amr_config")
  expect_equal(cfg$hai_cutoff, 3)
  expect_true(cfg$strict_validation)
})

test_that("amr_config: edge cases", {
  cfg <- amr_config(age_bins = "geriatric", mdr_definition = 1)
  expect_true(is.character(cfg$age_bins))
  expect_true("90+" %in% cfg$age_bins)
  expect_equal(cfg$mdr_definition, 1)
})

test_that("amr_config: missingness", {
  cfg <- amr_config(column_mappings = NULL)
  expect_true(is.list(cfg$column_mappings))
  expect_true(length(cfg$date_columns) > 0)
})

test_that("amr_config: error handling", {
  expect_error(amr_config(age_bins = "wrong_mode_name"), "Unknown age bin type")
})

# validate_config --------------------------------------------------------------

test_that("validate_config: valid behavior", {
  expect_true(validate_config(amr_config()))
})

test_that("validate_config: edge cases", {
  cfg <- amr_config(hai_cutoff = 0, event_gap_days = 1, mortality_window = 1, mdr_definition = 1)
  expect_true(validate_config(cfg))
})

test_that("validate_config: missingness", {
  cfg <- amr_config()
  cfg$age_bins <- c("<1")
  expect_error(validate_config(cfg), "age_bins")
})

test_that("validate_config: error handling", {
  expect_error(validate_config(list()), "amr_config object")
  expect_error(validate_config(amr_config(contaminant_method = "bad_mode")), "contaminant_method")
  expect_error(validate_config(amr_config(mdr_definition = 0)), "mdr_definition")
})
