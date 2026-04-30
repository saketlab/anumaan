# Category 2: prep_hai_cai.R

# prep_derive_hai_cai ----------------------------------------------------------

test_that("prep_derive_hai_cai: valid behavior", {
  dat <- cat2_hai_gap_logic()[c(1, 4), c("date_of_admission", "date_of_culture", "infection_type"), drop = FALSE]
  out <- prep_derive_hai_cai(dat)
  expect_equal(out$infection_type, c("CAI", "HAI"))
})

test_that("prep_derive_hai_cai: edge cases", {
  dat <- cat2_hai_gap_logic()[c(3, 2), c("date_of_admission", "date_of_culture"), drop = FALSE]
  dat$infection_type <- c("CAI", "HAI")
  out <- prep_derive_hai_cai(dat, overwrite = FALSE)
  expect_equal(out$infection_type, c("CAI", "HAI"))
})

test_that("prep_derive_hai_cai: missingness", {
  dat <- cat2_hai_gap_logic()[c(5, 6), c("date_of_admission", "date_of_culture"), drop = FALSE]
  out <- prep_derive_hai_cai(dat)
  expect_true(any(is.na(out$infection_type)))
})

test_that("prep_derive_hai_cai: error handling", {
  dat <- data.frame(x = 1)
  out <- prep_derive_hai_cai(dat)
  expect_true("infection_type" %in% names(out))
})

# prep_flag_hai_inferred -------------------------------------------------------

test_that("prep_flag_hai_inferred: valid behavior", {
  dat <- cat2_hai_observed_labels()[c(1, 2), c("infection_type", "infection_type_method"), drop = FALSE]
  dat$infection_type_method[1] <- "inferred_2day_cutoff"
  out <- prep_flag_hai_inferred(dat)
  expect_equal(out$infection_type_src, c("inferred", "observed"))
})

test_that("prep_flag_hai_inferred: edge cases", {
  dat <- cat2_hai_observed_labels()[c(2, 5), c("infection_type", "infection_type_method"), drop = FALSE]
  dat$infection_type_method[] <- NA_character_
  out <- prep_flag_hai_inferred(dat)
  expect_equal(out$infection_type_src, c("observed", "unknown"))
})

test_that("prep_flag_hai_inferred: missingness", {
  dat <- cat2_hai_observed_labels()[c(3, 4), "infection_type", drop = FALSE]
  out <- prep_flag_hai_inferred(dat)
  expect_true(all(out$infection_type_src == "unknown"))
})

test_that("prep_flag_hai_inferred: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_flag_hai_inferred(dat), "not found")
  expect_true(all(is.na(out$infection_type_src)))
})

# prep_reconcile_hai_observed_inferred -----------------------------------------

test_that("prep_reconcile_hai_observed_inferred: valid behavior", {
  dat <- cat2_hai_observed_vs_inferred()[c(1, 3), , drop = FALSE]
  expect_warning(out <- prep_reconcile_hai_observed_inferred(dat), "discordant")
  expect_equal(out$infection_type, c("HAI", "HAI"))
})

test_that("prep_reconcile_hai_observed_inferred: edge cases", {
  dat <- cat2_hai_observed_vs_inferred()[4, , drop = FALSE]
  out <- prep_reconcile_hai_observed_inferred(dat)
  expect_false(out$hai_discordant)
})

test_that("prep_reconcile_hai_observed_inferred: missingness", {
  dat <- cat2_hai_observed_vs_inferred()[c(3, 4), , drop = FALSE]
  dat$infection_type_observed[2] <- NA_character_
  dat$infection_type_inferred[2] <- NA_character_
  out <- prep_reconcile_hai_observed_inferred(dat)
  expect_equal(out$infection_type[1], "HAI")
})

test_that("prep_reconcile_hai_observed_inferred: error handling", {
  dat <- data.frame(x = 1)
  out <- prep_reconcile_hai_observed_inferred(dat)
  expect_identical(out, dat)
})

# prep_derive_icu_flag ---------------------------------------------------------

test_that("prep_derive_icu_flag: valid behavior", {
  dat <- cat2_icu_location_variants()[c(1, 2), , drop = FALSE]
  out <- prep_derive_icu_flag(dat)
  expect_equal(out$icu_flag, c(TRUE, FALSE))
})

test_that("prep_derive_icu_flag: edge cases", {
  dat <- cat2_icu_location_variants()[c(3, 4), , drop = FALSE]
  out <- prep_derive_icu_flag(dat)
  expect_equal(out$icu_flag, c(TRUE, TRUE))
})

test_that("prep_derive_icu_flag: missingness", {
  dat <- cat2_icu_location_variants()[7, , drop = FALSE]
  dat <- dat[rep(1, 2), , drop = FALSE]
  out <- prep_derive_icu_flag(dat)
  expect_false(any(out$icu_flag))
})

test_that("prep_derive_icu_flag: error handling", {
  out <- prep_derive_icu_flag(data.frame(x = 1), ward_col = "ward_icu", dept_col = "hospital_department")
  expect_true("icu_flag" %in% names(out))
})
