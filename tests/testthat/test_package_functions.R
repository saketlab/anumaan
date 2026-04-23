# anumaan Package - Core Function Tests

test_that("prep_get_contaminant_list returns expected contaminants", {
  contaminants <- prep_get_contaminant_list(
    syndrome = "Bloodstream infections",
    return_all = FALSE
  )

  expect_type(contaminants, "list")
  expect_true(length(contaminants$names) > 0)
})

test_that("prep_is_contaminant identifies known contaminants", {
  test_organisms <- c(
    "staphylococcus epidermidis",
    "escherichia coli",
    "klebsiella pneumoniae"
  )
  results <- prep_is_contaminant(
    organism_name = test_organisms,
    syndrome = "Bloodstream infections"
  )
  expect_length(results, 3)
  expect_type(results, "logical")
})

test_that("prep_fill_age calculates age from DOB", {
  test_data <- data.frame(
    DOB = as.Date("1990-01-01"),
    date_of_culture = as.Date("2020-06-15"),
    stringsAsFactors = FALSE
  )
  result <- prep_fill_age(test_data)
  expect_true("Age" %in% names(result))
  expect_true(result$Age > 30 && result$Age < 31)
})

test_that("prep_assign_age_bins creates age bins", {
  test_data <- data.frame(Age = c(2, 15, 35, 65, 90))
  result <- prep_assign_age_bins(test_data, bins = "GBD_standard")
  expect_true("Age_bin" %in% names(result))
  expect_s3_class(result$Age_bin, "factor")
  expect_equal(sum(!is.na(result$Age_bin)), 5)
})

test_that("prep_extract_genus extracts first word", {
  test_data <- data.frame(
    organism_normalized = c("escherichia coli", "klebsiella pneumoniae")
  )
  result <- prep_extract_genus(test_data)
  expect_equal(result$org_genus, c("escherichia", "klebsiella"))
})

test_that("prep_extract_species extracts second word", {
  test_data <- data.frame(
    organism_normalized = c("escherichia coli", "klebsiella pneumoniae")
  )
  result <- prep_extract_species(test_data)
  expect_equal(result$org_species, c("coli", "pneumoniae"))
})

test_that("prep_derive_los_from_dates computes length of stay", {
  test_data <- data.frame(
    date_of_admission = as.Date(c("2020-01-01", "2020-02-01")),
    date_of_final_outcome = as.Date(c("2020-01-10", "2020-02-05"))
  )
  result <- prep_derive_los_from_dates(
    test_data,
    admission_col = "date_of_admission",
    outcome_col   = "date_of_final_outcome"
  )
  expect_equal(result$los_days, c(9, 4))
})
