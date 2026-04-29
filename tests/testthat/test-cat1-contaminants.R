# Category 1: prep_contaminants.R (target helper subset)

# prep_get_contaminant_list ----------------------------------------------------

test_that("prep_get_contaminant_list: valid behavior", {
  out <- prep_get_contaminant_list(syndrome = "Bloodstream infections")
  expect_type(out, "list")
  expect_true(all(c("names", "patterns") %in% names(out)))
  expect_true(length(out$names) > 0)
})

test_that("prep_get_contaminant_list: edge cases", {
  all_out <- prep_get_contaminant_list(return_all = TRUE)
  bsi_out <- prep_get_contaminant_list(syndrome = "Bloodstream infections")
  expect_gte(length(all_out$names), length(bsi_out$names))
  expect_gt(length(all_out$patterns), 0)
})

test_that("prep_get_contaminant_list: missingness", {
  out <- prep_get_contaminant_list(syndrome = NA_character_)
  expect_type(out, "list")
  expect_true(length(out$names) >= 0)
})

test_that("prep_get_contaminant_list: error handling", {
  out <- prep_get_contaminant_list(syndrome = "NONEXISTENT_CENTRE_NAME", return_all = FALSE)
  expect_length(out$names, 0)
  expect_length(out$patterns, 0)
})

# prep_is_contaminant ----------------------------------------------------------

test_that("prep_is_contaminant: valid behavior", {
  out <- prep_is_contaminant(
    organism_name = c("Staphylococcus epidermidis", "Escherichia coli"),
    syndrome = "Bloodstream infections"
  )
  expect_equal(out, c(TRUE, FALSE))
})

test_that("prep_is_contaminant: edge cases", {
  out <- prep_is_contaminant(
    organism_name = c("S. EPIDERMIDIS", "staphylococcus epidermidis", "Staphylococcus epidermidis"),
    syndrome = "Bloodstream infections"
  )
  expect_true(all(out))
})

test_that("prep_is_contaminant: missingness", {
  out <- prep_is_contaminant(
    organism_name = c(NA, "", "   "),
    syndrome = "Bloodstream infections"
  )
  expect_equal(out, c(FALSE, FALSE, FALSE))
})

test_that("prep_is_contaminant: error handling", {
  out <- prep_is_contaminant(
    organism_name = c("Staphylococcus epidermidis", "Escherichia coli"),
    syndrome = "SYNDROME_THAT_DOES_NOT_EXIST"
  )
  expect_equal(out, c(FALSE, FALSE))
})
