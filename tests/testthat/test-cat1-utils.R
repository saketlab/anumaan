# Category 1: zzz_utils_internal.R

# parse_age_bin_labels ---------------------------------------------------------

test_that("parse_age_bin_labels: valid behavior", {
  out <- anumaan:::parse_age_bin_labels(c("<1", "1-5", "5-10", "85+"))
  expect_type(out, "list")
  expect_equal(out$labels, c("<1", "1-5", "5-10", "85+"))
  expect_equal(out$breaks, c(-Inf, 1, 5, 10, 85, Inf))
})

test_that("parse_age_bin_labels: edge cases", {
  out <- anumaan:::parse_age_bin_labels(c("0-0.08", "0.08-1", "18+"))
  expect_true(is.numeric(out$breaks))
  expect_true(is.infinite(tail(out$breaks, 1)))
})

test_that("parse_age_bin_labels: missingness", {
  expect_error(anumaan:::parse_age_bin_labels(c("", NA_character_)), "Cannot parse age bin label")
})

test_that("parse_age_bin_labels: error handling", {
  expect_error(anumaan:::parse_age_bin_labels("nonsense"), "Cannot parse age bin label")
})

# round_to_sum -----------------------------------------------------------------

test_that("round_to_sum: valid behavior", {
  x <- c(3.3, 3.3, 3.4)
  out <- round_to_sum(x, target = 10)
  expect_length(out, 3)
  expect_equal(sum(out), 10)
})

test_that("round_to_sum: edge cases", {
  expect_equal(round_to_sum(5.7, target = 6), 6)
  expect_equal(sum(round_to_sum(c(0, 0, 0), target = 0)), 0)
})

test_that("round_to_sum: missingness", {
  out <- round_to_sum(c(1.2, NA_real_), target = 1)
  expect_true(any(is.na(out)))
})

test_that("round_to_sum: error handling", {
  expect_error(round_to_sum(c("a", "b"), target = 1), NA)
})

# shorten_drug_class -----------------------------------------------------------

test_that("shorten_drug_class: valid behavior", {
  out <- shorten_drug_class(c("Third-generation-cephalosporins", "Beta-lactam/beta-lactamase-inhibitor"))
  expect_equal(out, c("3GC", "BL/BLI"))
})

test_that("shorten_drug_class: edge cases", {
  out <- shorten_drug_class(c("Third-generation-cephalosporins_R", "Anti-pseudomonal-penicillins"))
  expect_equal(out, c("3GC", "Anti-pseudo-PCN"))
})

test_that("shorten_drug_class: missingness", {
  out <- shorten_drug_class(c(NA_character_, ""))
  expect_true(is.na(out[1]))
  expect_equal(out[2], "")
})

test_that("shorten_drug_class: error handling", {
  expect_error(shorten_drug_class(1:3), NA)
})

# normalize_join ---------------------------------------------------------------

test_that("normalize_join: valid behavior", {
  out <- anumaan:::normalize_join(c("  Blood-Culture ", "CSF (sample)"))
  expect_equal(out, c("bloodculture", "csf sample"))
})

test_that("normalize_join: edge cases", {
  out <- anumaan:::normalize_join(c("E. coli", "a   b   c", "###Name###"))
  expect_equal(out, c("e coli", "a b c", "name"))
})

test_that("normalize_join: missingness", {
  out <- anumaan:::normalize_join(c(NA_character_, ""))
  expect_true(is.na(out[1]))
  expect_equal(out[2], "")
})

test_that("normalize_join: error handling", {
  expect_error(anumaan:::normalize_join(list("a")), NA)
})
