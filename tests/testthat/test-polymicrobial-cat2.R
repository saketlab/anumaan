# Category 2: prep_polymicrobial.R

# prep_flag_polymicrobial ------------------------------------------------------

test_that("prep_flag_polymicrobial: valid behavior", {
  dat <- data.frame(patient_id = c("p1", "p1", "p2"), organism_normalized = c("e coli", "k pn", "e coli"), stringsAsFactors = FALSE)
  out <- prep_flag_polymicrobial(dat)
  expect_true(all(c("n_organisms", "is_polymicrobial") %in% names(out)))
  expect_equal(out$is_polymicrobial, c(1L, 1L, 0L))
})

test_that("prep_flag_polymicrobial: edge cases", {
  dat <- data.frame(patient_id = c("p1", "p1", "p1"), organism_normalized = c("e coli", "e coli", "k pn"), facility = c("A", "A", "B"), stringsAsFactors = FALSE)
  out <- prep_flag_polymicrobial(dat, facility_col = "facility")
  expect_true(all(out$is_polymicrobial %in% c(0L, 1L)))
})

test_that("prep_flag_polymicrobial: missingness", {
  dat <- data.frame(patient_id = c("p1", NA), organism_normalized = c("e coli", NA), stringsAsFactors = FALSE)
  out <- prep_flag_polymicrobial(dat)
  expect_true("is_polymicrobial" %in% names(out))
})

test_that("prep_flag_polymicrobial: error handling", {
  expect_error(prep_flag_polymicrobial(data.frame(x = 1)), "Missing required columns")
  dat <- data.frame(patient_id = "p1", organism_normalized = "e coli", stringsAsFactors = FALSE)
  expect_error(prep_flag_polymicrobial(dat, facility_name = "NONEXISTENT_CENTRE"), "facility_col must be supplied")
})

# prep_compute_poly_weights ----------------------------------------------------

test_that("prep_compute_poly_weights: valid behavior", {
  dat <- data.frame(episode_id = c("e1", "e1", "e2"), organism_normalized = c("e coli", "k pn", "e coli"), is_polymicrobial = c(1, 1, 0), stringsAsFactors = FALSE)
  out <- prep_compute_poly_weights(dat, method = "equal")
  expect_true(all(c("poly_weight", "weight_method", "weight_confidence") %in% names(out)))
  expect_true(all(out$poly_weight > 0))
})

test_that("prep_compute_poly_weights: edge cases", {
  dat <- data.frame(episode_id = c("e1", "e1"), organism_normalized = c("e coli", "k pn"), is_polymicrobial = c(1, 1), stringsAsFactors = FALSE)
  out <- prep_compute_poly_weights(dat, method = "manual", weight_map = c("e coli" = 0.7, "k pn" = 0.3))
  expect_true(all(out$poly_weight > 0))
})

test_that("prep_compute_poly_weights: missingness", {
  dat <- data.frame(episode_id = c("e1", "e1"), organism_normalized = c("unknown1", "unknown2"), is_polymicrobial = c(1, 1), stringsAsFactors = FALSE)
  out <- prep_compute_poly_weights(dat, method = "monomicrobial_proportion")
  expect_true(all(out$poly_weight > 0))
})

test_that("prep_compute_poly_weights: error handling", {
  expect_error(prep_compute_poly_weights(data.frame(x = 1)), "Missing required columns")
  dat <- data.frame(episode_id = "e1", organism_normalized = "e coli", is_polymicrobial = 1)
  expect_error(prep_compute_poly_weights(dat, method = "manual"), "weight_map")
})

# prep_split_poly_episode ------------------------------------------------------

test_that("prep_split_poly_episode: valid behavior", {
  dat <- data.frame(is_polymicrobial = c(1, 0, 1), x = 1:3)
  out <- prep_split_poly_episode(dat, strategy = "exclude")
  expect_equal(nrow(out), 1)
})

test_that("prep_split_poly_episode: edge cases", {
  dat <- data.frame(is_polymicrobial = c(1, 0, 1), x = 1:3)
  out <- prep_split_poly_episode(dat, strategy = "separate")
  poly <- attr(out, "poly_data")
  expect_s3_class(poly, "data.frame")
  expect_equal(nrow(poly), 2)
})

test_that("prep_split_poly_episode: missingness", {
  dat <- data.frame(is_polymicrobial = c(NA, 0, 1), x = 1:3)
  out <- prep_split_poly_episode(dat, strategy = "fractional")
  expect_equal(nrow(out), 3)
})

test_that("prep_split_poly_episode: error handling", {
  dat <- data.frame(x = 1:2)
  expect_warning(out <- prep_split_poly_episode(dat), "not found")
  expect_identical(out, dat)
})
