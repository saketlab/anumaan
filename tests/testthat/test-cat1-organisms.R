# Category 1: prep_standardize_organisms.R (target helper subset)

# prep_standardize_organisms ---------------------------------------------------

test_that("prep_standardize_organisms: clinical text variants map as expected", {
  dat <- data.frame(
    organism_name = c(
      "Non fermenting Gram negative bacilli",
      "No growth in culture.",
      "Bacillus spp. grown.",
      "CONS (Coagulase Negative Staphylococci)",
      "MRSA(Methicillin resistant staphylococcus aureus)"
    ),
    stringsAsFactors = FALSE
  )
  out <- prep_standardize_organisms(dat, organism_col = "organism_name")

  expect_equal(out$organism_normalized[1], "non-fermenting gram-negative bacilli")
  expect_true(is.na(out$organism_normalized[2]))
  expect_equal(out$organism_normalized[3], "bacillus spp.")
  expect_equal(out$organism_normalized[4], "coagulase-negative staphylococci")
  expect_equal(out$organism_normalized[5], "staphylococcus aureus")
  expect_equal(out$is_MRSA[5], 1L)
})

test_that("prep_standardize_organisms + prep_is_contaminant: CONS and Bacillus flagged", {
  dat <- data.frame(
    organism_name = c(
      "CONS (Coagulase Negative Staphylococci)",
      "Bacillus spp. grown.",
      "MRSA(Methicillin resistant staphylococcus aureus)"
    ),
    stringsAsFactors = FALSE
  )
  out <- prep_standardize_organisms(dat, organism_col = "organism_name")
  is_contam <- prep_is_contaminant(
    organism_name = out$organism_normalized,
    syndrome = "Bloodstream infections"
  )

  expect_equal(is_contam, c(TRUE, TRUE, FALSE))
})

# prep_extract_genus -----------------------------------------------------------

test_that("prep_extract_genus: valid behavior", {
  dat <- data.frame(organism_normalized = c("escherichia coli", "klebsiella pneumoniae"), stringsAsFactors = FALSE)
  out <- prep_extract_genus(dat)
  expect_equal(out$org_genus, c("escherichia", "klebsiella"))
})

test_that("prep_extract_genus: edge cases", {
  dat <- data.frame(organism_normalized = c("staphylococcus aureus", "candida spp."), stringsAsFactors = FALSE)
  out <- prep_extract_genus(dat)
  expect_equal(out$org_genus, c("staphylococcus", "candida"))
})

test_that("prep_extract_genus: missingness", {
  dat <- data.frame(organism_normalized = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_extract_genus(dat)
  expect_true(all(is.na(out$org_genus)))
})

test_that("prep_extract_genus: error handling", {
  expect_error(prep_extract_genus(data.frame(x = 1), organism_col = "organism_normalized"), "not found")
})

# prep_extract_species ---------------------------------------------------------

test_that("prep_extract_species: valid behavior", {
  dat <- data.frame(organism_normalized = c("escherichia coli", "klebsiella pneumoniae"), stringsAsFactors = FALSE)
  out <- prep_extract_species(dat)
  expect_equal(out$org_species, c("coli", "pneumoniae"))
})

test_that("prep_extract_species: edge cases", {
  dat <- data.frame(organism_normalized = c("candida spp.", "enterococcus spp"), stringsAsFactors = FALSE)
  out <- prep_extract_species(dat)
  expect_equal(out$org_species, c("species", "species"))
})

test_that("prep_extract_species: missingness", {
  dat <- data.frame(organism_normalized = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_extract_species(dat)
  expect_true(all(is.na(out$org_species)))
})

test_that("prep_extract_species: error handling", {
  expect_error(prep_extract_species(data.frame(x = 1), organism_col = "organism_normalized"), "not found")
})

# prep_normalize_sp_variants ---------------------------------------------------

test_that("prep_normalize_sp_variants: valid behavior", {
  out <- anumaan:::prep_normalize_sp_variants(c("Acinetobacter sp", "E. coli"))
  expect_equal(out, c("Acinetobacter spp.", "Escherichia coli"))
})

test_that("prep_normalize_sp_variants: edge cases", {
  out <- anumaan:::prep_normalize_sp_variants(c("K. pneumoniae", "P. aeruginosa", "non fermanting bacilli"))
  expect_equal(out[1], "Klebsiella pneumoniae")
  expect_equal(out[2], "Pseudomonas aeruginosa")
  expect_true(grepl("Non-fermenting", out[3]))
})

test_that("prep_normalize_sp_variants: missingness", {
  out <- anumaan:::prep_normalize_sp_variants(c(NA_character_, ""))
  expect_true(is.na(out[1]))
  expect_equal(out[2], "")
})

test_that("prep_normalize_sp_variants: error handling", {
  bad_input <- function() "E. coli"
  expect_error(anumaan:::prep_normalize_sp_variants(bad_input), "cannot coerce")
})

# get_organism_taxonomy --------------------------------------------------------

test_that("get_organism_taxonomy: valid behavior", {
  tx <- anumaan:::get_organism_taxonomy()
  expect_s3_class(tx, "data.frame")
  expect_true(all(c("organism_name", "org_group") %in% names(tx)))
})

test_that("get_organism_taxonomy: edge cases", {
  tx <- anumaan:::get_organism_taxonomy()
  expect_true(nrow(tx) > 0)
  expect_true(any(tolower(tx$organism_name) == "escherichia coli"))
})

test_that("get_organism_taxonomy: missingness", {
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(organism_name = c(NA, "Escherichia coli"), organism_group = c("Other", "Enterobacterales")),
    tmp,
    row.names = FALSE
  )
  testthat::with_mocked_bindings(
    {
      tx <- anumaan:::get_organism_taxonomy()
      expect_true(is.na(tx$organism_name[1]))
      expect_equal(tx$org_group[2], "Enterobacterales")
    },
    find_extdata_file = function(filename) tmp,
    .package = "anumaan"
  )
})

test_that("get_organism_taxonomy: error handling", {
  testthat::with_mocked_bindings(
    {
      expect_warning(tx <- anumaan:::get_organism_taxonomy(), "not found")
      expect_equal(nrow(tx), 0)
      expect_true(all(c("organism_name", "org_group") %in% names(tx)))
    },
    find_extdata_file = function(filename) "",
    .package = "anumaan"
  )
})
