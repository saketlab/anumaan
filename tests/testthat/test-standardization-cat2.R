# Category 2: organism/antibiotic standardization

# prep_standardize_organisms ---------------------------------------------------

test_that("prep_standardize_organisms: valid behavior", {
  dat <- cat2_standardization_organism_abx()[c(1, 2), "organism_name", drop = FALSE]
  out <- prep_standardize_organisms(dat, organism_col = "organism_name")
  expect_true(all(c("organism_normalized", "organism_group") %in% names(out)))
})

test_that("prep_standardize_organisms: edge cases", {
  dat <- cat2_standardization_organism_abx()[c(3, 4, 6), "organism_name", drop = FALSE]
  out <- prep_standardize_organisms(dat, organism_col = "organism_name")
  expect_true(any(grepl("coagulase-negative", out$organism_normalized, ignore.case = TRUE)))
})

test_that("prep_standardize_organisms: missingness", {
  dat <- data.frame(organism_name = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_standardize_organisms(dat, organism_col = "organism_name")
  expect_true(all(is.na(out$organism_normalized)))
})

test_that("prep_standardize_organisms: error handling", {
  expect_error(prep_standardize_organisms(data.frame(x = 1), organism_col = "organism_name"), "not found")
})

# prep_assign_organism_group ---------------------------------------------------

test_that("prep_assign_organism_group: valid behavior", {
  dat <- data.frame(
    organism_normalized = c("Escherichia coli", "Unknown environmental isolate"),
    stringsAsFactors = FALSE
  )
  out <- prep_assign_organism_group(dat)
  expect_true("org_group" %in% names(out))
})

test_that("prep_assign_organism_group: edge cases", {
  dat <- data.frame(organism_normalized = c("Escherichia coli", "Escherichia coli"), stringsAsFactors = FALSE)
  out <- prep_assign_organism_group(dat)
  expect_true(all(out$org_group == out$org_group[1]))
})

test_that("prep_assign_organism_group: missingness", {
  dat <- data.frame(organism_normalized = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_assign_organism_group(dat)
  expect_true(all(out$org_group == "Other" | is.na(out$org_group)))
})

test_that("prep_assign_organism_group: error handling", {
  expect_error(prep_assign_organism_group(data.frame(x = 1)), "not found")
})

# prep_flag_organism_unmatched -------------------------------------------------

test_that("prep_flag_organism_unmatched: valid behavior", {
  dat <- data.frame(organism_normalized = c("Escherichia coli", "nonexistent organism"), stringsAsFactors = FALSE)
  out <- prep_flag_organism_unmatched(dat)
  expect_true("is_organism_unmatched" %in% names(out))
  expect_false(out$is_organism_unmatched[1])
  expect_true(out$is_organism_unmatched[2])
})

test_that("prep_flag_organism_unmatched: edge cases", {
  dat <- data.frame(organism_normalized = c("escherichia coli", "Escherichia coli"), stringsAsFactors = FALSE)
  out <- prep_flag_organism_unmatched(dat)
  expect_false(any(out$is_organism_unmatched))
})

test_that("prep_flag_organism_unmatched: missingness", {
  dat <- data.frame(organism_normalized = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_flag_organism_unmatched(dat)
  expect_false(any(out$is_organism_unmatched, na.rm = TRUE))
})

test_that("prep_flag_organism_unmatched: error handling", {
  expect_error(prep_flag_organism_unmatched(data.frame(x = 1)), "not found")
})

# prep_standardize_antibiotics -------------------------------------------------

test_that("prep_standardize_antibiotics: valid behavior", {
  dat <- cat2_standardization_organism_abx()[c(1, 5), "antibiotic_name", drop = FALSE]
  who <- data.frame(Antibiotic = c("Amikacin", "Ampicillin"), Class = c("Aminoglycosides", "Aminopenicillins"), Category = c("Access", "Access"), stringsAsFactors = FALSE)
  out <- prep_standardize_antibiotics(dat, who_table = who)
  expect_true(all(c("antibiotic_normalized", "antibiotic_class", "aware_category") %in% names(out)))
})

test_that("prep_standardize_antibiotics: edge cases", {
  dat <- data.frame(antibiotic_name = c("  amikacin ", "UnknownDrugZZ"), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = c("Amikacin"), Class = c("Aminoglycosides"), Category = c("Access"), stringsAsFactors = FALSE)
  out <- prep_standardize_antibiotics(dat, who_table = who)
  expect_equal(out$antibiotic_normalized[1], "amikacin")
  expect_equal(out$antibiotic_normalized[2], "unknowndrugzz")
})

test_that("prep_standardize_antibiotics: missingness", {
  dat <- data.frame(antibiotic_name = c(NA, ""), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = c("Amikacin"), Class = c("Aminoglycosides"), Category = c("Access"), stringsAsFactors = FALSE)
  out <- prep_standardize_antibiotics(dat, who_table = who)
  expect_true(all(is.na(out$antibiotic_normalized)))
})

test_that("prep_standardize_antibiotics: error handling", {
  expect_error(prep_standardize_antibiotics(data.frame(x = 1), antibiotic_col = "antibiotic_name"), "not found")
})

# prep_classify_antibiotic_class -----------------------------------------------

test_that("prep_classify_antibiotic_class: valid behavior", {
  dat <- data.frame(antibiotic_normalized = c("Amikacin", "Ampicillin"), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = c("Amikacin", "Ampicillin"), Class = c("Aminoglycosides", "Aminopenicillins"), stringsAsFactors = FALSE)
  out <- prep_classify_antibiotic_class(dat, who_table = who)
  expect_equal(out$antibiotic_class, c("Aminoglycosides", "Aminopenicillins"))
})

test_that("prep_classify_antibiotic_class: edge cases", {
  dat <- data.frame(antibiotic_normalized = c("Unknown", "Amikacin"), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = c("Amikacin"), Class = c("Aminoglycosides"), stringsAsFactors = FALSE)
  out <- prep_classify_antibiotic_class(dat, who_table = who)
  expect_true(is.na(out$antibiotic_class[1]))
  expect_equal(out$antibiotic_class[2], "Aminoglycosides")
})

test_that("prep_classify_antibiotic_class: missingness", {
  dat <- data.frame(antibiotic_normalized = c(NA, ""), stringsAsFactors = FALSE)
  who <- data.frame(Antibiotic = c("Amikacin"), Class = c("Aminoglycosides"), stringsAsFactors = FALSE)
  out <- prep_classify_antibiotic_class(dat, who_table = who)
  expect_true(all(is.na(out$antibiotic_class)))
})

test_that("prep_classify_antibiotic_class: error handling", {
  dat <- data.frame(antibiotic_normalized = "Amikacin", stringsAsFactors = FALSE)
  expect_warning(out <- prep_classify_antibiotic_class(dat, who_table = NULL), "not provided")
  expect_equal(out$class_source, "needs_who_table")
})
