# Category 2: prep_ast.R

# prep_clean_ast_values --------------------------------------------------------

test_that("prep_clean_ast_values: valid behavior", {
  dat <- data.frame(antibiotic_value = c("R", "S", "I"), stringsAsFactors = FALSE)
  out <- prep_clean_ast_values(dat)
  expect_equal(out$antibiotic_value, c("R", "S", "I"))
})

test_that("prep_clean_ast_values: edge cases", {
  dat <- data.frame(antibiotic_value = c("Resistant", "Susceptible", "Intermediate", "R   e coli"), stringsAsFactors = FALSE)
  out <- prep_clean_ast_values(dat)
  expect_equal(out$antibiotic_value, c("R", "S", "I", "R"))
})

test_that("prep_clean_ast_values: missingness", {
  dat <- data.frame(antibiotic_value = c(NA, "", "-"), stringsAsFactors = FALSE)
  out <- prep_clean_ast_values(dat)
  expect_true(all(is.na(out$antibiotic_value)))
})

test_that("prep_clean_ast_values: error handling", {
  expect_error(prep_clean_ast_values(data.frame(x = 1), value_col = "antibiotic_value"), "not found")
})

# prep_recode_intermediate_ast -------------------------------------------------

test_that("prep_recode_intermediate_ast: valid behavior", {
  dat <- data.frame(antibiotic_name = c("colistin", "amikacin"), antibiotic_value = c("I", "I"), stringsAsFactors = FALSE)
  out <- prep_recode_intermediate_ast(dat)
  expect_equal(out$antibiotic_value, c("S", "R"))
})

test_that("prep_recode_intermediate_ast: edge cases", {
  dat <- data.frame(antibiotic_name = c("colistin", "amikacin"), antibiotic_value = c("I", "I"), stringsAsFactors = FALSE)
  out <- prep_recode_intermediate_ast(dat, colistin_to_s = FALSE, others_to_r = TRUE)
  expect_equal(out$antibiotic_value, c("I", "R"))
})

test_that("prep_recode_intermediate_ast: missingness", {
  dat <- data.frame(antibiotic_name = c("colistin", NA), antibiotic_value = c(NA, "I"), stringsAsFactors = FALSE)
  out <- prep_recode_intermediate_ast(dat)
  expect_true(is.na(out$antibiotic_value[1]))
})

test_that("prep_recode_intermediate_ast: error handling", {
  expect_error(prep_recode_intermediate_ast(data.frame(x = 1), antibiotic_col = "antibiotic_name"), "not found")
})

# prep_harmonize_ast -----------------------------------------------------------

test_that("prep_harmonize_ast: valid behavior", {
  dat <- data.frame(antibiotic_name_std = c("colistin", "amikacin"), ast_value_raw = c("I", "Resistant"), stringsAsFactors = FALSE)
  out <- prep_harmonize_ast(dat)
  expect_true(all(c("ast_value_harmonized", "intermediate_recoded") %in% names(out)))
  expect_equal(out$ast_value_harmonized, c("S", "R"))
})

test_that("prep_harmonize_ast: edge cases", {
  dat <- data.frame(ast_value_raw = c("S", "R"), stringsAsFactors = FALSE)
  out <- prep_harmonize_ast(dat)
  expect_equal(out$ast_value_harmonized, c("S", "R"))
})

test_that("prep_harmonize_ast: missingness", {
  dat <- data.frame(antibiotic_name_std = c("colistin", "amikacin"), ast_value_raw = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_harmonize_ast(dat)
  expect_true(all(is.na(out$ast_value_harmonized)))
})

test_that("prep_harmonize_ast: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_harmonize_ast(dat, ast_col = "ast_value_raw"), "not found")
  expect_true(all(is.na(out$ast_value_harmonized)))
})

# prep_check_organism_ast_consistency ------------------------------------------

test_that("prep_check_organism_ast_consistency: valid behavior", {
  dat <- data.frame(organism_normalized = "escherichia coli", antibiotic_normalized = "amikacin", organism_group = "Enterobacterales", stringsAsFactors = FALSE)
  out <- prep_check_organism_ast_consistency(dat)
  expect_true("is_ast_inconsistent" %in% names(out))
})

test_that("prep_check_organism_ast_consistency: edge cases", {
  dat <- data.frame(organism_normalized = c("escherichia coli", "escherichia coli"), antibiotic_normalized = c("amikacin", "ampicillin"), organism_group = c("Enterobacterales", "Enterobacterales"), stringsAsFactors = FALSE)
  out <- prep_check_organism_ast_consistency(dat)
  expect_equal(nrow(out), 2)
})

test_that("prep_check_organism_ast_consistency: missingness", {
  dat <- data.frame(organism_normalized = c(NA, "escherichia coli"), antibiotic_normalized = c("amikacin", NA), stringsAsFactors = FALSE)
  out <- prep_check_organism_ast_consistency(dat)
  expect_true("is_ast_inconsistent" %in% names(out))
})

test_that("prep_check_organism_ast_consistency: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_check_organism_ast_consistency(dat), "not found")
  expect_true(all(is.na(out$is_ast_inconsistent)))
})

# prep_flag_invalid_ast --------------------------------------------------------

test_that("prep_flag_invalid_ast: valid behavior", {
  dat <- data.frame(ast_value_harmonized = c("S", "I", "R"), stringsAsFactors = FALSE)
  out <- prep_flag_invalid_ast(dat)
  expect_false(any(out$is_ast_invalid))
})

test_that("prep_flag_invalid_ast: edge cases", {
  dat <- data.frame(ast_value_harmonized = c("R", "X", "S"), stringsAsFactors = FALSE)
  out <- prep_flag_invalid_ast(dat)
  expect_true(out$is_ast_invalid[2])
})

test_that("prep_flag_invalid_ast: missingness", {
  dat <- data.frame(ast_value_harmonized = c(NA, "", "R"), stringsAsFactors = FALSE)
  out <- prep_flag_invalid_ast(dat)
  expect_true(out$is_ast_invalid[2])
})

test_that("prep_flag_invalid_ast: error handling", {
  expect_warning(out <- prep_flag_invalid_ast(data.frame(x = 1), col = "ast_value_harmonized"), "not found")
  expect_true(all(is.na(out$is_ast_invalid)))
})

# prep_deduplicate_ast ---------------------------------------------------------
# Fixture: cat2_ast_long_df() = cat2_ast_duplicates_same_day() — 5 rows:
#   rows 1-2: pt_020 | e.coli | amikacin | 2025-01-15 | S  (exact duplicate)
#   row  3:   pt_020 | e.coli | amikacin | 2025-01-15 | R  (conflict with row 1 after exact-dedup)
#   row  4:   pt_020 | e.coli | cefotaxime | 2025-01-15 | I
#   row  5:   pt_021 | k.pneumo | ampicillin | 2025-01-16 | S

test_that("prep_deduplicate_ast: detect mode removes exact dups then flags conflicts", {
  dat <- cat2_ast_long_df()
  out <- prep_deduplicate_ast(dat, mode = "detect")
  expect_true("is_ast_duplicate" %in% names(out))
  # row 2 is an exact dup of row 1 — removed in step 1 → 4 rows remain
  expect_equal(nrow(out), 4L)
  # amikacin S vs R is still a conflict after exact-dedup
  expect_true(any(out$is_ast_duplicate))
})

test_that("prep_deduplicate_ast: remove mode resolves both exact dups and conflicts", {
  dat <- cat2_ast_long_df()
  out <- prep_deduplicate_ast(dat, mode = "remove", strategy = "resistant_wins")
  expect_false("is_ast_duplicate" %in% names(out))
  # amikacin conflict: R wins; cefotaxime and ampicillin have no conflict → 3 rows
  expect_equal(nrow(out), 3L)
})

test_that("prep_deduplicate_ast: exact duplicates only — no conflicts flagged", {
  dat <- cat2_ast_exact_duplicates()
  out <- prep_deduplicate_ast(dat, mode = "detect")
  # rows 1-2 are exact dups → reduced to 1 row before conflict detection
  expect_equal(nrow(out), 2L)
  expect_false(any(out$is_ast_duplicate))
})

test_that("prep_deduplicate_ast: remove mode on exact-dup-only data is a no-op for conflicts", {
  dat <- cat2_ast_exact_duplicates()
  out <- prep_deduplicate_ast(dat, mode = "remove", strategy = "resistant_wins")
  expect_equal(nrow(out), 2L)
})

test_that("prep_deduplicate_ast: missingness — NA value breaks exact-dup grouping", {
  dat <- cat2_ast_long_df()
  dat$ast_value_harmonized[1] <- NA
  # rows 1 (NA) and 2 (S) are no longer exact dups → all 5 rows survive exact-dedup
  # amikacin group {NA, S, R}: n_distinct(na.rm=TRUE) = 2 (S and R) → conflict flagged
  out <- prep_deduplicate_ast(dat, mode = "detect")
  expect_true("is_ast_duplicate" %in% names(out))
  expect_equal(nrow(out), 5L)
})

test_that("prep_deduplicate_ast: error handling — warns and returns unchanged on missing columns", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_deduplicate_ast(dat, mode = "detect"), "not found")
  expect_identical(out, dat)
})

# prep_pivot_ast_wide_to_long --------------------------------------------------

test_that("prep_pivot_ast_wide_to_long: valid behavior", {
  dat <- data.frame(patient_id = c("p1", "p2"), AMK = c("S", "R"), CIP = c("R", "S"), stringsAsFactors = FALSE)
  out <- prep_pivot_ast_wide_to_long(dat, antibiotic_cols = c("AMK", "CIP"), remove_missing = TRUE)
  expect_true(all(c("antibiotic_name", "antibiotic_value") %in% names(out)))
  expect_equal(nrow(out), 4)
})

test_that("prep_pivot_ast_wide_to_long: edge cases", {
  dat <- data.frame(patient_id = c("p1", "p2"), AMK = c("S", NA), stringsAsFactors = FALSE)
  out <- prep_pivot_ast_wide_to_long(dat, antibiotic_cols = c("AMK"), remove_missing = FALSE, create_event_id = TRUE)
  expect_true("event_id" %in% names(out))
})

test_that("prep_pivot_ast_wide_to_long: missingness", {
  dat <- data.frame(patient_id = c("p1", "p2"), AMK = c("", "-"), stringsAsFactors = FALSE)
  out <- prep_pivot_ast_wide_to_long(dat, antibiotic_cols = c("AMK"), remove_missing = TRUE)
  expect_equal(nrow(out), 0)
})

test_that("prep_pivot_ast_wide_to_long: error handling", {
  dat <- data.frame(patient_id = "p1", AMK = "S")
  expect_error(prep_pivot_ast_wide_to_long(dat, antibiotic_cols = NULL, pattern = NULL), "must be provided")
})

# prep_create_wide_ast_matrix --------------------------------------------------

test_that("prep_create_wide_ast_matrix: valid behavior", {
  dat <- data.frame(event_id = c("e1", "e1", "e2"), antibiotic_normalized = c("amikacin", "ciprofloxacin", "amikacin"), antibiotic_value = c("S", "R", "R"), patient_id = c("p1", "p1", "p2"), organism_normalized = c("e coli", "e coli", "k pn"), date_of_culture = as.Date(c("2020-01-01", "2020-01-01", "2020-01-02")), stringsAsFactors = FALSE)
  out <- prep_create_wide_ast_matrix(dat)
  expect_true("event_id" %in% names(out))
  expect_true(any(grepl("abx_", names(out))))
})

test_that("prep_create_wide_ast_matrix: edge cases", {
  dat <- data.frame(event_id = c("e1", "e1"), antibiotic_normalized = c("amikacin", "amikacin"), antibiotic_value = c("S", "R"), patient_id = c("p1", "p1"), organism_normalized = c("e coli", "e coli"), date_of_culture = as.Date(c("2020-01-01", "2020-01-01")), stringsAsFactors = FALSE)
  out <- prep_create_wide_ast_matrix(dat)
  expect_equal(nrow(out), 1)
})

test_that("prep_create_wide_ast_matrix: missingness", {
  dat <- data.frame(event_id = c("e1", "e2"), antibiotic_normalized = c("amikacin", "ciprofloxacin"), antibiotic_value = c(NA, "R"), patient_id = c("p1", "p2"), organism_normalized = c("e coli", "k pn"), date_of_culture = as.Date(c("2020-01-01", "2020-01-02")), stringsAsFactors = FALSE)
  out <- prep_create_wide_ast_matrix(dat)
  expect_equal(nrow(out), 2)
})

test_that("prep_create_wide_ast_matrix: error handling", {
  expect_error(prep_create_wide_ast_matrix(data.frame(x = 1)), "not found")
})
