# Category 2: prep_events.R

# prep_create_event_ids --------------------------------------------------------

test_that("prep_create_event_ids: valid behavior", {
  dat <- cat2_events_df()
  out <- prep_create_event_ids(dat)
  expect_true("event_id" %in% names(out))
  expect_true(all(!is.na(out$event_id)))
})

test_that("prep_create_event_ids: edge cases", {
  dat <- cat2_events_df()
  out <- prep_create_event_ids(dat, gap_days = 0)
  expect_true(dplyr::n_distinct(out$event_id) >= 2)
})

test_that("prep_create_event_ids: missingness", {
  dat <- cat2_events_df()
  dat$date_of_culture[1] <- as.Date(NA)
  out <- prep_create_event_ids(dat)
  expect_true("event_id" %in% names(out))
})

test_that("prep_create_event_ids: error handling", {
  expect_error(prep_create_event_ids(data.frame(x = 1)), "Missing required columns")
})

# prep_deduplicate_events ------------------------------------------------------

test_that("prep_deduplicate_events: valid behavior", {
  dat <- data.frame(event_id = c("e1", "e1", "e2"), organism_normalized = c("e coli", "e coli", "k pn"), antibiotic_normalized = c("amikacin", "amikacin", "ampicillin"), x = 1:3, stringsAsFactors = FALSE)
  out <- prep_deduplicate_events(dat, keep = "first")
  expect_equal(nrow(out), 2)
})

test_that("prep_deduplicate_events: edge cases", {
  dat <- data.frame(a = c(1, 1, 2), b = c("x", "x", "y"))
  out <- prep_deduplicate_events(dat, key_cols = c("a", "b"), keep = "none")
  expect_equal(nrow(out), 1)
})

test_that("prep_deduplicate_events: missingness", {
  dat <- data.frame(event_id = c("e1", "e1"), organism_normalized = c("e coli", "e coli"), antibiotic_normalized = c(NA, NA), stringsAsFactors = FALSE)
  out <- prep_deduplicate_events(dat, keep = "all")
  expect_equal(nrow(out), 2)
})

test_that("prep_deduplicate_events: error handling", {
  expect_error(prep_deduplicate_events(data.frame(x = 1), keep = "bad_mode"), "keep must be")
})

# prep_flag_readmission --------------------------------------------------------

test_that("prep_flag_readmission: valid behavior", {
  dat <- data.frame(patient_id = c("p1", "p1", "p1"), admission_date = as.Date(c("2020-01-01", "2020-01-10", "2020-05-01")), stringsAsFactors = FALSE)
  out <- prep_flag_readmission(dat, gap_linked_days = 30, gap_new_days = 90)
  expect_equal(out$readmission_class, c("index", "linked_readmission", "late_readmission"))
})

test_that("prep_flag_readmission: edge cases", {
  dat <- data.frame(patient_id = c("p1", "p1", "p1"), admission_date = as.Date(c("2020-01-01", "2020-02-10", "2020-03-20")), stringsAsFactors = FALSE)
  out <- prep_flag_readmission(dat, gap_linked_days = 30, gap_new_days = 90)
  expect_true(any(out$readmission_class == "new_readmission"))
})

test_that("prep_flag_readmission: missingness", {
  dat <- data.frame(patient_id = c("p1", "p1"), admission_date = as.Date(c(NA, "2020-01-10")), stringsAsFactors = FALSE)
  out <- prep_flag_readmission(dat)
  expect_true("readmission_class" %in% names(out))
})

test_that("prep_flag_readmission: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_flag_readmission(dat), "not found")
  expect_true(all(is.na(out$readmission_class)))
})

# prep_classify_readmission ----------------------------------------------------

test_that("prep_classify_readmission: valid behavior", {
  dat <- data.frame(readmission_class = c("index", "linked", "new", "late"), stringsAsFactors = FALSE)
  out <- prep_classify_readmission(dat)
  expect_equal(out$readmission_class, c("index", "linked_readmission", "new_readmission", "late_readmission"))
})

test_that("prep_classify_readmission: edge cases", {
  dat <- data.frame(readmission_class = c("FIRST", "SAME_EPISODE", "READMIT", "LATE"), stringsAsFactors = FALSE)
  out <- prep_classify_readmission(dat)
  expect_equal(out$readmission_class, c("index", "linked_readmission", "new_readmission", "late_readmission"))
})

test_that("prep_classify_readmission: missingness", {
  dat <- data.frame(readmission_class = c(NA, ""), stringsAsFactors = FALSE)
  out <- prep_classify_readmission(dat)
  expect_true(all(is.na(out$readmission_class)))
})

test_that("prep_classify_readmission: error handling", {
  dat <- data.frame(x = 1)
  expect_warning(out <- prep_classify_readmission(dat), "not found")
  expect_identical(out, dat)
})
