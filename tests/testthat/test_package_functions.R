# Core package behavior outside Category 1 (kept minimal)

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

test_that("prep_derive_los_from_dates computes length of stay", {
  test_data <- data.frame(
    date_of_admission = as.Date(c("2020-01-01", "2020-02-01")),
    date_of_final_outcome = as.Date(c("2020-01-10", "2020-02-05"))
  )
  result <- prep_derive_los_from_dates(
    test_data,
    admission_col = "date_of_admission",
    outcome_col = "date_of_final_outcome"
  )
  expect_equal(result$los_days, c(9, 4))
})
