# Tests for prep_eda_plots.R

make_amr_data <- function() {
  data.frame(
    PatientInformation_id = c("P1", "P1", "P2", "P2", "P3", "P4", "P5", "P6"),
    organism_name         = c("e.coli", "e.coli", "k.pneumoniae", "k.pneumoniae",
                              "s.aureus", "e.coli", "k.pneumoniae", "s.aureus"),
    antibiotic_name       = c("amoxicillin", "ciprofloxacin", "amoxicillin",
                              "ciprofloxacin", "amoxicillin", "ciprofloxacin",
                              "amoxicillin", "ciprofloxacin"),
    antibiotic_value      = c("R", "S", "R", "R", "S", "R", "S", "R"),
    center_name           = c("Centre A", "Centre A", "Centre A", "Centre A",
                              "Centre B", "Centre B", "Centre B", "Centre B"),
    final_outcome         = c("Discharged", "Discharged", "Death", "Death",
                              "LAMA", "Discharged", "Referred", "Death"),
    sample_type           = c("Blood", "Blood", "Urine", "Urine",
                              "Blood", "Urine", "Blood", "Urine"),
    Age_bin               = factor(
      c("<5", "18-44", "45-64", "65+", "<5", "18-44", "45-64", "65+"),
      levels = c("<5", "18-44", "45-64", "65+")
    ),
    is_polymicrobial      = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE),
    admission_date        = as.Date(c("2023-01-01", "2023-01-01", "2023-01-05",
                                      "2023-01-05", "2023-02-10", "2023-02-10",
                                      "2023-03-01", "2023-03-01")),
    culture_date          = as.Date(c("2023-01-03", "2023-01-03", "2023-01-04",
                                      "2023-01-04", "2023-02-15", "2023-02-12",
                                      "2023-03-02", "2023-03-04")),
    final_outcome_date    = as.Date(c("2023-01-10", "2023-01-10", "2023-01-20",
                                      "2023-01-20", "2023-02-20", "2023-02-18",
                                      "2023-03-08", "2023-03-09")),
    location              = c("ICU", "Ward", "ICU", "Ward",
                               "ICU", "Ward", "Other", "ICU"),
    patient_id            = c("P1", "P1", "P2", "P2", "P3", "P4", "P5", "P6"),
    age_years             = c(3, 30, 55, 70, 2, 25, 50, 68),
    stringsAsFactors      = FALSE
  )
}

test_that("plot_top_organisms returns ggplot and errors on missing column", {
  expect_s3_class(plot_top_organisms(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$organism_name <- NULL
  expect_error(plot_top_organisms(bad, mode = "overall"), "Column\\(s\\) not found")
})

test_that("plot_abx_susceptibility returns ggplot and errors on missing column", {
  expect_s3_class(plot_abx_susceptibility(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$antibiotic_value <- NULL
  expect_error(plot_abx_susceptibility(bad, mode = "overall"), "Column\\(s\\) not found")
})

test_that("plot_abx_heatmap returns ggplot and errors when center missing in single mode", {
  expect_s3_class(plot_abx_heatmap(make_amr_data(), mode = "all"), "ggplot")
  expect_error(plot_abx_heatmap(make_amr_data(), mode = "single"), "'center' must be provided")
})

test_that("plot_outcome_distribution returns ggplot and errors on missing column", {
  expect_s3_class(plot_outcome_distribution(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$final_outcome <- NULL
  expect_error(plot_outcome_distribution(bad), "Column\\(s\\) not found")
})

test_that("plot_outcome_by_organism returns ggplot and errors on missing column", {
  expect_s3_class(plot_outcome_by_organism(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$antibiotic_value <- NULL
  expect_error(plot_outcome_by_organism(bad, mode = "overall"), "Column\\(s\\) not found")
})

test_that("plot_death_discharged returns ggplot and errors on missing column", {
  expect_s3_class(plot_death_discharged(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$final_outcome <- NULL
  expect_error(plot_death_discharged(bad), "Column\\(s\\) not found")
})

test_that("plot_resistance_by_sample returns ggplot and errors on missing column", {
  expect_s3_class(plot_resistance_by_sample(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$sample_type <- NULL
  expect_error(plot_resistance_by_sample(bad, mode = "overall"), "Column\\(s\\) not found")
})

test_that("plot_outcome_by_agebin returns ggplot and errors on missing column", {
  expect_s3_class(plot_outcome_by_agebin(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$Age_bin <- NULL
  expect_error(plot_outcome_by_agebin(bad), "Column\\(s\\) not found")
})

test_that("plot_mono_poly_by_facility returns ggplot and errors on missing column", {
  expect_s3_class(plot_mono_poly_by_facility(make_amr_data()), "ggplot")
  bad <- make_amr_data(); bad$is_polymicrobial <- NULL
  expect_error(plot_mono_poly_by_facility(bad), "Column\\(s\\) not found")
})

test_that("plot_hai_cai_by_facility returns ggplot and errors on missing column", {
  expect_s3_class(plot_hai_cai_by_facility(make_amr_data()), "ggplot")
  bad <- make_amr_data(); bad$admission_date <- NULL
  expect_error(plot_hai_cai_by_facility(bad), "Column\\(s\\) not found")
})

test_that("plot_location_by_facility returns ggplot and errors on missing column", {
  expect_s3_class(plot_location_by_facility(make_amr_data()), "ggplot")
  bad <- make_amr_data(); bad$location <- NULL
  expect_error(plot_location_by_facility(bad), "Column\\(s\\) not found")
})

test_that("plot_los_ridge returns ggplot and errors when center missing in single mode", {
  expect_s3_class(plot_los_ridge(make_amr_data(), mode = "all"), "ggplot")
  expect_error(plot_los_ridge(make_amr_data(), mode = "single"), "'center' must be provided")
})

test_that("plot_age_ridge returns ggplot and errors on missing column", {
  expect_s3_class(plot_age_ridge(make_amr_data(), mode = "all"), "ggplot")
  bad <- make_amr_data(); bad$age_years <- NULL
  expect_error(plot_age_ridge(bad), "required columns are missing")
})

test_that("plot_los_by_agebin returns ggplot and errors on missing column", {
  expect_s3_class(plot_los_by_agebin(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$admission_date <- NULL
  expect_error(plot_los_by_agebin(bad), "required columns are missing")
})

test_that("plot_outcome_by_year returns ggplot and errors on missing column", {
  expect_s3_class(plot_outcome_by_year(make_amr_data(), mode = "overall"), "ggplot")
  bad <- make_amr_data(); bad$final_outcome_date <- NULL
  expect_error(plot_outcome_by_year(bad), "Column\\(s\\) not found")
})
