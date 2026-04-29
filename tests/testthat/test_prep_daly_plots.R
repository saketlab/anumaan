# Tests for prep_daly_plots.R

# -- Fixtures ------------------------------------------------------------------

make_hospital_counts <- function() {
  data.frame(
    center_name  = c("Centre A", "Centre B", "Centre C"),
    deaths_h     = c(50L, 42L, 58L),
    discharged_h = c(200L, 175L, 225L),
    cases_h      = c(250L, 217L, 283L),
    stringsAsFactors = FALSE
  )
}

make_hospital_daly <- function() {
  compute_hospital_daly(
    hospital_counts  = make_hospital_counts(),
    total_deaths     = 150L,
    total_discharged = 600L,
    yll_base         = 500,
    yll_associated   = 350,
    yll_attributable = 95,
    yld_base         = 80,
    yld_associated   = 55,
    yld_attributable = 15
  )
}

make_organism_yll <- function() {
  data.frame(
    organism_name      = c("e.coli", "k.pneumoniae", "s.aureus",
                           "a.baumannii", "p.aeruginosa"),
    YLL_associated_k   = c(52.3, 38.1, 29.5, 22.8, 18.2),
    YLL_attributable_k = c(14.5, 10.2,  7.8,  6.1,  4.9),
    stringsAsFactors   = FALSE
  )
}

make_yll_heatmap <- function() {
  data.frame(
    pathogen  = rep(c("e.coli", "k.pneumoniae", "s.aureus"), each = 3),
    profile   = rep(c("3GC_R", "FQ_R", "CARB_R"), times = 3),
    yll_value = c(12.5, 8.3, 5.1, 9.8, 7.2, 3.4, 6.7, 4.5, 2.1),
    stringsAsFactors = FALSE
  )
}

make_yld <- function() {
  data.frame(
    pathogen         = c("e.coli", "k.pneumoniae", "s.aureus", "a.baumannii"),
    YLD_associated   = c(8.4, 6.2, 4.8, 3.1),
    YLD_attributable = c(2.1, 1.6, 1.2, 0.8),
    stringsAsFactors = FALSE
  )
}

# -- compute_hospital_daly() ---------------------------------------------------

test_that("compute_hospital_daly returns data frame with per-1000 columns", {
  result <- make_hospital_daly()
  expect_s3_class(result, "data.frame")
  expect_true(all(c("YLL_assoc_per_1000", "YLL_attr_per_1000",
                    "YLD_assoc_per_1000", "YLD_attr_per_1000",
                    "DALY_assoc_per_1000", "DALY_attr_per_1000") %in% names(result)))
  expect_equal(nrow(result), 3L)
})

# -- plot_burden_by_hospital() -------------------------------------------------

test_that("plot_burden_by_hospital returns ggplot and errors on missing column", {
  expect_s3_class(plot_burden_by_hospital(make_hospital_daly(), metric = "YLL"), "ggplot")
  bad <- make_hospital_daly(); bad$YLL_assoc_per_1000 <- NULL
  expect_error(plot_burden_by_hospital(bad, metric = "YLL"), "Column\\(s\\) not found")
})

# -- plot_burden_by_organism() -------------------------------------------------

test_that("plot_burden_by_organism returns ggplot and errors on missing column", {
  expect_s3_class(plot_burden_by_organism(make_organism_yll(), metric = "YLL"), "ggplot")
  bad <- make_organism_yll(); bad$YLL_associated_k <- NULL
  expect_error(plot_burden_by_organism(bad, metric = "YLL"), "Column\\(s\\) not found")
})

test_that("plot_burden_by_organism normalises when n_admissions is provided", {
  p <- plot_burden_by_organism(make_organism_yll(), metric = "YLL", n_admissions = 500)
  expect_s3_class(p, "ggplot")
})

# -- plot_yll_heatmap() --------------------------------------------------------

test_that("plot_yll_heatmap errors when value_col is missing and returns ggplot when supplied", {
  expect_error(plot_yll_heatmap(make_yll_heatmap(), type = "associated",
                                n_admissions = 500),
               "'value_col' is required")
  expect_s3_class(plot_yll_heatmap(make_yll_heatmap(), type = "associated",
                                   n_admissions = 500, value_col = "yll_value"),
                  "ggplot")
})

test_that("plot_yll_heatmap errors on missing value_col column in data", {
  bad <- make_yll_heatmap()
  expect_error(plot_yll_heatmap(bad, type = "associated",
                                n_admissions = 500, value_col = "nonexistent"),
               "Column\\(s\\) not found")
})

# -- plot_yld_heatmap() --------------------------------------------------------

test_that("plot_yld_heatmap returns ggplot and errors on missing column", {
  expect_s3_class(plot_yld_heatmap(make_yld()), "ggplot")
  bad <- make_yld(); bad$YLD_associated <- NULL
  expect_error(plot_yld_heatmap(bad), "Column\\(s\\) not found")
})

test_that("plot_yld_heatmap normalises when n_admissions is provided", {
  p <- plot_yld_heatmap(make_yld(), n_admissions = 400)
  expect_s3_class(p, "ggplot")
})
