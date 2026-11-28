# Organism, drug and syndrome names are compared through a lower-case key, so
# user input in any case matches reference tables and user-supplied maps.

test_that(".name_key lower-cases, trims and collapses whitespace", {
  expect_equal(
    .name_key(c(" Klebsiella  Pneumoniae ", "E. COLI", NA)),
    c("klebsiella pneumoniae", "e. coli", NA)
  )
})

test_that("organism group and RR pathogen resolve regardless of input case", {
  d <- data.frame(
    organism_normalized = c("Escherichia coli", "KLEBSIELLA PNEUMONIAE", "escherichia coli"),
    antibiotic_class = "Carbapenems"
  )
  g <- suppressMessages(prep_assign_organism_group(d))
  expect_equal(g$org_group, rep("Enterobacterales", 3))

  rr <- suppressMessages(daly_add_rr_mappings(d))
  expect_equal(rr$rr_pathogen, c("escherichia coli", "klebsiella pneumoniae", "escherichia coli"))
  expect_true(all(rr$rr_pathogen %in% suppressMessages(daly_load_rr_reference())$pathogen))
})

test_that("run_preprocess output carries org_group and rr_pathogen", {
  raw <- data.frame(
    PID = c("p1", "p2"), Organism = c("E. coli", "Klebsiella pneumoniae"),
    Antibiotic = c("Meropenem", "Ceftriaxone"), Result = c("R", "S"),
    CultureDate = "2021-01-03", Sample_type1_name = "Blood"
  )
  r <- suppressWarnings(suppressMessages(run_preprocess(raw, verbose = FALSE, generate_report = FALSE)))
  expect_equal(r$data$org_group, c("Enterobacterales", "Enterobacterales"))
  expect_false(anyNA(r$data$rr_pathogen))
})

test_that("panel_map names match pathogen values case-insensitively", {
  skip_if_not_installed("quadprog")
  m <- data.frame(
    pathogen = "klebsiella pneumoniae", antibiotic_class = c("A", "B"),
    n_tested = 100L, n_resistant = c(30L, 50L)
  )
  m$marginal_resistance <- m$n_resistant / m$n_tested
  out <- suppressMessages(estimate_profiles_convex(m,
    panel_map = list("Klebsiella Pneumoniae" = c("A", "B")), solver = "quadprog"
  ))
  expect_equal(nrow(out), 4L)
  expect_equal(unique(out$pathogen), "klebsiella pneumoniae")
})

test_that("syndrome and pathogen filters are case-agnostic", {
  fac <- data.frame(
    patient_id = sprintf("P%02d", 1:4), syndrome = "Bloodstream infections",
    organism = c("E. coli", "E. coli", "K. pneumoniae", "K. pneumoniae"),
    specimen = "Blood", poly = FALSE,
    outcome = c("Death", "Alive", "Death", "Death")
  )
  n <- suppressMessages(daly_count_non_fatal_incident_cases(fac, "syndrome", "BLOODSTREAM INFECTIONS",
    "patient_id", "outcome", discharged_value = "alive"
  ))
  expect_equal(n$N_NF_L, 1L)

  cfr <- suppressMessages(daly_calc_case_fatality(fac, "syndrome", "bloodstream infections",
    "specimen", "Blood", "poly", "patient_id", "organism", "outcome",
    death_value = "Death", pathogen_name = "e. coli"
  ))$cfr_table
  expect_equal(cfr$CFR_LK, 0.5)

  adj <- utils::read.csv(system.file("extdata", "adjustment_for_CFR.csv", package = "anumaan"),
    check.names = FALSE
  )
  crl <- suppressMessages(daly_calc_cr_l(transform(fac, visit = "Inpatient"), "syndrome",
    "bloodstream infections", "patient_id", "visit",
    adjustment_ref = adj
  ))
  expect_equal(crl$adjustment_type, "None") # found in the reference despite case
})

test_that("case-insensitive panel matching keeps the data's own pathogen labels", {
  skip_if_not_installed("quadprog")
  m <- data.frame(
    pathogen = "Klebsiella pneumoniae", antibiotic_class = c("A", "B"),
    n_tested = 100L, n_resistant = c(30L, 50L)
  )
  m$marginal_resistance <- m$n_resistant / m$n_tested
  out <- suppressMessages(estimate_profiles_convex(m,
    panel_map = list("KLEBSIELLA PNEUMONIAE" = c("A", "B")), solver = "quadprog"
  ))
  expect_equal(unique(out$pathogen), "Klebsiella pneumoniae")
})
