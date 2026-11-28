# DALY burden components: deaths, fractions, CFR, incidence, PAF, YLL

vital_reg <- function() {
  data.frame(
    cause_of_death     = c("A41.9", "A41.9", "J18.9", "J18.9", "I21.9", "N39.0"),
    is_infection_death = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    syndrome           = c("BSI", "BSI", "LRI", "LRI", NA, "UTI")
  )
}

facility <- function() {
  data.frame(
    patient_id = sprintf("P%02d", 1:10),
    syndrome = "BSI",
    specimen = "Blood",
    poly = FALSE,
    organism = rep(c("E. coli", "K. pneumoniae"), c(6, 4)),
    outcome = c(
      "Death", "Alive", "Alive", "Death", "Alive", "Alive",
      "Death", "Death", "Alive", "Alive"
    )
  )
}

test_that("population death components follow D_J, S_J, M_LJ and D_L definitions", {
  vr <- vital_reg()
  d_j <- suppressMessages(daly_calc_deaths_by_cause(vr))
  expect_equal(
    setNames(d_j$D_J, d_j$cause_of_death),
    c(A41.9 = 2L, I21.9 = 1L, J18.9 = 2L, N39.0 = 1L)
  )

  s_j <- suppressMessages(daly_calc_infection_fraction(vr))
  expect_equal(s_j$S_J[s_j$cause_of_death == "I21.9"], 0)
  expect_equal(s_j$S_J[s_j$cause_of_death == "A41.9"], 1)

  m_lj <- suppressMessages(daly_calc_syndrome_fraction(vr))
  expect_true(all(m_lj$M_LJ >= 0 & m_lj$M_LJ <= 1))

  d_l <- suppressMessages(daly_calc_deaths_by_syndrome(pop_data = vr))
  expect_equal(setNames(d_l$D_L, d_l$syndrome), c(BSI = 2, LRI = 2, UTI = 1))
  # infection deaths are conserved across syndromes
  expect_equal(sum(d_l$D_L), sum(vr$is_infection_death))
})

test_that("daly_calc_deaths_by_cause honours a pre-aggregated deaths column", {
  agg <- data.frame(cause_of_death = c("A41.9", "J18.9"), n = c(10, 4))
  d_j <- suppressMessages(daly_calc_deaths_by_cause(agg, deaths_col = "n"))
  expect_equal(setNames(d_j$D_J, d_j$cause_of_death), c(A41.9 = 10, J18.9 = 4))
})

test_that("facility counts, top pathogens, CFR and fatal pathogen fractions are correct", {
  fac <- facility()
  top <- suppressMessages(daly_get_top_pathogens(fac, "organism", n = 1))
  expect_equal(top$organism, "E. coli")

  n <- suppressMessages(daly_count_non_fatal_incident_cases(fac, "syndrome", "BSI", "patient_id",
    "outcome", discharged_value = "Alive"
  ))
  expect_equal(n$N_NF_L, 6L) # 10 cases, 4 deaths

  cfr <- suppressMessages(daly_calc_case_fatality(fac, "syndrome", "BSI", "specimen", "Blood",
    "poly", "patient_id", "organism", "outcome",
    death_value = "Death"
  ))$cfr_table
  expect_equal(setNames(cfr$CFR_LK, cfr$organism), c(`E. coli` = 2 / 6, `K. pneumoniae` = 2 / 4))

  pf <- suppressMessages(daly_calc_pathogen_fraction_fatal(fac, "syndrome", "BSI",
    polymicrobial_col = "poly", patient_col = "patient_id", pathogen_col = "organism",
    outcome_col = "outcome"
  ))$P_Lk_fatal
  expect_equal(sum(pf$P_Lk_fatal), 1)
  expect_equal(pf$N_F_LK, c(2, 2))
})

test_that("daly_calc_cr_l implements the three adjustment types", {
  adj <- utils::read.csv(system.file("extdata", "adjustment_for_CFR.csv", package = "anumaan"),
    check.names = FALSE
  )
  expect_equal(nrow(adj), 12L) # guards against unquoted commas in the reference file

  fac <- data.frame(
    patient_id = c("P1", "P1", "P2", "P3", "P4"),
    visit = c("Outpatient", "Inpatient", "Outpatient", "Outpatient", "Inpatient")
  )
  crl <- function(syn) {
    fac$syndrome <- syn
    suppressWarnings(suppressMessages(daly_calc_cr_l(fac, "syndrome", syn, "patient_id", "visit",
      adjustment_ref = adj
    )))$CR_L
  }
  expect_equal(crl("Bloodstream infections"), 1)
  # patients with >= 1 inpatient visit / all patients = 2 / 4
  expect_equal(crl("Urinary tract infections and pyelonephritis – community acquired"), 0.5)
  # OP patients who also had IP / OP patients = 1 / 3
  expect_equal(crl("Sexually transmitted infections"), 1 / 3)
  expect_warning(daly_calc_cr_l(transform(fac, syndrome = "X"), "syndrome", "X", "patient_id",
    "visit",
    adjustment_ref = adj
  ), "not found")
})

test_that("daly_calc_incidence_from_cfr divides deaths by the pathogen-weighted CFR", {
  out <- suppressMessages(daly_calc_incidence_from_cfr(
    deaths_L = 10,
    cfr_lk_tbl = data.frame(pathogen = c("E", "K"), CFR_LK = c(0.2, 0.5)),
    P_Lk_prime_tbl = data.frame(pathogen = c("E", "K"), P_Lk_prime = c(0.6, 0.4))
  ))
  expect_equal(out$CFR_L, 0.6 * 0.2 + 0.4 * 0.5)
  expect_equal(out$I_L, 10 / 0.32)
})

test_that("PAF mortality and fatal resistance prevalence match the closed-form expressions", {
  pr <- list(K = data.frame(profile = c("S", "R"), probability = c(0.6, 0.4), RR_LOS_profile = c(1, 2)))
  paf <- suppressMessages(daly_calc_paf_mortality(pr))
  expect_equal(paf$K$PAF_k_mort, round(0.4 / 1.4, 6))
  expect_equal(paf$K$denominator_mort, 1.4)

  prev <- suppressMessages(daly_calc_resistance_prevalence_fatal(pr))
  expect_equal(prev$K$E_OR_k, 1.4)
  expect_equal(prev$K$R_k, round(0.8 / 1.4, 6))

  bad <- list(K = data.frame(profile = "S", probability = 1, RR_LOS_profile = -1))
  expect_warning(out <- suppressMessages(daly_calc_paf_mortality(bad)), "denominator")
  expect_length(out, 0L)
})

test_that("YLL associated uses bundled life expectancy and YLL attributable applies PAF", {
  cohort <- data.frame(
    patient_id = sprintf("P%02d", 1:6),
    final_outcome = c("Death", "Death", "Discharged", "Death", "Discharged", "Death"),
    organism_name = rep(c("Escherichia coli", "Klebsiella pneumoniae"), 3),
    age_bin = c("40-45", "60-65", "20-25", "0-1", "70-75", "85+"),
    sex = c("Male", "Female", "Male", "Female", "Male", "Female")
  )
  yll <- suppressMessages(daly_calc_yll_associated(cohort,
    outcome_col = "final_outcome", pathogen_col = "organism_name",
    patient_col = "patient_id", age_bin_col = "age_bin", sex_col = "sex"
  ))
  # only the four deaths contribute; YLL = sex-specific life expectancy at death
  # (inst/extdata/life_expectancy_all.csv: Female 0-1 = 72.5, Female 85+ = 6.4)
  pd <- yll$patient_data
  expect_setequal(pd$patient_id, c("P01", "P02", "P04", "P06"))
  expect_equal(pd$yll_contribution[pd$patient_id == "P04"], 72.5)
  expect_equal(pd$yll_contribution[pd$patient_id == "P06"], 6.4)
  expect_equal(yll$total, sum(pd$yll_contribution))
  expect_equal(sum(yll$per_pathogen$n_patients), 4L)

  paf_mort <- list(
    "Escherichia coli" = list(PAF_k_mort = 0.2),
    "Klebsiella pneumoniae" = list(PAF_k_mort = 0.5)
  )
  attr <- suppressMessages(daly_calc_yll_attributable(yll$patient_data, paf_mort,
    pathogen_col = "organism_name", patient_col = "patient_id", age_bin_col = "age_bin"
  ))
  pd <- attr$patient_data
  expected <- pd$yll_contribution * ifelse(pd$organism_name == "Escherichia coli", 0.2, 0.5)
  expect_equal(pd$YLL_attributable_contribution, expected)
})

test_that("daly_calc_pathogen_fraction_nonfatal weights discharged patients only", {
  fac <- facility()
  fac$outcome[fac$outcome == "Alive"] <- "Discharged"
  out <- suppressMessages(daly_calc_pathogen_fraction_nonfatal(fac, "syndrome", "BSI",
    polymicrobial_col = "poly", patient_col = "patient_id", pathogen_col = "organism",
    outcome_col = "outcome"
  ))$P_Lk_prime
  # 6 discharged: 4 E. coli, 2 K. pneumoniae
  expect_equal(setNames(out$P_Lk_prime, out$organism), c(`E. coli` = 4 / 6, `K. pneumoniae` = 2 / 6))
  expect_equal(sum(out$P_Lk_prime), 1)
})

test_that("daly_count_non_fatal_incident_cases counts distinct discharged admissions per hospital x syndrome", {
  adm <- data.frame(
    readmission_id = c("A1", "A1", "A1", "A2", "A3", "A4", "A5", "A6", "A7"),
    hospital = c("H1", "H1", "H1", "H1", "H1", "H2", "H2", "H2", "H2"),
    syndrome = c("BSI", "BSI", "BSI", "BSI", "UTI", "bsi", "", NA, "UTI"),
    outcome = c("Discharged", "Discharged", "Discharged", "Death", "Discharged",
                "Discharged", "Discharged", "Discharged", "LAMA")
  )
  expect_message(
    out <- daly_count_non_fatal_incident_cases(adm, "syndrome", NULL, "readmission_id", "outcome",
      facility_col = "hospital"
    ),
    "2 discharged admission\\(s\\) have no syndrome"
  )
  out <- out[order(out$hospital, out$syndrome), ]
  # A1 (3 rows) counts once; A2 death, A7 LAMA excluded; A5/A6 have no syndrome
  expect_equal(out$hospital, c("H1", "H1", "H2"))
  expect_equal(tolower(out$syndrome), c("bsi", "uti", "bsi"))
  expect_equal(out$N_NF_L, c(1L, 1L, 1L))

  bsi <- suppressMessages(daly_count_non_fatal_incident_cases(adm, "syndrome", "BSI", "readmission_id",
    "outcome", facility_col = "hospital", facility_name = "H1"
  ))
  expect_equal(bsi$N_NF_L, 1L)
  pooled <- suppressMessages(daly_count_non_fatal_incident_cases(adm, "syndrome", "bsi", "readmission_id", "outcome"))
  expect_equal(pooled$N_NF_L, 2L) # A1 and A4, case-insensitive syndrome match

  two <- rbind(adm, data.frame(readmission_id = "A1", hospital = "H1", syndrome = "UTI", outcome = "Discharged"))
  expect_warning(
    suppressMessages(daly_count_non_fatal_incident_cases(two, "syndrome", NULL, "readmission_id", "outcome")),
    "more than one syndrome"
  )
  expect_error(daly_count_non_fatal_incident_cases(adm, "syndrome", NULL, "readmission_id", "missing_col"), "Missing column")
})
