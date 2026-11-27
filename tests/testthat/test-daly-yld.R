# Tests for R/daly_yld.R and the LOS-propagation extension to
# daly_assign_rr_to_profiles() in R/daly_rr_and_los.R.
#
# Covers:
#   1. daly_assign_rr_to_profiles() backward compatibility (no los_r_col/
#      los_s_col -> output identical to before this task's extension).
#   2. daly_assign_rr_to_profiles() LOS propagation: LOSR_years/LOSS_years
#      assigned to the profile's dominant (max-RR) class, NA for the
#      all-susceptible profile, and days -> years conversion.
#   3. daly_calc_yld_associated()/daly_calc_yld_attributable() against a
#      small, hand-computable synthetic example, verifying the direct
#      equations exactly:
#        YLD_assoc   = I_L * P'_Lk * DW_L * sum_resistant[R' * LOSR]
#        YLD_attrib  = I_L * P'_Lk * DW_L * sum_resistant[R' * (LOSR - LOSS)]
#   4. Neither primary function requires a pre-computed baseline YLD or an
#      associated-fraction/PAF object.
#   5. daly_calc_paf_los() still works standalone (unchanged, decoupled).
#   6. daly_calc_yld_baseline() and daly_calc_fraction_associated_yld() no
#      longer exist (removed as pure intermediates of the old pathway).

.synthetic_profiles_output <- function() {
  # 2-class pathogen: classes "A", "B". Profiles (alphabetical class order,
  # matching compute_resistance_profiles()'s convention): SS, SB (B-resistant
  # only), SA (A-resistant only), AB (both resistant).
  profiles <- data.frame(
    profile = c("SS", "SB", "SA", "AB"),
    A = c(0L, 0L, 1L, 1L),
    B = c(0L, 1L, 0L, 1L),
    probability = c(0.50, 0.20, 0.20, 0.10),
    stringsAsFactors = FALSE
  )
  list(Test_Pathogen = list(profiles = profiles, classes = c("A", "B")))
}

.synthetic_rr_table_days <- function() {
  # Class A: larger RR (dominant whenever resistant to both A and B).
  data.frame(
    pathogen = c("Test_Pathogen", "Test_Pathogen"),
    antibiotic_class = c("A", "B"),
    RR_LOS = c(2.0, 1.5),
    median_LOS_R = c(20, 15), # days
    median_LOS_S = c(10, 10), # days
    stringsAsFactors = FALSE
  )
}

test_that("daly_assign_rr_to_profiles() is unchanged when los_r_col/los_s_col are not supplied", {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()

  out <- daly_assign_rr_to_profiles(profiles_output, rr_table)

  expect_false("LOSR_years" %in% names(out$Test_Pathogen))
  expect_false("LOSS_years" %in% names(out$Test_Pathogen))
  expect_equal(out$Test_Pathogen$dominant_class, c("all_susceptible", "B", "A", "A"))
  expect_equal(out$Test_Pathogen$RR_LOS_profile, c(1.0, 1.5, 2.0, 2.0))
})

test_that("daly_assign_rr_to_profiles() propagates absolute LOSR/LOSS to the dominant class, NA for all-susceptible", {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()

  out <- daly_assign_rr_to_profiles(
    profiles_output, rr_table,
    los_r_col = "median_LOS_R", los_s_col = "median_LOS_S", los_unit = "days"
  )
  df <- out$Test_Pathogen

  # SS: all-susceptible -> no dominant class -> NA
  expect_true(is.na(df$LOSR_years[df$profile == "SS"]))
  expect_true(is.na(df$LOSS_years[df$profile == "SS"]))

  # SB: dominant class B -> LOSR = 15/365, LOSS = 10/365
  expect_equal(df$LOSR_years[df$profile == "SB"], round(15 / 365, 6L))
  expect_equal(df$LOSS_years[df$profile == "SB"], round(10 / 365, 6L))

  # SA: dominant class A -> LOSR = 20/365, LOSS = 10/365
  expect_equal(df$LOSR_years[df$profile == "SA"], round(20 / 365, 6L))
  expect_equal(df$LOSS_years[df$profile == "SA"], round(10 / 365, 6L))

  # AB: both resistant, A has higher RR (2.0 > 1.5) -> dominant class A
  expect_equal(df$dominant_class[df$profile == "AB"], "A")
  expect_equal(df$LOSR_years[df$profile == "AB"], round(20 / 365, 6L))
  expect_equal(df$LOSS_years[df$profile == "AB"], round(10 / 365, 6L))
})

test_that("daly_assign_rr_to_profiles() accepts LOS already in years without conversion", {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()
  rr_table$median_LOS_R_years <- rr_table$median_LOS_R / 365
  rr_table$median_LOS_S_years <- rr_table$median_LOS_S / 365

  out <- daly_assign_rr_to_profiles(
    profiles_output, rr_table,
    los_r_col = "median_LOS_R_years", los_s_col = "median_LOS_S_years", los_unit = "years"
  )
  df <- out$Test_Pathogen
  expect_equal(df$LOSR_years[df$profile == "SA"], round(20 / 365, 6L))
})

test_that("daly_assign_rr_to_profiles() requires los_r_col and los_s_col together", {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()
  expect_error(
    daly_assign_rr_to_profiles(profiles_output, rr_table, los_r_col = "median_LOS_R"),
    "los_r_col and los_s_col must both be supplied"
  )
})

.synthetic_profiles_with_los <- function() {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()
  daly_assign_rr_to_profiles(
    profiles_output, rr_table,
    los_r_col = "median_LOS_R", los_s_col = "median_LOS_S", los_unit = "days"
  )
}

test_that("daly_calc_yld_associated() matches the direct equation I_L * P'_Lk * DW_L * sum(R' * LOSR)", {
  profiles_with_los <- .synthetic_profiles_with_los()

  P_Lk_prime_tbl <- data.frame(
    pathogen = "Test_Pathogen",
    P_Lk_prime = 0.4,
    N_NF_L = 100,
    stringsAsFactors = FALSE
  )

  out <- daly_calc_yld_associated(
    profiles_with_los = profiles_with_los,
    P_Lk_prime_tbl = P_Lk_prime_tbl,
    DW_sepsis = 0.15
  )

  # Hand-computed expected value:
  # resistant profiles only: SB (p=0.20, LOSR=15/365), SA (p=0.20, LOSR=20/365),
  # AB (p=0.10, LOSR=20/365, dominant class A).
  losr_sb <- 15 / 365
  losr_sa <- 20 / 365
  losr_ab <- 20 / 365
  expected_profile_term <- round(0.20 * losr_sb + 0.20 * losr_sa + 0.10 * losr_ab, 6L)
  expected_yld <- 100 * 0.4 * 0.15 * expected_profile_term

  expect_equal(out$I_L, 100)
  expect_equal(out$DW_L, 0.15)
  expect_equal(out$profile_term, expected_profile_term)
  expect_equal(out$YLD_associated, expected_yld, tolerance = 1e-8)
})

test_that("daly_calc_yld_attributable() matches the direct equation I_L * P'_Lk * DW_L * sum(R' * (LOSR - LOSS))", {
  profiles_with_los <- .synthetic_profiles_with_los()

  P_Lk_prime_tbl <- data.frame(
    pathogen = "Test_Pathogen",
    P_Lk_prime = 0.4,
    N_NF_L = 100,
    stringsAsFactors = FALSE
  )

  out <- daly_calc_yld_attributable(
    profiles_with_los = profiles_with_los,
    P_Lk_prime_tbl = P_Lk_prime_tbl,
    DW_sepsis = 0.15
  )

  loss_sb <- 10 / 365
  loss_sa <- 10 / 365
  loss_ab <- 10 / 365
  losr_sb <- 15 / 365
  losr_sa <- 20 / 365
  losr_ab <- 20 / 365
  expected_profile_term <- round(
    0.20 * (losr_sb - loss_sb) +
      0.20 * (losr_sa - loss_sa) +
      0.10 * (losr_ab - loss_ab),
    6L
  )
  expected_yld <- 100 * 0.4 * 0.15 * expected_profile_term

  expect_equal(out$profile_term, expected_profile_term)
  expect_equal(out$YLD_attributable, expected_yld, tolerance = 1e-8)

  # Attributable must be strictly less than associated (LOSS > 0 always
  # subtracts a positive amount from every resistant profile's term).
  assoc <- daly_calc_yld_associated(profiles_with_los, P_Lk_prime_tbl, DW_sepsis = 0.15)
  expect_lt(out$YLD_attributable, assoc$YLD_associated)
})

test_that("daly_calc_yld_associated()/attributable() do not require a pre-computed baseline YLD or fraction/PAF object", {
  # The only inputs are profiles_with_los and P_Lk_prime_tbl (+ a DW source).
  # This test simply confirms the call succeeds without ever constructing
  # any kind of "baseline YLD" or "fraction associated" intermediate object.
  profiles_with_los <- .synthetic_profiles_with_los()
  P_Lk_prime_tbl <- data.frame(
    pathogen = "Test_Pathogen", P_Lk_prime = 0.4, N_NF_L = 100,
    stringsAsFactors = FALSE
  )
  expect_no_error(
    daly_calc_yld_associated(profiles_with_los, P_Lk_prime_tbl, DW_sepsis = 0.15)
  )
  expect_no_error(
    daly_calc_yld_attributable(profiles_with_los, P_Lk_prime_tbl, DW_sepsis = 0.15)
  )
})

test_that("daly_calc_yld_associated() warns when DW_L is sourced from yld_ref instead of a scalar", {
  profiles_with_los <- .synthetic_profiles_with_los()
  P_Lk_prime_tbl <- data.frame(
    pathogen = "Test_Pathogen", P_Lk_prime = 0.4, N_NF_L = 100,
    stringsAsFactors = FALSE
  )
  yld_ref <- data.frame(
    location_name = c("India", "Kerala"),
    DW_sepsis = c(0.36, 0.55),
    stringsAsFactors = FALSE
  )
  expect_warning(
    daly_calc_yld_associated(profiles_with_los, P_Lk_prime_tbl, yld_ref = yld_ref),
    "double-counting"
  )
})

test_that("daly_calc_yld_associated() warns and returns NA for pathogens missing from profiles_with_los", {
  profiles_with_los <- .synthetic_profiles_with_los()
  P_Lk_prime_tbl <- data.frame(
    pathogen = c("Test_Pathogen", "Unknown_Pathogen"),
    P_Lk_prime = c(0.4, 0.3),
    N_NF_L = c(100, 50),
    stringsAsFactors = FALSE
  )
  expect_warning(
    out <- daly_calc_yld_associated(profiles_with_los, P_Lk_prime_tbl, DW_sepsis = 0.15),
    "no matching entry in profiles_with_los"
  )
  expect_false(is.na(out$YLD_associated[out$pathogen == "Test_Pathogen"]))
  expect_true(is.na(out$YLD_associated[out$pathogen == "Unknown_Pathogen"]))
})

test_that("daly_calc_paf_los() still works standalone and is not required by the primary YLD pathway", {
  profiles_with_los <- .synthetic_profiles_with_los()

  # Works on its own, decoupled from daly_calc_yld_associated()/attributable().
  paf_out <- daly_calc_paf_los(profiles_with_los)
  expect_true("Test_Pathogen" %in% names(paf_out))
  expect_true(is.numeric(paf_out$Test_Pathogen$PAF_k))
  expect_gt(paf_out$Test_Pathogen$PAF_k, 0)
})

test_that("daly_calc_yld_baseline() and daly_calc_fraction_associated_yld() have been removed", {
  expect_false(exists("daly_calc_yld_baseline", where = asNamespace("anumaan"), inherits = FALSE))
  expect_false(exists("daly_calc_fraction_associated_yld", where = asNamespace("anumaan"), inherits = FALSE))
})

# ---- Profile strings, hospital-specific LOS and draw-based intervals -------

.fake_los_fit <- function() {
  # shared RR per class; hospital-specific adjusted LOS in days
  rr <- data.frame(
    pathogen = "K", antibiotic_class = c("A", "B"), hospital = NA_character_,
    RR_LOS = c(1.2, 1.5), CI_lower = c(1.1, 1.3), CI_upper = c(1.3, 1.7),
    LOSR = c(24, 30), LOSS = c(20, 20), stringsAsFactors = FALSE
  )
  hosp <- data.frame(
    pathogen = "K", antibiotic_class = rep(c("A", "B"), each = 2),
    hospital = rep(c("H1", "H2"), 2),
    LOSR = c(24, 12, 30, 15), LOSS = c(20, 10, 20, 10), stringsAsFactors = FALSE
  )
  set.seed(1)
  B <- 200
  mk <- function(r, s) {
    list(
      RR_LOS = rep(1, B),
      LOSR = cbind(H1 = r[1] + stats::rnorm(B), H2 = r[2] + stats::rnorm(B)),
      LOSS = cbind(H1 = s[1] + stats::rnorm(B), H2 = s[2] + stats::rnorm(B)),
      LOSR_all = rep(mean(r), B), LOSS_all = rep(mean(s), B)
    )
  }
  draws <- list("K | A" = mk(c(24, 12), c(20, 10)), "K | B" = mk(c(30, 15), c(20, 10)))
  list(rr = rr, hospital = hosp, draws = draws, skipped = NULL)
}

.fake_string_profiles <- function() {
  data.frame(
    group = rep(c("H1", "H2"), each = 4), pathogen = "K",
    profile_delta = rep(c("SS", "RS", "SR", "RR"), 2),
    R_mean = c(0.4, 0.2, 0.3, 0.1, 0.5, 0.2, 0.2, 0.1),
    stringsAsFactors = FALSE
  )
}

test_that("daly_assign_rr_to_profiles() handles profile strings with hospital-specific LOS", {
  fit <- .fake_los_fit()
  out <- suppressMessages(daly_assign_rr_to_profiles(
    .fake_string_profiles(), fit, profile_class_order = c("A", "B")
  ))$K
  expect_equal(out$dominant_class, rep(c("all_susceptible", "A", "B", "B"), 2))
  expect_equal(out$RR_LOS_profile, rep(c(1, 1.2, 1.5, 1.5), 2))
  # RR profile takes its own hospital's LOS of the dominant class, in years
  expect_equal(out$LOSR_years[out$hospital == "H2" & out$profile_delta == "RR"], round(15 / 365, 6))
  expect_equal(out$LOSS_years[out$hospital == "H1" & out$profile_delta == "RS"], round(20 / 365, 6))
  expect_true(all(is.na(out$LOSR_years[out$dominant_class == "all_susceptible"])))
  expect_true(all(out$los_level[out$dominant_class != "all_susceptible"] == "hospital"))
  expect_equal(out$probability, out$R_mean)
})

test_that("daly_assign_rr_to_profiles() breaks ties by class order and checks profile length", {
  fit <- .fake_los_fit()
  fit$rr$RR_LOS <- c(1.5, 1.5)
  expect_message(
    out <- daly_assign_rr_to_profiles(.fake_string_profiles(), fit, profile_class_order = c("A", "B")),
    "tied"
  )
  expect_equal(out$K$dominant_class[out$K$profile_delta == "RR"], c("A", "A"))
  expect_error(
    suppressMessages(daly_assign_rr_to_profiles(.fake_string_profiles(), fit, profile_class_order = c("A", "B", "C"))),
    "does not match"
  )
  expect_error(daly_assign_rr_to_profiles(.fake_string_profiles(), fit), "profile_class_order")
})

test_that("daly_calc_yld_attributable() matches each hospital to its own profiles and LOS", {
  fit <- .fake_los_fit()
  prof <- suppressMessages(daly_assign_rr_to_profiles(.fake_string_profiles(), fit, profile_class_order = c("A", "B")))
  plk <- data.frame(pathogen = "K", center = c("H1", "H2"), P_Lk_prime = 0.5, N_NF_L = c(100, 200))
  out <- suppressMessages(daly_calc_yld_attributable(prof, plk, DW_sepsis = 0.1, facility_col = "center"))
  # H1: 0.2*(24-20) + 0.3*(30-20) + 0.1*(30-20) = 4.8 days
  expect_equal(out$YLD_attributable[1], 100 * 0.5 * 0.1 * 4.8 / 365, tolerance = 1e-3)
  # H2: 0.2*(12-10) + 0.2*(15-10) + 0.1*(15-10) = 1.9 days
  expect_equal(out$YLD_attributable[2], 200 * 0.5 * 0.1 * 1.9 / 365, tolerance = 1e-3)
  expect_error(
    suppressMessages(daly_calc_yld_attributable(prof, plk, DW_sepsis = 0.1)),
    "facility_col"
  )
})

test_that("daly_calc_yld_*() return draw-based intervals from LOS and profile draws", {
  fit <- .fake_los_fit()
  prof <- suppressMessages(daly_assign_rr_to_profiles(.fake_string_profiles(), fit, profile_class_order = c("A", "B")))
  plk <- data.frame(pathogen = "K", center = c("H1", "H2"), P_Lk_prime = 0.5, N_NF_L = c(100, 200))
  a <- suppressMessages(daly_calc_yld_associated(prof, plk, DW_sepsis = 0.1, facility_col = "center",
                                                 los_draws = fit$draws))
  expect_true(all(a$YLD_associated_lower < a$YLD_associated & a$YLD_associated < a$YLD_associated_upper))
  expect_equal(a$n_draws, c(200L, 200L))
  expect_equal(a$uncertainty_sources, c("LOS only", "LOS only"))

  sp <- .fake_string_profiles()
  pdraws <- do.call(rbind, lapply(1:50, function(b) {
    data.frame(pathogen = "K", hospital = sp$group, profile = sp$profile_delta, draw = b,
               probability = sp$R_mean)
  }))
  b <- suppressMessages(daly_calc_yld_attributable(prof, plk, DW_sepsis = 0.1, facility_col = "center",
                                                   los_draws = fit$draws, profile_draws = pdraws))
  expect_equal(b$n_draws, c(50L, 50L))
  expect_equal(b$uncertainty_sources, c("LOS + profile", "LOS + profile"))
})
