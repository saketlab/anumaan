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
    mean_LOS_R = c(20, 15), # days
    mean_LOS_S = c(10, 10), # days
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
    los_r_col = "mean_LOS_R", los_s_col = "mean_LOS_S", los_unit = "days"
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
  rr_table$mean_LOS_R_years <- rr_table$mean_LOS_R / 365
  rr_table$mean_LOS_S_years <- rr_table$mean_LOS_S / 365

  out <- daly_assign_rr_to_profiles(
    profiles_output, rr_table,
    los_r_col = "mean_LOS_R_years", los_s_col = "mean_LOS_S_years", los_unit = "years"
  )
  df <- out$Test_Pathogen
  expect_equal(df$LOSR_years[df$profile == "SA"], round(20 / 365, 6L))
})

test_that("daly_assign_rr_to_profiles() requires los_r_col and los_s_col together", {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()
  expect_error(
    daly_assign_rr_to_profiles(profiles_output, rr_table, los_r_col = "mean_LOS_R"),
    "los_r_col and los_s_col must both be supplied"
  )
})

.synthetic_profiles_with_los <- function() {
  profiles_output <- .synthetic_profiles_output()
  rr_table <- .synthetic_rr_table_days()
  daly_assign_rr_to_profiles(
    profiles_output, rr_table,
    los_r_col = "mean_LOS_R", los_s_col = "mean_LOS_S", los_unit = "days"
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
