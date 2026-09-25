# Tests for the HAI/CAI classification fix in daly_fit_los_rr() and
# daly_fit_los_rr_distribution() (R/daly_rr_and_los.R).
#
# Both functions used to unconditionally call daly_derive_hai_cai_for_los(),
# re-deriving HAI/CAI from the admission-to-culture date gap even when a
# final classification (e.g. from prep_derive_hai_cai() in the
# anumaan-analysis preprocessing pipeline) was already present in the data.
# They now read infection_type_col directly and no longer call
# daly_derive_hai_cai_for_los() at all.
#
# Verification strategy: build patients with an IDENTICAL admission-to-culture
# date gap (so date-based re-derivation, if it still ran, would classify every
# patient the same way), but explicitly set infection_type_col to different
# values across two runs. Since HAI's LOS clock starts at the culture date and
# CAI's starts at the admission date, and the two dates differ by a known,
# fixed number of days, the resulting mean LOS must shift by exactly that gap
# between the two runs if (and only if) infection_type_col is being honoured
# directly rather than being overridden by date-based re-derivation (which
# would produce IDENTICAL results in both runs, since the date gap never
# changes).

.synthetic_los_data <- function(infection_type_value) {
  set.seed(42)
  n <- 30
  admission <- as.Date("2026-01-01")
  culture <- admission + 7 # fixed 7-day gap, every patient
  # Jittered discharge offsets, independent of R/S assignment below, so each
  # group has non-degenerate within-group variance for the gamma fit.
  discharge <- culture + sample(15:25, n, replace = TRUE)
  data.frame(
    PatientInformation_id = as.character(seq_len(n)),
    center_name = "Hosp1",
    organism_name = "Test_Pathogen",
    syndrome = "BSI",
    infection_type = infection_type_value, # pre-computed, already final
    antibiotic_class = "A",
    antibiotic_name = "drugA",
    antibiotic_value = sample(c("R", "S"), n, replace = TRUE),
    date_of_admission = admission,
    date_of_first_positive_culture = culture,
    final_outcome_date = discharge,
    final_outcome = "Discharged",
    stringsAsFactors = FALSE
  )
}

test_that("daly_fit_los_rr_distribution() honours infection_type_col directly (HAI vs CAI shifts mean LOS by exactly the admission-culture gap)", {
  hai_data <- .synthetic_los_data("HAI") # clock starts at culture date
  cai_data <- .synthetic_los_data("CAI") # clock starts at admission date

  common_args <- list(
    patient_id_col = "PatientInformation_id",
    facility_col = "center_name",
    organism_col = "organism_name",
    syndrome_col = "syndrome",
    infection_type_col = "infection_type",
    antibiotic_class_col = "antibiotic_class",
    antibiotic_name_col = "antibiotic_name",
    antibiotic_value_col = "antibiotic_value",
    date_admission_col = "date_of_admission",
    date_discharge_col = "final_outcome_date",
    date_culture_col = "date_of_first_positive_culture",
    final_outcome_col = "final_outcome",
    min_n = 5
  )

  out_hai <- do.call(daly_fit_los_rr_distribution, c(list(data = hai_data), common_args))
  out_cai <- do.call(daly_fit_los_rr_distribution, c(list(data = cai_data), common_args))

  expect_equal(nrow(out_hai), 1L)
  expect_equal(nrow(out_cai), 1L)

  # Admission is 7 days before culture, so CAI's clock (starts at admission)
  # yields LOS values exactly 7 days longer than HAI's clock (starts at
  # culture) for the SAME underlying discharge dates. If infection_type were
  # still being re-derived from the (identical, in this fixture) date gap
  # instead of read directly, out_hai and out_cai would be identical.
  expect_equal(out_cai$mean_LOS_R - out_hai$mean_LOS_R, 7, tolerance = 0.5)
  expect_equal(out_cai$mean_LOS_S - out_hai$mean_LOS_S, 7, tolerance = 0.5)
})

test_that("daly_fit_los_rr_distribution() defaults infection_type_col to 'infection_type'", {
  expect_equal(
    formals(daly_fit_los_rr_distribution)$infection_type_col,
    "infection_type"
  )
})

test_that("daly_fit_los_rr() defaults infection_type_col to 'infection_type'", {
  expect_equal(
    formals(daly_fit_los_rr)$infection_type_col,
    "infection_type"
  )
})

test_that("daly_fit_los_rr_distribution() errors clearly when infection_type_col is missing", {
  bad_data <- .synthetic_los_data("HAI")
  bad_data$infection_type <- NULL
  expect_error(
    daly_fit_los_rr_distribution(
      bad_data,
      patient_id_col = "PatientInformation_id",
      facility_col = "center_name",
      organism_col = "organism_name",
      syndrome_col = "syndrome",
      infection_type_col = "infection_type",
      antibiotic_class_col = "antibiotic_class",
      antibiotic_name_col = "antibiotic_name",
      antibiotic_value_col = "antibiotic_value",
      date_admission_col = "date_of_admission",
      date_discharge_col = "final_outcome_date",
      date_culture_col = "date_of_first_positive_culture",
      final_outcome_col = "final_outcome"
    ),
    "infection_type"
  )
})
