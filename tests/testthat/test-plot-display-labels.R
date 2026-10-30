# Tests for R/plot_display_labels.R -- generic display-label helpers used by
# PPC/validation reports (hospital/class name formatting). The PPC-plot
# structural tests that used to live below these (exercising the probit
# posterior-predictive plotting internals in R/probit_posterior_predictive.R)
# were removed along with that file.

test_that("hospital_display_label() replaces underscores with spaces, nothing else", {
  expect_identical(hospital_display_label("AIIMS_trauma_center"), "AIIMS trauma center")
  expect_identical(hospital_display_label("Assam_Medical_College"), "Assam Medical College")
  expect_identical(hospital_display_label("SGRH"), "SGRH")
})

test_that("class_display_label() looks up known classes and falls back gracefully for unknown ones", {
  expect_identical(class_display_label("Third_generation_cephalosporins"), "Third-generation cephalosporins")
  expect_identical(class_display_label("Aminoglycosides"), "Aminoglycosides")
  # Unmapped class: falls back to underscore-to-space, never errors.
  expect_identical(class_display_label("Some_Future_Class"), "Some Future Class")
})

test_that("class_short_label() looks up known abbreviations and falls back gracefully", {
  expect_identical(class_short_label("Third_generation_cephalosporins"), "3GC")
  expect_identical(class_short_label("Aminoglycosides"), "AMG")
  expect_identical(class_short_label("Carbapenems"), "CARB")
  expect_identical(class_short_label("Fluoroquinolones"), "FQ")
  expect_identical(class_short_label("Penicillins"), "PEN")
  # Unmapped: 4-char upper-cased fallback, never errors.
  expect_identical(class_short_label("Zzznotarealclass"), "ZZZN")
})

test_that("class_pair_label() combines two classes with 'x' using short labels by default", {
  expect_identical(class_pair_label("Aminoglycosides", "Carbapenems"), "AMG x CARB")
  expect_identical(
    class_pair_label("Aminoglycosides", "Carbapenems", short = FALSE),
    "Aminoglycosides x Carbapenems"
  )
})
