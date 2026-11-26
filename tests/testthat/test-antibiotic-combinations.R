# prep_standardize_antibiotics: combination drugs must not collapse to one component

combo_who <- function() {
  data.frame(
    Antibiotic = c("Piperacillin", "Piperacillin/Tazobactam", "Trimethoprim",
                   "Trimethoprim/Sulfamethoxazole", "Sulfamethoxazole", "Ampicillin",
                   "Ampicillin/Sulbactam", "Sulbactam", "Cefoperazone", "Gentamicin",
                   "Benzylpenicillin", "Fusidic-acid", "Imipenem/cilastatin"),
    Class = c("Penicillins", "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal",
              "Trimethoprim-derivatives", "Sulfonamide-trimethoprim-combinations", "Sulfonamides",
              "Aminopenicillins", "Beta-lactam/beta-lactamase-inhibitor", "Beta-lactamase-inhibitors",
              "Third-generation-cephalosporins", "Aminoglycosides", "Penicillins",
              "Steroid antibacterials", "Carbapenems"),
    Category = c("Watch", "Watch", "Access", "Access", "Access", "Access", "Access", "Access",
                 "Watch", "Access", "Access", "Watch", "Watch"),
    stringsAsFactors = FALSE
  )
}

std <- function(x, who = combo_who()) {
  suppressMessages(prep_standardize_antibiotics(
    data.frame(antibiotic_name = x, stringsAsFactors = FALSE), who_table = who))
}

test_that("hyphen, plus and slash combinations match the combination row", {
  out <- std(c("piperacillin-tazobactam", "Piperacillin/Tazobactam", "PIPERACILLIN+TAZOBACTAM",
               "trimethoprim-sulfamethoxazole", "ampicillin-sulbactam"))
  expect_equal(out$antibiotic_normalized,
               c(rep("piperacillin/tazobactam", 3), "trimethoprim/sulfamethoxazole", "ampicillin/sulbactam"))
  expect_equal(out$antibiotic_class,
               c(rep("Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal", 3),
                 "Sulfonamide-trimethoprim-combinations", "Beta-lactam/beta-lactamase-inhibitor"))
})

test_that("component order and one-letter misspellings still match the combination", {
  out <- std(c("sulfamethoxazole/trimethoprim", "Piperacillin+Tazobactum", "tazobactam and piperacillin"))
  expect_equal(out$antibiotic_normalized,
               c("trimethoprim/sulfamethoxazole", "piperacillin/tazobactam", "piperacillin/tazobactam"))
})

test_that("single drugs keep mapping to their own rows", {
  out <- std(c("piperacillin", "trimethoprim", "ampicillin", "gentamicin", "penicillin", "imipenem"))
  expect_equal(out$antibiotic_normalized,
               c("piperacillin", "trimethoprim", "ampicillin", "gentamicin", "benzylpenicillin", "imipenem/cilastatin"))
  expect_equal(out$antibiotic_class[1], "Penicillins")
})

test_that("high-level gentamicin maps to gentamicin", {
  out <- std(c("gentamicin hl", "gentamicin_hl"))
  expect_equal(out$antibiotic_normalized, c("gentamicin", "gentamicin"))
  expect_equal(out$antibiotic_class, c("Aminoglycosides", "Aminoglycosides"))
})

test_that("generic word 'acid' does not drive a match", {
  out <- std("nalidixic acid")
  expect_false(identical(out$antibiotic_normalized, "fusidic-acid"))
})

test_that("a combination with no reference row keeps its full name, NA class, and warns", {
  expect_warning(
    out <- suppressMessages(prep_standardize_antibiotics(
      data.frame(antibiotic_name = "cefoperazone-sulbactam"), who_table = combo_who())),
    "no reference row"
  )
  expect_equal(out$antibiotic_normalized, "cefoperazone/sulbactam")
  expect_true(is.na(out$antibiotic_class))
  expect_true(is.na(out$aware_category))
})

test_that("default reference includes the combination supplement", {
  skip_if(find_extdata_file("WHO_aware_class.csv") == "")
  out <- suppressMessages(prep_standardize_antibiotics(data.frame(
    antibiotic_name = c("cefoperazone-sulbactam", "ticarcillin-clavulanic acid",
                        "ceftazidime-avibactam", "meropenem-vaborbactam"))))
  expect_equal(out$antibiotic_normalized,
               c("cefoperazone/sulbactam", "ticarcillin/clavulanic acid",
                 "ceftazidime/avibactam", "meropenem/vaborbactam"))
  expect_equal(out$antibiotic_class,
               c("Third-generation-cephalosporins", "Beta-lactam/beta-lactamase-inhibitor",
                 "Third-generation-cephalosporins", "Carbapenems"))
  expect_equal(out$aware_category[c(1, 3, 4)], c("Not recommended", "Reserve", "Reserve"))
  expect_equal(nrow(out), 4L)
})
