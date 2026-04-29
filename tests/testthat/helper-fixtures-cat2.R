# Category 2 shared fixtures
#
# These fixtures are synthetic and source-inspired. They mimic the shapes and
# signal combinations seen in the temporary hospital extracts under `Data/`
# without copying real patient IDs, dates, or centre-specific identifiers.


# ---------------------------------------------------------------------------
# Intake / schema
# ---------------------------------------------------------------------------

cat2_intake_minimal_valid <- function() {
  data.frame(
    patient_id = paste0("pt_", sprintf("%03d", 1:20)),
    date_of_admission = as.Date(c(
      "2025-01-01", "2025-01-03", "2025-01-05", "2025-01-06", "2025-01-08",
      "2025-01-10", "2025-01-12", "2025-01-15", "2025-01-18", "2025-01-20",
      "2025-01-22", "2025-01-25", "2025-01-28", "2025-02-01", "2025-02-03",
      "2025-02-05", "2025-02-08", "2025-02-10", "2025-02-12", "2025-02-15"
    )),
    date_of_culture = as.Date(c(
      "2025-01-01", "2025-01-05", "2025-01-06", "2025-01-08", "2025-01-09",
      "2025-01-11", "2025-01-14", "2025-01-17", "2025-01-19", "2025-01-22",
      "2025-01-24", "2025-01-26", "2025-01-30", "2025-02-02", "2025-02-05",
      "2025-02-06", "2025-02-09", "2025-02-12", "2025-02-13", "2025-02-17"
    )),
    diagnosis_time_admission = as.POSIXct(c(
      "2025-01-01 08:00:00", "2025-01-03 09:30:00", "2025-01-05 10:15:00",
      "2025-01-06 14:00:00", "2025-01-08 11:20:00", "2025-01-10 07:45:00",
      "2025-01-12 13:10:00", "2025-01-15 16:30:00", "2025-01-18 06:50:00",
      "2025-01-20 12:00:00", "2025-01-22 08:40:00", "2025-01-25 09:05:00",
      "2025-01-28 15:15:00", "2025-02-01 10:50:00", "2025-02-03 11:35:00",
      "2025-02-05 17:25:00", "2025-02-08 08:55:00", "2025-02-10 14:10:00",
      "2025-02-12 09:45:00", "2025-02-15 07:35:00"
    ), tz = "UTC"),
    organism_name = c(
      "Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus",
      "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterococcus faecalis",
      "Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus",
      "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterobacter cloacae",
      "Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus",
      "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Enterococcus faecium",
      "Escherichia coli", "Streptococcus pneumoniae"
    ),
    antibiotic_name = c(
      "Amikacin", "Cefotaxime", "Vancomycin", "Colistin", "Piperacillin/Tazobactam",
      "Ampicillin", "Meropenem", "Ciprofloxacin", "Linezolid", "Tigecycline",
      "Gentamicin", "Cefepime", "Amikacin", "Ceftriaxone", "Vancomycin",
      "Colistin", "Levofloxacin", "Teicoplanin", "Meropenem", "Penicillin"
    ),
    antibiotic_value = c("R", "S", "S", "I", "R", "S", "R", "S", "S", "I",
                         "R", "S", "S", "R", "S", "I", "R", "S", "S", "S"),
    final_outcome = c(
      "Survived", "Died", "Survived", "Survived", "Died",
      "Survived", "Survived", "Died", "Survived", "Survived",
      "Died", "Survived", "Survived", "Died", "Survived",
      "Survived", "Died", "Survived", "Survived", "Died"
    ),
    final_outcome_date = as.Date(c(
      "2025-01-05", "2025-01-08", "2025-01-11", "2025-01-12", "2025-01-15",
      "2025-01-16", "2025-01-18", "2025-01-22", "2025-01-24", "2025-01-27",
      "2025-01-29", "2025-02-01", "2025-02-03", "2025-02-06", "2025-02-08",
      "2025-02-10", "2025-02-13", "2025-02-15", "2025-02-17", "2025-02-20"
    )),
    location = c(
      "ICU", "Ward", "OPD", "ICU", "Ward",
      "OPD", "ICU", "Ward", "OPD", "ICU",
      "Ward", "OPD", "ICU", "Ward", "OPD",
      "ICU", "Ward", "OPD", "ICU", "Ward"
    ),
    previous_hospitalisation = c(
      "Yes", "No", "No", "Yes", "No",
      "Unknown", "Yes", "No", "Yes", "No",
      "No", "Yes", "No", "Unknown", "Yes",
      "No", "Yes", "No", "No", "Yes"
    ),
    stringsAsFactors = FALSE
  )
}

cat2_intake_bad_keys <- function() {
  data.frame(
    patient_id = c("pt_001", "pt_001", "", "NA", NA_character_),
    date_of_culture = c("2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04", NA_character_),
    organism_name = c("Escherichia coli", "Escherichia coli", "Klebsiella pneumoniae", "", NA_character_),
    stringsAsFactors = FALSE
  )
}

cat2_schema_aliases <- function() {
  data.frame(
    PID = c("pt_001", "pt_002", "pt_003"),
    CultureDate = c("2025-01-01", "2025-01-05", "2025-01-06"),
    Organism = c("E. coli", "K. pneumoniae", "MRSA"),
    Sample_Type = c("Blood culture", "Urine culture / sensitivity", "Blood"),
    Ward_ICU = c("Ward", "ICU", "critical care"),
    stringsAsFactors = FALSE
  )
}

cat2_intake_date_edge_cases <- function() {
  data.frame(
    patient_id = c("pt_edge_001", "pt_edge_002", "pt_edge_003", "pt_edge_004"),
    date_of_admission = c("43831", "20201301", "205-19-10", "2025-02"),
    date_of_culture = c("43832", "20201401", "2025/01/20", "2025-02-15"),
    final_outcome_date = c("43840", "20201501", "", "2025-02-20"),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Standardization
# ---------------------------------------------------------------------------

cat2_standardization_organism_abx <- function() {
  data.frame(
    patient_id = paste0("pt_", sprintf("%03d", 1:7)),
    organism_name = c(
      "E. coli",
      "MRSA(Methicillin resistant staphylococcus aureus)",
      "CONS (Coagulase Negative Staphylococci)",
      "Non fermenting Gram negative bacilli",
      "No growth in culture.",
      "Bacillus spp. grown.",
      "Unknown environmental isolate"
    ),
    specimen_type = c("Blood", "Blood", "Blood", "ETA", "Blood", "Blood", "Urine"),
    antibiotic_name = c(
      "Amikacin",
      "Vancomycin",
      "Oxacillin",
      "Colistin",
      "Ampicillin",
      "Ciproflox",
      "UnknownDrugZZ"
    ),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# AST
# ---------------------------------------------------------------------------

cat2_ast_long_clean <- function() {
  data.frame(
    patient_id = c("pt_001", "pt_001", "pt_002", "pt_003"),
    event_id = c("ev_001", "ev_001", "ev_002", "ev_003"),
    culture_date = as.Date(c("2025-01-01", "2025-01-01", "2025-01-05", "2025-01-06")),
    organism_normalized = c(
      "escherichia coli",
      "escherichia coli",
      "klebsiella pneumoniae",
      "staphylococcus aureus"
    ),
    organism_group = c("Enterobacterales", "Enterobacterales", "Enterobacterales", "Gram-positive cocci"),
    antibiotic_normalized = c("amikacin", "cefotaxime", "ciprofloxacin", "vancomycin"),
    antibiotic_name = c("Amikacin", "Cefotaxime", "Ciprofloxacin", "Vancomycin"),
    antibiotic_class = c("Aminoglycosides", "Third-generation-cephalosporins", "Fluoroquinolones", "Glycopeptides"),
    antibiotic_value = c("R", "S", "S", "S"),
    ast_value_raw = c("Resistant", "Susceptible", "S", "Sensitive"),
    ast_value_harmonized = c("R", "S", "S", "S"),
    stringsAsFactors = FALSE
  )
}

cat2_ast_long_dirty <- function() {
  data.frame(
    patient_id = c("pt_010", "pt_010", "pt_011", "pt_012", "pt_013", "pt_014"),
    event_id = c("ev_010", "ev_010", "ev_011", "ev_012", "ev_013", "ev_014"),
    culture_date = as.Date(c("2025-01-08", "2025-01-08", "2025-01-09", "2025-01-10", "2025-01-11", "2025-01-12")),
    organism_normalized = c(
      "escherichia coli",
      "escherichia coli",
      "klebsiella pneumoniae",
      "acinetobacter baumannii",
      "staphylococcus aureus",
      "enterococcus faecium"
    ),
    organism_group = c("Enterobacterales", "Enterobacterales", "Enterobacterales", "Gram-negative bacilli", "Gram-positive cocci", "Gram-positive cocci"),
    antibiotic_normalized = c("amikacin", "colistin", "ciprofloxacin", "cefotaxime", "oxacillin", "linezolid"),
    antibiotic_name = c("Amikacin", "Colistin", "Ciprofloxacin", "Cefotaxime", "Oxacillin", "Linezolid"),
    antibiotic_class = c("Aminoglycosides", "Polymyxins", "Fluoroquinolones", "Third-generation-cephalosporins", "Penicillins", "Oxazolidinones"),
    antibiotic_value = c("Resistant", "Intermediate", "No conclusion", "1", "", "X"),
    ast_value_raw = c("Resistant", "Intermediate", "No conclusion", "1", "", "X"),
    ast_value_harmonized = c("R", "I", NA_character_, NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
}

cat2_ast_duplicates_same_day <- function() {
  data.frame(
    patient_id = c("pt_020", "pt_020", "pt_020", "pt_020", "pt_021"),
    event_id = c("ev_020", "ev_020", "ev_020", "ev_020", "ev_021"),
    culture_date = as.Date(c("2025-01-15", "2025-01-15", "2025-01-15", "2025-01-15", "2025-01-16")),
    organism_normalized = c(
      "escherichia coli",
      "escherichia coli",
      "escherichia coli",
      "escherichia coli",
      "klebsiella pneumoniae"
    ),
    organism_group = c("Enterobacterales", "Enterobacterales", "Enterobacterales", "Enterobacterales", "Enterobacterales"),
    antibiotic_normalized = c("amikacin", "amikacin", "amikacin", "cefotaxime", "ampicillin"),
    antibiotic_name = c("Amikacin", "Amikacin", "Amikacin", "Cefotaxime", "Ampicillin"),
    antibiotic_class = c("Aminoglycosides", "Aminoglycosides", "Aminoglycosides", "Third-generation-cephalosporins", "Aminopenicillins"),
    susceptibility = c("DISK", "MIC", "DISK", "DISK", "DISK"),
    antibiotic_value = c("S", "S", "R", "I", "S"),
    ast_value_raw = c("S", "S", "R", "I", "S"),
    ast_value_harmonized = c("S", "S", "R", "I", "S"),
    stringsAsFactors = FALSE
  )
}

cat2_ast_wide_matrix <- function() {
  data.frame(
    patient_id = c("pt_001", "pt_002", "pt_003"),
    organism_name = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus"),
    AMK = c("R", "S", NA),
    CIP = c("S", "R", NA),
    CTX = c("I", NA, NA),
    COL = c(NA, "S", NA),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Enrichment
# ---------------------------------------------------------------------------

cat2_enrichment_age_variants <- function() {
  data.frame(
    patient_id = c("pt_101", "pt_102", "pt_103", "pt_104"),
    Age = c(NA_real_, NA_real_, 5, NA_real_),
    DOB = as.Date(c("2010-01-15", NA, "2020-06-01", NA)),
    date_of_culture = as.Date(c("2025-01-15", "2025-01-20", "2025-01-25", "2025-01-25")),
    age_years = c(NA_real_, 2, NA_real_, 0),
    age_months = c(NA_real_, 6, NA_real_, 3),
    age_days = c(NA_real_, 0, NA_real_, 10),
    specimen_type = c("blood", "blood", "peritoneal abscess", "csf"),
    diagnosis_1 = c("sepsis", "", "abdomen pain", "meningitis"),
    hospital_department = c(NA_character_, NA_character_, "General Surgery", NA_character_),
    admission_date = as.Date(c("2025-01-10", "2025-01-18", "2025-01-20", "2025-01-20")),
    stringsAsFactors = FALSE
  )
}

cat2_enrichment_department_variants <- function() {
  data.frame(
    Age = c(5, 40, 67, 22, NA),
    specimen_type = c("blood", "peritoneal abscess", "urine", "csf", "blood"),
    diagnosis_1 = c("", "", "urosepsis", "meningitis", ""),
    hospital_department = c(NA_character_, NA_character_, "MEDICINE", NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
}

cat2_enrichment_los_variants <- function() {
  data.frame(
    admission_date = as.Date(c("2025-01-01", "2025-01-05", NA, "2025-01-08")),
    outcome_date = as.Date(c("2025-01-03", "2025-01-05", NA, "2025-01-06")),
    unit_admission_date = as.Date(c(NA, NA, "2025-01-09", NA)),
    unit_duration_days = c(NA_real_, NA_real_, 3, NA_real_),
    stringsAsFactors = FALSE
  )
}

cat2_enrichment_los_date_edge_cases <- function() {
  data.frame(
    patient_id = c("pt_los_001", "pt_los_002", "pt_los_003", "pt_los_004", "pt_los_005"),
    admission_date = c("43831", "20201301", "205-19-10", "2020-01-05", "2020-01-10"),
    outcome_date = c("43833", "20201701", "2020-01-20", "2020-01-04", "2020-01-10"),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# HAI / CAI
# ---------------------------------------------------------------------------

cat2_hai_gap_logic <- function() {
  data.frame(
    patient_id = c("pt_201", "pt_202", "pt_203", "pt_204", "pt_205", "pt_206"),
    date_of_admission = as.Date(c("2025-01-01", "2025-01-01", "2025-01-01", "2025-01-01", NA, "2025-01-01")),
    date_of_culture = as.Date(c("2025-01-01", "2025-01-02", "2025-01-03", "2025-01-06", "2025-01-05", NA)),
    infection_type = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
}

cat2_hai_observed_labels <- function() {
  data.frame(
    patient_id = c("pt_211", "pt_212", "pt_213", "pt_214", "pt_215"),
    type_of_infection = c(
      "Community acquired infection",
      "Health care associated infection",
      "Not known",
      "NULL",
      ""
    ),
    infection_type = c("CAI", "HAI", NA_character_, NA_character_, NA_character_),
    infection_type_method = c("provided", "provided", NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
}

cat2_hai_observed_vs_inferred <- function() {
  data.frame(
    patient_id = c("pt_221", "pt_222", "pt_223", "pt_224"),
    infection_type_observed = c("HAI", "CAI", NA_character_, "CAI"),
    infection_type_inferred = c("CAI", "HAI", "HAI", "CAI"),
    stringsAsFactors = FALSE
  )
}

cat2_icu_location_variants <- function() {
  data.frame(
    ward_icu = c("ICU", "Ward", "critical care", NA, "Multi disciplinary-ICU", "4F WARD", NA),
    hospital_department = c("", "General", NA, "Intensive Therapy", NA, "Medicine", NA),
    stringsAsFactors = FALSE
  )
}

cat2_hai_realistic_integrated <- function() {
  data.frame(
    patient_id = c("pt_231", "pt_232", "pt_233", "pt_234", "pt_235", "pt_236"),
    PatientInformation_id = 501:506,
    type_of_infection = c(
      "Community acquired infection",
      "Health care associated infection",
      "Not known",
      "NULL",
      "",
      "Health care associated infection"
    ),
    infection_type = c("CAI", "HAI", NA_character_, NA_character_, NA_character_, "HAI"),
    date_of_admission = as.Date(c("2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04", "2025-01-05", "2025-01-06")),
    date_of_culture = as.Date(c("2025-01-01", "2025-01-07", "2025-01-04", "2025-01-08", "2025-01-05", "2025-01-10")),
    date_HAI = as.Date(c(NA, "2025-01-06", NA, "2025-01-07", NA, "2025-01-09")),
    nature_of_hai = c(
      NA_character_,
      "Catheter associated blood stream infection",
      NA_character_,
      "Ventilator Associated Pneumonia",
      NA_character_,
      "Catheter associated urinary tract infection"
    ),
    hai_yes_no = c(NA_character_, "Yes", NA_character_, "Yes", NA_character_, "Yes"),
    any_device_inserted = c(
      NA_character_,
      "Central line catheter",
      NA_character_,
      "Ventilator",
      NA_character_,
      "Urinary catheter"
    ),
    central_line_days = c(NA_real_, 4, NA_real_, NA_real_, NA_real_, NA_real_),
    ventilator_days = c(NA_real_, NA_real_, NA_real_, 6, NA_real_, NA_real_),
    urinary_days = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 5),
    location = c("Ward", "ICU", "Ward", "ICU", "Ward", "Ward"),
    nameicu = c(NA_character_, "Multi disciplinary-ICU", NA_character_, "Medicine-ICU", NA_character_, NA_character_),
    nameofward = c("Male ward", NA_character_, "4F WARD", NA_character_, "Female ward", "WARD-03"),
    hospital_department = c("Medicine", "Medicine", "General Surgery", "Pulmonology", "Pediatrics", "Urology"),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Events / polymicrobial / outputs
# ---------------------------------------------------------------------------

cat2_events_episode_timing <- function() {
  data.frame(
    patient_id = c("pt_301", "pt_301", "pt_301", "pt_302", "pt_302"),
    date_of_culture = as.Date(c("2025-01-01", "2025-01-01", "2025-01-20", "2025-01-05", NA)),
    admission_date = as.Date(c("2025-01-01", "2025-01-10", "2025-05-01", "2025-01-05", "2025-03-20")),
    organism_normalized = c("escherichia coli", "escherichia coli", "escherichia coli", "klebsiella pneumoniae", "klebsiella pneumoniae"),
    specimen_type = c("blood", "blood", "blood", "urine", "urine"),
    antibiotic_name = c("amikacin", "amikacin", "amikacin", "ampicillin", "ampicillin"),
    antibiotic_value = c("S", "S", "R", "S", "S"),
    stringsAsFactors = FALSE
  )
}

cat2_polymicrobial_episode_mix <- function() {
  data.frame(
    episode_id = c("ep_401", "ep_402", "ep_402", "ep_403", "ep_403", "ep_403"),
    patient_id = c("pt_401", "pt_402", "pt_402", "pt_403", "pt_403", "pt_403"),
    organism_normalized = c(
      "escherichia coli",
      "klebsiella pneumoniae",
      "enterococcus faecalis",
      "acinetobacter baumannii",
      "pseudomonas aeruginosa",
      "staphylococcus aureus"
    ),
    is_polymicrobial = c(0, 1, 1, 1, 1, 1),
    stringsAsFactors = FALSE
  )
}

cat2_outputs_ready <- function() {
  data.frame(
    patient_id = c("pt_601", "pt_602", "pt_603", "pt_604"),
    culture_date = as.Date(c("2025-01-01", "2025-01-02", NA, "2025-01-04")),
    organism_name = c("Escherichia coli", "Klebsiella pneumoniae", "", "Staphylococcus aureus"),
    antibiotic_name = c("Amikacin", "Ampicillin", "Amikacin", ""),
    ast_value_harmonized = c("R", "S", NA, "I"),
    final_outcome = c("Died", "Survived", "Survived", "Survived"),
    contaminant_flag = c(FALSE, FALSE, TRUE, FALSE),
    event_id = c("ev_601", "ev_602", "ev_603", "ev_604"),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Backward-compatible wrappers used by the current Cat 2 test scripts
# ---------------------------------------------------------------------------

cat2_intake_df <- function() {
  cat2_intake_minimal_valid()
}

cat2_schema_df <- function() {
  cat2_schema_aliases()
}

cat2_ast_long_df <- function() {
  cat2_ast_duplicates_same_day()
}

cat2_enrich_df <- function() {
  dat <- cat2_enrichment_age_variants()[1:2, , drop = FALSE]
  dat$DOB[2] <- as.Date(NA)
  dat$date_of_culture[2] <- as.Date("2025-01-20")
  dat$age_years[2] <- 35
  dat$age_months[2] <- 0
  dat$age_days[2] <- 0
  dat
}

cat2_events_df <- function() {
  cat2_events_episode_timing()
}

cat2_outputs_df <- function() {
  cat2_outputs_ready()[1:3, , drop = FALSE]
}
