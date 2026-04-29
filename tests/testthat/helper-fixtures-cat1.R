# Category 1 shared fixtures

cat1_base_data <- function() {
  # 4 patients with isolate counts: p1=3, p2=5, p3=1, p4=6
  patient_id <- c(
    rep("p1", 3),
    rep("p2", 5),
    "p3",
    rep("p4", 6)
  )

  isolate_id <- c(
    sprintf("p1_iso_%02d", 1:3),
    sprintf("p2_iso_%02d", 1:5),
    "p3_iso_01",
    sprintf("p4_iso_%02d", 1:6)
  )

  isolate_df <- data.frame(
    patient_id = patient_id,
    isolate_id = isolate_id,
    stringsAsFactors = FALSE
  )

  isolate_df$admission_date <- as.Date("2020-01-01") + seq_len(nrow(isolate_df)) - 1
  isolate_df$culture_date <- isolate_df$admission_date + 1
  isolate_df$outcome_date <- isolate_df$culture_date + 5

  isolate_df$dob <- as.Date(c(
    rep("1988-05-01", 3),
    rep("1975-03-12", 5),
    "2019-08-15",
    rep("1968-10-22", 6)
  ))
  isolate_df$age <- c(
    rep(31, 3),
    rep(44, 5),
    0.2,
    rep(51, 6)
  )
  # Keep explicit age component columns for age-binning scenarios.
  isolate_df$age_years <- isolate_df$age
  isolate_df$age_months <- c(
    rep(0, 8),
    2,  # 0.2 years-ish pediatric example when interpreted as months/day components
    rep(0, 6)
  )
  isolate_df$age_days <- c(
    rep(0, 8),
    15,
    rep(0, 6)
  )
  # Aliases requested by user: year/months columns.
  isolate_df$year <- isolate_df$age_years
  isolate_df$months <- isolate_df$age_months
  isolate_df$final_outcome <- c(
    rep("Survived", 3),
    rep("Died", 5),
    "Survived",
    rep("Survived", 6)
  )
  isolate_df$infection_type <- c(
    rep("HAI", 3),
    rep("CAI", 5),
    "CAI",
    rep("HAI", 6)
  )
  isolate_df$specimen_type <- c(
    rep("Blood culture", 3),
    rep("Urine culture", 5),
    "CSF",
    rep("Blood culture", 6)
  )
  # Raw organism names intentionally include spelling/abbreviation variants.
  isolate_df$organism_name_raw <- c(
    rep("E.coli", 3),
    c("Kleb pneumo", "Staph Aureus", "Acinetabecter sp..", "CONS", "K. pneumoniae"),
    "E. coli",
    c("P. aeruginosa", "Escherichia  coli", "staph. aureus", "Enterococcus spp", "Acineto bacter sp", "Unknown bug")
  )
  isolate_df$organism_normalized <- c(
    rep("Escherichia coli", 3),
    c("Klebsiella pneumoniae", "Staphylococcus aureus", "Acinetobacter spp.", "Staphylococcus epidermidis", "Klebsiella pneumoniae"),
    "Escherichia coli",
    c("Pseudomonas aeruginosa", "Escherichia coli", "Staphylococcus aureus", "Enterococcus spp.", "Acinetobacter spp.", NA_character_)
  )

  # Patient 1 requested pattern: 3 isolates x 4 antibiotics = 12 rows (long format)
  p1_iso <- isolate_df[isolate_df$patient_id == "p1", , drop = FALSE]
  p1_long <- p1_iso[rep(seq_len(nrow(p1_iso)), each = 4), , drop = FALSE]
  p1_long$antibiotic_name_raw <- rep(c("AMK", "WrongDrugZZ", "CIP", "Amekacin"), times = nrow(p1_iso))
  p1_long$antibiotic_value <- rep(c("I", "S", "R", "Resistant"), times = nrow(p1_iso))

  other_iso <- isolate_df[isolate_df$patient_id != "p1", , drop = FALSE]
  other_long <- other_iso
  other_long$antibiotic_name_raw <- c(
    "Colistin",
    "Nalidixic acid",
    "Ciproflox",
    "Amikacin",
    "AMP",
    "TZP",
    NA_character_,
    "NULL",
    "UnknownAbx",
    "CIP",
    "AMK",
    "Nalidixic acid"
  )
  ast_levels <- c(
    "Susceptible", "Intermediate", "No conclusion",
    NA_character_, "Unknown", "Not known", "NULL",
    "S", "I", "R", "Resistant", "Susceptible"
  )
  other_long$antibiotic_value <- ast_levels[seq_len(nrow(other_iso))]

  out <- rbind(p1_long, other_long)

  abx_map <- c(
    "AMK" = "Amikacin",
    "CIP" = "Ciprofloxacin",
    "Ciproflox" = "Ciprofloxacin",
    "AMP" = "Ampicillin",
    "TZP" = "Piperacillin-Tazobactam",
    "Amikacin" = "Amikacin",
    "Amekacin" = "Amikacin",
    "Colistin" = "Colistin",
    "Nalidixic acid" = "Nalidixic acid"
  )
  out$antibiotic_name_std <- unname(abx_map[out$antibiotic_name_raw])
  out$antibiotic_normalized <- tolower(out$antibiotic_name_std)

  out
}

cat1_with_missing <- function() {
  out <- cat1_base_data()

  out$final_outcome[c(2, 7)] <- c(NA_character_, "")
  out$infection_type[c(3, 12)] <- NA_character_
  out$organism_normalized[11] <- NA_character_
  out$antibiotic_name_raw[c(5, 9)] <- c(NA_character_, "")
  out$antibiotic_name_std[c(5, 9)] <- NA_character_
  out$antibiotic_normalized[c(5, 9)] <- NA_character_

  out
}

cat1_with_invalid_dates <- function() {
  out <- cat1_base_data()

  out$admission_date <- format(out$admission_date, "%d-%m-%Y")
  out$culture_date <- as.character(out$culture_date)
  out$outcome_date <- format(out$outcome_date, "%Y/%m/%d")
  out$dob <- format(out$dob, "%Y%m%d")

  # Mix of encoded / serial / swapped-digit styles for parser tests.
  out$culture_date[1] <- "43831"                 # Excel serial
  out$culture_date[2] <- "1577923200000"         # Unix timestamp (ms)
  out$culture_date[3] <- "20201301"              # reversed YYYYDDMM
  out$culture_date[4] <- "323032302d30312d3034"  # encoded-like undecodable string

  out
}

cat1_edge_neonatal <- function() {
  out <- cat1_base_data()

  out$age_years <- floor(out$age)
  out$age_months <- 0
  out$age_days <- 0

  # First 3 rows represent neonatal/fine-grained ages.
  out$age[1:3] <- c(0.00, 0.04, 0.30)
  out$age_years[1:3] <- c(0, 0, 0)
  out$age_months[1:3] <- c(0, 0, 3)
  out$age_days[1:3] <- c(3, 15, 0)

  out
}

# Scenario A: only DOB is provided; age should be derived after date parsing/coercion.
cat1_age_from_dob_only <- function() {
  out <- cat1_base_data()

  out$age <- NA_real_
  out$age_years <- NA_real_
  out$year <- NA_real_
  out$age_months <- NA_real_
  out$months <- NA_real_
  out$age_days <- NA_real_

  out$dob <- as.character(out$dob)
  out$culture_date <- as.character(out$culture_date)
  out$admission_date <- as.character(out$admission_date)

  # Mixed DOB/date formats including coded/faulty values.
  out$dob[1:8] <- c(
    "19880501",         # compact ymd
    "12-03-1975",       # dmy
    "1988/05/01",       # slash format
    "19750312",         # compact ymd
    "19881301",         # reversed YYYYDDMM-style
    "313938382d30352d3031", # encoded-like invalid
    NA_character_,
    ""
  )
  out$culture_date[1:8] <- c(
    "43831",            # Excel serial
    "1577923200000",    # Unix ms
    "2020-01-03",       # ISO
    "202-10-13",        # malformed year format
    "20201301",         # reversed YYYYDDMM
    "323032302d30312d3034", # encoded-like invalid
    "N/A",
    "NULL"
  )
  out$admission_date[1:8] <- c(
    "2020-01-01",       # ISO
    "01-02-2020",       # mdy
    "20200103",         # compact ymd
    "43830",            # Excel serial
    "1577836800000",    # Unix ms
    "202-10-13",        # malformed year format
    "N/A",
    "NULL"
  )

  out
}

# Scenario B: age in years has invalid values (negative/too large/decimal).
cat1_age_years_edge <- function() {
  out <- cat1_base_data()

  edge_years <- c(
    -2, -0.1, 0, 0.2, 0.8, 1.5, 2.9, 5.2, 12.7, 17.9, 65.5, 90.3,
    121, 130, NA, 0.04, 0.30, 18, 44, 51, 0.02, 0.08, 0.25, 1.0
  )
  out$age <- edge_years
  out$age_years <- edge_years
  out$year <- edge_years

  out
}

# Scenario C: derive decimal age from years + months + days components.
cat1_age_from_components <- function() {
  out <- cat1_base_data()

  out$age <- NA_real_
  out$age_years <- c(0, 0, 0, 0, 1, 5, 10, 18, 30, 45, 60, 75, 0, 0, 2, 4, 8, 16, 25, 40, 50, 65, 80, NA)
  out$age_months <- c(0, 0, 1, 3, 6, 0, 0, 0, 0, 0, 0, 0, NA, 11, 0, 6, 0, 0, 0, 0, 0, 0, 0, 2)
  out$age_days <- c(3, 15, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 28, 0, 5, 0, 2, 0, 0, 0, 0, 0, 0, NA)
  out$year <- out$age_years
  out$months <- out$age_months

  out
}

cat1_min_age <- function() {
  data.frame(
    Age        = c(0.5, 5, 10, 40, 90, NA),
    stringsAsFactors = FALSE
  )
}

# Compound age columns: years + months + days (for neonatal / fine-grained testing)
cat1_compound_age <- function() {
  data.frame(
    age_years  = c(0,  0,  0,  1,  30, NA),
    age_months = c(0,  0,  3,  6,  0,  NA),
    age_days   = c(3,  15, 0,  0,  0,  0),
    stringsAsFactors = FALSE
  )
}

cat1_min_dates <- function() {
  data.frame(
    admission_date = as.Date(c("2020-01-01", "2020-01-03")),
    culture_date = as.Date(c("2020-01-02", "2020-01-04")),
    outcome_date = as.Date(c("2020-01-05", "2020-01-06")),
    dob = as.Date(c("1980-01-01", "1990-01-01")),
    age = c(40, 30),
    final_outcome = c("Survived", "Died"),
    stringsAsFactors = FALSE
  )
}

cat1_min_outcome <- function() {
  data.frame(
    final_outcome = c("alive", "expired", "LAMA", NA, ""),
    stringsAsFactors = FALSE
  )
}

cat1_min_infection <- function() {
  data.frame(
    infection_type = c("hospital acquired", "community acquired", "other", NA),
    stringsAsFactors = FALSE
  )
}
