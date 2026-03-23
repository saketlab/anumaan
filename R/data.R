# data.R
# Data objects and mappings for AMR preprocessing

#' Default Column Name Mappings
#'
#' Named list mapping standard column names to common aliases found in
#' AMR surveillance datasets from different sources.
#'
#' @format A named list with 15 elements
#' @export
default_column_mappings <- list(
  patient_id = c(
    "PatientInformation_id", "PatientInformation", "patient_ID",
    "Patient_ID", "PatientID", "Subject_ID", "MRN",
    "medical_record_number", "UHID", "patient_no", "Case_ID"
  ),
  gender = c(
    "Gender", "sex", "Sex", "patient_gender", "Patient_Gender",
    "gender_code", "sex_code"
  ),
  state = c(
    "State", "patient_state", "state_name", "region",
    "State_Name", "state_code"
  ),
  location = c(
    "Location", "city", "hospital_city", "hospital_location",
    "site", "City", "Hospital_Location"
  ),
  DOB = c(
    "date_of_birth", "DateOfBirth", "birth_date", "dob", "DOB",
    "Date_of_Birth", "BirthDate"
  ),
  Age = c("age", "patient_age", "Age_Years", "age_years", "AGE"),
  date_of_admission = c(
    "admission_date", "Date.of.admission",
    "Date_of_admission_in_hospital", "AdmissionDate",
    "hospital_admission_date", "admit_date",
    "Date.of.admission.in.hospital"
  ),
  date_of_culture = c(
    "date_of_event", "Date.of.event", "culture_date",
    "collection_date", "sample_date", "specimen_date",
    "event_date", "culture_collection_date",
    "Date_of_culture", "CultureDate"
  ),
  date_of_final_outcome = c(
    "Date.of.14.day.outcome", "outcome_date",
    "discharge_date", "death_date", "Date_of_outcome",
    "final_outcome_date", "DateOfOutcome",
    "date_of_discharge"
  ),
  final_outcome = c(
    "Final.outcome", "outcome", "patient_outcome",
    "status", "final_status", "discharge_status",
    "Outcome", "Status"
  ),
  organism_name = c(
    "Organism", "organism", "pathogen", "pathogen_name",
    "bacteria", "microorganism", "organism_identified",
    "OrganismName", "isolated_organism"
  ),
  antibiotic_name = c(
    "antibiotic", "drug", "drug_name", "antimicrobial",
    "antimicrobial_name", "antibiotic_tested",
    "drug_tested", "Antibiotic", "AntibioticName"
  ),
  antibiotic_value = c(
    "antibiotic_result", "susceptibility", "resistance",
    "result", "remarks", "antibiotic_remarks",
    "susceptibility_result", "Remarks", "Result",
    "susceptibility_status"
  ),
  specimen_type = c(
    "specimen", "sample_type", "sample", "Sample_type1_name",
    "source", "specimen_source", "sample_source",
    "culture_source", "SpecimenType", "SampleType"
  ),
  diagnosis = c(
    "Diagnosis", "diagnosis_1", "primary_diagnosis",
    "clinical_diagnosis", "Diagnosis_1", "ICD_code",
    "diagnosis_code"
  )
)


#' Beta-Lactam Class Hierarchy
#'
#' Ordered vector of beta-lactam classes for resistance class selection.
#' Order represents clinical hierarchy (most to least important).
#'
#' @return Character vector
#' @export
get_beta_lactam_hierarchy <- function() {
  c(
    "Carbapenems",
    "Fourth-generation-cephalosporins",
    "Third-generation-cephalosporins",
    "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal",
    "Beta-lactam/beta-lactamase-inhibitor",
    "Aminopenicillins",
    "Penicillins"
  )
}


#' Standard Age Bins (GBD Compatible)
#'
#' Returns age bin boundaries for stratification, compatible with GBD
#' (Global Burden of Disease) methodology.
#'
#' @param type Character. "GBD_standard" for 5-year bins, "pediatric" for
#'   child-focused bins, or "geriatric" for elderly-focused bins.
#'
#' @return Character vector of age bin labels
#' @export
get_age_bins <- function(type = "GBD_standard") {
  switch(type,
    GBD_standard = c(
      "<1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30",
      "30-35", "35-40", "40-45", "45-50", "50-55", "55-60",
      "60-65", "65-70", "70-75", "75-80", "80-85", "85+"
    ),
    pediatric = c(
      "0-0.08", "0.08-1", "1-2", "2-5", "5-10", "10-15", "15-18", "18+"
    ),
    geriatric = c(
      "0-50", "50-60", "60-65", "65-70", "70-75", "75-80", "80-85",
      "85-90", "90+"
    ),
    stop("Unknown age bin type. Use 'GBD_standard', 'pediatric', or 'geriatric'")
  )
}


#' Get Magiorakos MDR/XDR Thresholds
#'
#' Returns pathogen-specific MDR/XDR classification criteria based on
#' Magiorakos et al. 2012 (Clin Microbiol Infect).
#'
#' @return Data frame with MDR/XDR thresholds per organism group
#' @export
#' @references
#' Magiorakos AP, Srinivasan A, Carey RB, et al. Multidrug-resistant,
#' extensively drug-resistant and pandrug-resistant bacteria: an international
#' expert proposal for interim standard definitions for acquired resistance.
#' Clin Microbiol Infect. 2012;18(3):268-281.
get_magiorakos_thresholds <- function() {
  tibble::tribble(
    ~organism_group, ~mdr_threshold, ~xdr_threshold, ~total_categories,
    "Enterobacterales", 3, "all_but_2", 9,
    "Pseudomonas aeruginosa", 3, "all_but_2", 10,
    "Acinetobacter spp", 3, "all_but_2", 9,
    "Staphylococcus aureus", 3, "all_but_2", 9,
    "Enterococcus spp", 3, "all_but_2", 7,
    "Streptococcus pneumoniae", 3, "all_but_2", 5
  )
}


#' Get Antimicrobial Categories for MDR/XDR Classification
#'
#' Returns antimicrobial categories used in Magiorakos MDR/XDR definitions.
#' Categories are pathogen-specific.
#'
#' @param organism_group Character. Organism group name.
#' @return Character vector of antimicrobial categories for that organism
#' @export
get_antimicrobial_categories <- function(organism_group = "Enterobacterales") {
  categories <- list(
    Enterobacterales = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (3rd gen)",
      "Cephalosporins (4th gen)",
      "Fluoroquinolones",
      "Monobactams",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Tigecycline"
    ),
    "Pseudomonas aeruginosa" = c(
      "Aminoglycosides",
      "Antipseudomonal carbapenems",
      "Antipseudomonal cephalosporins",
      "Antipseudomonal fluoroquinolones",
      "Antipseudomonal penicillins + beta-lactamase inhibitors",
      "Monobactams",
      "Phosphonic acids",
      "Polymyxins",
      "Ceftazidime-avibactam",
      "Ceftolozane-tazobactam"
    ),
    "Acinetobacter spp" = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (extended-spectrum)",
      "Fluoroquinolones",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Sulbactam",
      "Tetracyclines",
      "Tigecycline"
    ),
    "Staphylococcus aureus" = c(
      "Aminoglycosides",
      "Fluoroquinolones",
      "Folate pathway inhibitors",
      "Fusidic acid",
      "Glycopeptides",
      "Linezolid",
      "Mupirocin",
      "Rifampicin",
      "Tetracyclines"
    )
  )

  if (organism_group %in% names(categories)) {
    return(categories[[organism_group]])
  } else {
    warning(sprintf("No category list for '%s', using Enterobacterales default", organism_group))
    return(categories$Enterobacterales)
  }
}


#' Get Organism Taxonomy Mapping
#'
#' Reads the organism taxonomy from inst/extdata/organisms.csv and returns
#' a data frame mapping organism names to organism groups.
#'
#' @return Data frame with columns organism_name and org_group
#' @keywords internal
get_organism_taxonomy <- function() {
  file_path <- find_extdata_file("organisms.csv")
  if (file_path == "") {
    warning("organisms.csv not found. Returning empty taxonomy.")
    return(data.frame(
      organism_name = character(),
      org_group = character(),
      stringsAsFactors = FALSE
    ))
  }
  taxonomy <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  if ("organism_group" %in% names(taxonomy) && !"org_group" %in% names(taxonomy)) {
    names(taxonomy)[names(taxonomy) == "organism_group"] <- "org_group"
  }
  taxonomy
}


#' Get RR Pathogen Mapping
#'
#' Returns a mapping from normalized organism names to RR pathogen categories
#' used in burden estimation. Reads from inst/extdata/organisms.csv.
#'
#' @return Data frame with columns organism_name and rr_pathogen
#' @keywords internal
get_rr_pathogen_map <- function() {
  file_path <- find_extdata_file("organisms.csv")
  if (file_path == "") {
    warning("organisms.csv not found. Returning empty RR pathogen map.")
    return(data.frame(
      organism_name = character(),
      rr_pathogen = character(),
      stringsAsFactors = FALSE
    ))
  }
  taxonomy <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  # Use organism_name as the rr_pathogen if no dedicated column exists
  if (!"rr_pathogen" %in% names(taxonomy)) {
    taxonomy$rr_pathogen <- taxonomy$organism_name
  }
  taxonomy[, c("organism_name", "rr_pathogen")]
}


#' Get Antibiotic Class to RR Drug Mapping
#'
#' Returns a mapping from WHO antibiotic classes to RR drug categories.
#' Reads from inst/extdata/WHO_aware_class.csv.
#'
#' @return Data frame with columns Class and rr_drug
#' @keywords internal
get_class_rr_map <- function() {
  file_path <- find_extdata_file("WHO_aware_class.csv")
  if (file_path == "") {
    warning("WHO_aware_class.csv not found. Returning empty class-RR map.")
    return(data.frame(
      Class = character(),
      rr_drug = character(),
      stringsAsFactors = FALSE
    ))
  }
  who <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  if (!"rr_drug" %in% names(who)) {
    who$rr_drug <- who$Class
  }
  unique(who[, c("Class", "rr_drug")])
}


#' Normalize String for Joining
#'
#' Lowercases, trims whitespace, and removes punctuation for fuzzy joins.
#'
#' @param x Character vector
#' @return Character vector of normalized strings
#' @keywords internal
normalize_join <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- gsub("\\s+", " ", x)
  x
}
