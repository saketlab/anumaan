# Preprocessing Workflow in anumaan

``` r
library(anumaan)
library(dplyr)
```

## Overview

This vignette focuses on the **preprocessing layer only**. It does not
cover DALY estimation, spatial workflows, or downstream visualization.

The examples are deliberately small and synthetic, but they mirror the
kinds of inputs covered by the package’s Category 1 and Category 2
preprocessing tests:

- schema alignment and column mapping
- date parsing and table validation
- organism, specimen, antibiotic, and AST cleaning
- age, length-of-stay, and HAI/CAI derivation
- facility type assignment
- event creation, deduplication, and readmission classification
- MDR/XDR classification
- contaminant and polymicrobial handling
- analysis-readiness filtering

## A Small Synthetic Raw Dataset

The main example starts with non-standard column names plus the kinds of
messy values that preprocessing is supposed to absorb.

``` r
raw_alias <- data.frame(
  PID = c("pt_001", "pt_001", "pt_002", "pt_003", "pt_003", "pt_004"),
  Date.of.admission = c("2021-01-01", "2021-01-01", "43831", "2021/01/10", "2021/01/10", "20211301"),
  CultureDate = c("2021-01-03", "2021-01-03", "43833", "2021/01/11", "2021/01/11", "2021-01-15"),
  Date.of.14.day.outcome = c("2021-01-08", "2021-01-08", "43840", "2021/01/18", "2021/01/18", "2021-01-25"),
  Final.outcome = c("alive", "alive", "expired", "Discharged Alive", "Died/Expired", "LAMA"),
  Organism = c(
    "E. coli",
    "E. coli",
    "MRSA(Methicillin resistant staphylococcus aureus)",
    "CONS (Coagulase Negative Staphylococci)",
    "Non fermenting Gram negative bacilli",
    "No growth in culture."
  ),
  Sample_type1_name = c(
    "blood culure",
    "blood culure",
    "Blood",
    "Brain abscess",
    "ETA",
    "Superficial Biopsy"
  ),
  Antibiotic = c("Amikacin", "Ciproflox", "Vancomycin", "Oxacillin", "Colistin", "UnknownDrugZZ"),
  Result = c("Resistant", "Susceptible", "S", "Intermediate", "1", ""),
  Gender = c("M", "M", "female", "WOMAN", "man", "unknown"),
  DOB = c("1988-05-01", "1988-05-01", "1975/03/12", "1990-03-20", "19881301", NA),
  infection_type = c("Hospital-Acquired", "Hospital-Acquired", "Community acquired", "", "Healthcare Associated", ""),
  stringsAsFactors = FALSE
)

raw_alias
#>      PID Date.of.admission CultureDate Date.of.14.day.outcome    Final.outcome
#> 1 pt_001        2021-01-01  2021-01-03             2021-01-08            alive
#> 2 pt_001        2021-01-01  2021-01-03             2021-01-08            alive
#> 3 pt_002             43831       43833                  43840          expired
#> 4 pt_003        2021/01/10  2021/01/11             2021/01/18 Discharged Alive
#> 5 pt_003        2021/01/10  2021/01/11             2021/01/18     Died/Expired
#> 6 pt_004          20211301  2021-01-15             2021-01-25             LAMA
#>                                            Organism  Sample_type1_name
#> 1                                           E. coli       blood culure
#> 2                                           E. coli       blood culure
#> 3 MRSA(Methicillin resistant staphylococcus aureus)              Blood
#> 4           CONS (Coagulase Negative Staphylococci)      Brain abscess
#> 5              Non fermenting Gram negative bacilli                ETA
#> 6                             No growth in culture. Superficial Biopsy
#>      Antibiotic       Result  Gender        DOB        infection_type
#> 1      Amikacin    Resistant       M 1988-05-01     Hospital-Acquired
#> 2     Ciproflox  Susceptible       M 1988-05-01     Hospital-Acquired
#> 3    Vancomycin            S  female 1975/03/12    Community acquired
#> 4     Oxacillin Intermediate   WOMAN 1990-03-20                      
#> 5      Colistin            1     man   19881301 Healthcare Associated
#> 6 UnknownDrugZZ              unknown       <NA>
```

## Schema Checks and Column Mapping

### Inspect what preprocessing is possible

``` r
prep_report_capabilities(raw_alias)
```

### Build a rename map explicitly

Use the staged API when you want the mapping to be visible and
reviewable.

``` r
column_map <- prep_build_column_map(
  raw_alias,
  column_map = c(
    patient_id = "PID",
    date_of_admission = "Date.of.admission",
    date_of_culture = "CultureDate",
    date_of_final_outcome = "Date.of.14.day.outcome",
    final_outcome = "Final.outcome",
    organism_name = "Organism",
    specimen_type = "Sample_type1_name",
    antibiotic_name = "Antibiotic",
    antibiotic_value = "Result",
    gender = "Gender",
    DOB = "DOB"
  )
)

column_map
#>               patient_id        date_of_admission          date_of_culture 
#>                    "PID"      "Date.of.admission"            "CultureDate" 
#>    date_of_final_outcome            final_outcome            organism_name 
#> "Date.of.14.day.outcome"          "Final.outcome"               "Organism" 
#>            specimen_type          antibiotic_name         antibiotic_value 
#>      "Sample_type1_name"             "Antibiotic"                 "Result" 
#>                   gender                      DOB 
#>                 "Gender"                    "DOB"
```

``` r
prepped <- prep_apply_column_map(raw_alias, column_map)
names(prepped)
#>  [1] "patient_id"            "date_of_admission"     "date_of_culture"      
#>  [4] "date_of_final_outcome" "final_outcome"         "organism_name"        
#>  [7] "specimen_type"         "antibiotic_name"       "antibiotic_value"     
#> [10] "gender"                "DOB"                   "infection_type"
```

### One-step convenience wrapper

If you do not need the staged rename,
[`prep_standardize_column_names()`](https://saketlab.github.io/anumaan/reference/prep_standardize_column_names.md)
wraps the same logic in one call.

``` r
std_cols <- prep_standardize_column_names(raw_alias, fuzzy_match = TRUE)
names(std_cols$data)
#>  [1] "PID"                   "date_of_admission"     "date_of_culture"      
#>  [4] "date_of_final_outcome" "final_outcome"         "organism_name"        
#>  [7] "specimen_type"         "antibiotic_name"       "antibiotic_value"     
#> [10] "gender"                "DOB"                   "infection_type"
```

### Assert the minimum standard names

``` r
prep_assert_standard_names(
  prepped,
  required_standard_names = c(
    "patient_id",
    "date_of_culture",
    "organism_name",
    "antibiotic_name",
    "antibiotic_value"
  ),
  strict = TRUE
)
```

## Dates and Table Validation

### Parse a single vector of mixed date encodings

``` r
prep_parse_date_column(
  c("43831", "1577836800000", "20211301", "2021/01/05", "bad-date"),
  col_name = "date_of_culture",
  table_label = "example_dates"
)
#> [1] "2020-01-01" "2020-01-01" "2021-01-13" "2021-01-05" NA
```

### Validate and coerce the working table

[`prep_validate_table()`](https://saketlab.github.io/anumaan/reference/prep_validate_table.md)
is a good checkpoint before running downstream logic.

``` r
validated <- prep_validate_table(
  prepped,
  required_cols = c("patient_id", "organism_name", "antibiotic_name", "antibiotic_value"),
  key_col = "patient_id",
  date_cols = c("date_of_admission", "date_of_culture", "date_of_final_outcome", "DOB"),
  table_label = "raw_alias",
  stop_on_missing = FALSE
)

prepped <- validated$data
prepped[, c("patient_id", "date_of_admission", "date_of_culture", "date_of_final_outcome", "DOB")]
#>   patient_id date_of_admission date_of_culture date_of_final_outcome        DOB
#> 1     pt_001        2021-01-01      2021-01-03            2021-01-08 1988-05-01
#> 2     pt_001        2021-01-01      2021-01-03            2021-01-08 1988-05-01
#> 3     pt_002        2020-01-01      2020-01-03            2020-01-10 1975-03-12
#> 4     pt_003        2021-01-10      2021-01-11            2021-01-18 1990-03-20
#> 5     pt_003        2021-01-10      2021-01-11            2021-01-18 1988-01-13
#> 6     pt_004        2021-01-13      2021-01-15            2021-01-25       <NA>
```

### Validate chronology

``` r
prep_validate_date_logic(
  prepped,
  admission_col = "date_of_admission",
  culture_col = "date_of_culture",
  outcome_col = "date_of_final_outcome",
  dob_col = "DOB"
)
```

## Standardize Core Clinical Fields

### Organisms

``` r
prepped <- prep_standardize_organisms(prepped, organism_col = "organism_name")
prepped <- prep_assign_organism_group(prepped, organism_col = "organism_normalized")
#> 
#> Other 
#>     6
prepped <- prep_extract_genus(prepped, organism_col = "organism_normalized")
prepped <- prep_extract_species(prepped, organism_col = "organism_normalized")
prepped <- prep_flag_organism_unmatched(prepped, organism_col = "organism_normalized")

prepped[, c(
  "organism_name",
  "organism_normalized",
  "organism_group",
  "org_genus",
  "org_species",
  "is_organism_unmatched"
)]
#>                                       organism_name
#> 1                                           E. coli
#> 2                                           E. coli
#> 3 MRSA(Methicillin resistant staphylococcus aureus)
#> 4           CONS (Coagulase Negative Staphylococci)
#> 5              Non fermenting Gram negative bacilli
#> 6                             No growth in culture.
#>                    organism_normalized        organism_group      org_genus
#> 1                     escherichia coli      Enterobacterales    escherichia
#> 2                     escherichia coli      Enterobacterales    escherichia
#> 3                staphylococcus aureus   Gram-positive cocci staphylococcus
#> 4     coagulase-negative staphylococci   Gram-positive cocci      coagulase
#> 5 non-fermenting gram-negative bacilli Gram-negative bacilli            non
#> 6                                 <NA>                  <NA>           <NA>
#>     org_species is_organism_unmatched
#> 1          coli                 FALSE
#> 2          coli                 FALSE
#> 3        aureus                 FALSE
#> 4 staphylococci                 FALSE
#> 5          gram                 FALSE
#> 6          <NA>                 FALSE
```

### Specimen, sex, outcome, and infection type

In version 0.1.0.9021, specimen standardisation keeps the detailed
`sample_category` from the CSV reference but collapses sterile-site
labels to three downstream analysis classes: `Sterile`, `Non-Sterile`,
and `Others/Ambiguous`. The reference also includes stewardship labels
that appear in 2021-era extracts, including `Brain abscess`,
`Instrument`, `Lung aspirate`, `Lymph node`, `Hair`, and
`Superficial Biopsy`.

``` r
prepped <- prep_standardize_specimens(prepped, specimen_col = "specimen_type")
prepped <- prep_standardize_sex(prepped, col = "gender")
prepped <- prep_standardize_final_outcome(prepped, col = "final_outcome")
#> 
#>     Died Survived     <NA> 
#>        2        3        1
prepped <- prep_standardize_infection_type(prepped, col = "infection_type")
#> 
#>  CAI  HAI <NA> 
#>    1    3    2

prepped[, c(
  "specimen_type",
  "specimen_normalized",
  "sample_category",
  "sterile_classification",
  "gender",
  "outcome_std",
  "infection_type"
)]
#>        specimen_type         specimen_normalized   sample_category
#> 1       blood culure                       Blood             Blood
#> 2       blood culure                       Blood             Blood
#> 3              Blood                       Blood             Blood
#> 4      Brain abscess               Brain abscess        CNS/Ocular
#> 5                ETA Endotracheal aspirate (ETA) Respiratory tract
#> 6 Superficial Biopsy          Superficial Biopsy            Others
#>   sterile_classification gender outcome_std infection_type
#> 1                Sterile      M    Survived            HAI
#> 2                Sterile      M    Survived            HAI
#> 3                Sterile      F        Died            CAI
#> 4                Sterile      F    Survived           <NA>
#> 5            Non-Sterile      M        Died            HAI
#> 6       Others/Ambiguous   <NA>        <NA>           <NA>
```

``` r
specimen_9021 <- data.frame(
  specimen_type = c(
    "Instrument", "Lung aspirate", "Brain abscess",
    "Superficial Biopsy", "Hair", "Lymph node",
    "Drain Fluid", "Pigtail fluid"
  ),
  stringsAsFactors = FALSE
)

prep_standardize_specimens(specimen_9021)[, c(
  "specimen_type",
  "specimen_normalized",
  "sample_category",
  "sterile_classification"
)]
#>        specimen_type specimen_normalized   sample_category
#> 1         Instrument          Instrument   Catheter/Device
#> 2      Lung aspirate       Lung aspirate Respiratory tract
#> 3      Brain abscess       Brain abscess        CNS/Ocular
#> 4 Superficial Biopsy  Superficial Biopsy            Others
#> 5               Hair                Hair  Skin/Soft tissue
#> 6         Lymph node          Lymph node    Sterile tissue
#> 7        Drain Fluid         Drain Fluid        Body fluid
#> 8      Pigtail fluid       Pigtail fluid        Body fluid
#>   sterile_classification
#> 1       Others/Ambiguous
#> 2                Sterile
#> 3                Sterile
#> 4       Others/Ambiguous
#> 5            Non-Sterile
#> 6                Sterile
#> 7            Non-Sterile
#> 8       Others/Ambiguous
```

### Antibiotics and AST values

``` r
prepped <- prep_standardize_antibiotics(prepped, antibiotic_col = "antibiotic_name")
prepped <- prep_clean_ast_values(prepped, value_col = "antibiotic_value")
#> 
#>    I    R    S <NA> 
#>    1    1    2    2
prepped <- prep_harmonize_ast(
  prepped,
  antibiotic_col = "antibiotic_normalized",
  ast_col = "antibiotic_value"
)
#> 
#>    I    R    S <NA> 
#>    1    1    2    2
prepped <- prep_check_organism_ast_consistency(
  prepped,
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized"
)
prepped <- prep_flag_invalid_ast(prepped, col = "ast_value_harmonized")

prepped[, c(
  "antibiotic_name",
  "antibiotic_normalized",
  "antibiotic_class",
  "aware_category",
  "antibiotic_value",
  "ast_value_harmonized",
  "is_ast_inconsistent",
  "is_ast_invalid"
)]
#>   antibiotic_name antibiotic_normalized antibiotic_class aware_category
#> 1        Amikacin              amikacin  Aminoglycosides         Access
#> 2       Ciproflox         ciprofloxacin Fluoroquinolones          Watch
#> 3      Vancomycin         vancomycin_iv       Vancomycin          Watch
#> 4       Oxacillin             oxacillin      Penicillins         Access
#> 5        Colistin           colistin_iv         Colistin        Reserve
#> 6   UnknownDrugZZ         unknowndrugzz             <NA>           <NA>
#>   antibiotic_value ast_value_harmonized is_ast_inconsistent is_ast_invalid
#> 1                R                    R               FALSE          FALSE
#> 2                S                    S               FALSE          FALSE
#> 3                S                    S                TRUE          FALSE
#> 4                I                    R                TRUE          FALSE
#> 5             <NA>                 <NA>                TRUE          FALSE
#> 6             <NA>                 <NA>               FALSE          FALSE
```

## Demographic and Clinical Derivations

### Fill age and assign age bins

``` r
prepped <- prep_fill_age(
  prepped,
  age_col = "Age",
  dob_col = "DOB",
  date_col = "date_of_culture",
  overwrite = FALSE
)
#>            age_method age_confidence n
#> 1 calculated_from_dob           high 5
#> 2                <NA>           <NA> 1

prepped <- prep_assign_age_bins(prepped, age_col = "Age", bins = "GBD_standard")

prepped[, c("patient_id", "DOB", "date_of_culture", "Age", "age_method", "age_confidence", "Age_bin")]
#>   patient_id        DOB date_of_culture      Age          age_method
#> 1     pt_001 1988-05-01      2021-01-03 32.67625 calculated_from_dob
#> 2     pt_001 1988-05-01      2021-01-03 32.67625 calculated_from_dob
#> 3     pt_002 1975-03-12      2020-01-03 44.81314 calculated_from_dob
#> 4     pt_003 1990-03-20      2021-01-11 30.81451 calculated_from_dob
#> 5     pt_003 1988-01-13      2021-01-11 32.99658 calculated_from_dob
#> 6     pt_004       <NA>      2021-01-15       NA                <NA>
#>   age_confidence Age_bin
#> 1           high   30-35
#> 2           high   30-35
#> 3           high   40-45
#> 4           high   30-35
#> 5           high   30-35
#> 6           <NA>    <NA>
```

### Derive length of stay

``` r
prepped <- prep_derive_los_from_dates(
  prepped,
  admission_col = "date_of_admission",
  outcome_col = "date_of_final_outcome",
  los_col = "los_days"
)
#>   mean_los median_los min_los max_los
#> 1      8.5          8       7      12

prepped[, c("patient_id", "date_of_admission", "date_of_final_outcome", "los_days")]
#>   patient_id date_of_admission date_of_final_outcome los_days
#> 1     pt_001        2021-01-01            2021-01-08        7
#> 2     pt_001        2021-01-01            2021-01-08        7
#> 3     pt_002        2020-01-01            2020-01-10        9
#> 4     pt_003        2021-01-10            2021-01-18        8
#> 5     pt_003        2021-01-10            2021-01-18        8
#> 6     pt_004        2021-01-13            2021-01-25       12
```

### Derive HAI/CAI and source flags

``` r
prepped <- prep_derive_hai_cai(
  prepped,
  infection_type_col = "infection_type",
  admission_col = "date_of_admission",
  culture_col = "date_of_culture",
  hai_cutoff = 2,
  overwrite = FALSE
)
#>   infection_type infection_type_method n
#> 1            HAI  inferred_2day_cutoff 4
#> 2            CAI  inferred_2day_cutoff 2

prepped <- prep_flag_hai_inferred(prepped)
#> 
#> inferred 
#>        6

prepped[, c("patient_id", "date_of_admission", "date_of_culture", "infection_type", "infection_type_method", "infection_type_src")]
#>   patient_id date_of_admission date_of_culture infection_type
#> 1     pt_001        2021-01-01      2021-01-03            HAI
#> 2     pt_001        2021-01-01      2021-01-03            HAI
#> 3     pt_002        2020-01-01      2020-01-03            CAI
#> 4     pt_003        2021-01-10      2021-01-11            CAI
#> 5     pt_003        2021-01-10      2021-01-11            HAI
#> 6     pt_004        2021-01-13      2021-01-15            HAI
#>   infection_type_method infection_type_src
#> 1  inferred_2day_cutoff           inferred
#> 2  inferred_2day_cutoff           inferred
#> 3  inferred_2day_cutoff           inferred
#> 4  inferred_2day_cutoff           inferred
#> 5  inferred_2day_cutoff           inferred
#> 6  inferred_2day_cutoff           inferred
```

### Assign facility type

[`prep_assign_facility_type()`](https://saketlab.github.io/anumaan/reference/prep_assign_facility_type.md)
is deliberately generic – it has no hardcoded facility names. You supply
a named mapping from facility value to the category you want
(e.g. Government/Private, Urban/Rural, or a region), and facilities not
found in the mapping get `default` (`NA` unless you set one).

``` r
facility_demo <- data.frame(
  center_name = c("General Hospital", "City Medical Center", "Unknown Hospital"),
  stringsAsFactors = FALSE
)

prep_assign_facility_type(
  facility_demo,
  mapping = c("General Hospital" = "Government", "City Medical Center" = "Private")
)
#>           center_name facility_type
#> 1    General Hospital    Government
#> 2 City Medical Center       Private
#> 3    Unknown Hospital          <NA>
```

## Event Creation and AST Reshaping

### Create event IDs from long-format isolates

``` r
prepped <- prep_create_event_ids(
  prepped,
  patient_col = "patient_id",
  date_col = "date_of_culture",
  organism_col = "organism_normalized",
  specimen_col = "specimen_normalized",
  antibiotic_col = "antibiotic_normalized",
  value_col = "ast_value_harmonized",
  gap_days = 14
)

prepped[, c("patient_id", "organism_normalized", "date_of_culture", "event_id")]
#>   patient_id                  organism_normalized date_of_culture
#> 1     pt_001                     escherichia coli      2021-01-03
#> 2     pt_001                     escherichia coli      2021-01-03
#> 3     pt_002                staphylococcus aureus      2020-01-03
#> 4     pt_003     coagulase-negative staphylococci      2021-01-11
#> 5     pt_003 non-fermenting gram-negative bacilli      2021-01-11
#> 6     pt_004                                 <NA>      2021-01-15
#>                                                                               event_id
#> 1                                           pt_001_blood_20210103_escherichia coli_001
#> 2                                           pt_001_blood_20210103_escherichia coli_001
#> 3                                      pt_002_blood_20200103_staphylococcus aureus_001
#> 4                   pt_003_brain abscess_20210111_coagulase-negative staphylococci_001
#> 5 pt_003_endotracheal aspirate (eta)_20210111_non-fermenting gram-negative bacilli_002
#> 6                                        pt_004_superficial biopsy_20210115___NA___001
```

### Deduplicate event-level rows

``` r
event_dedup <- prep_deduplicate_events(
  prepped,
  event_col = "event_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized"
)

nrow(prepped)
#> [1] 6
nrow(event_dedup)
#> [1] 6
```

### Classify readmissions from admission dates

[`prep_flag_readmission()`](https://saketlab.github.io/anumaan/reference/prep_flag_readmission.md)
classifies each admission relative to the patient’s previous one, based
on the gap in days:

- **`index`** – the patient’s first admission
- **`linked_readmission`** – within `gap_linked_days` (default 30) of
  the previous admission; treated as the same clinical episode
- **`new_readmission`** – between `gap_linked_days` and `gap_new_days`
  (default 90)
- **`late_readmission`** – more than `gap_new_days` later; a fully new
  event

``` r
readmit_demo <- data.frame(
  patient_id     = rep("pt_101", 4),
  admission_date = as.Date(c("2021-01-01", "2021-01-10", "2021-03-05", "2021-09-20")),
  stringsAsFactors = FALSE
)

prep_flag_readmission(readmit_demo, patient_col = "patient_id", admission_col = "admission_date")
#> 
#>              index   late_readmission linked_readmission    new_readmission 
#>                  1                  1                  1                  1
#> # A tibble: 4 × 3
#>   patient_id admission_date readmission_class 
#>   <chr>      <date>         <chr>             
#> 1 pt_101     2021-01-01     index             
#> 2 pt_101     2021-01-10     linked_readmission
#> 3 pt_101     2021-03-05     new_readmission   
#> 4 pt_101     2021-09-20     late_readmission
```

### Standardize non-conforming readmission labels

If a `readmission_class` column already exists but uses non-standard
labels (e.g. from an external system),
[`prep_classify_readmission()`](https://saketlab.github.io/anumaan/reference/prep_classify_readmission.md)
remaps common synonyms to the controlled vocabulary above. Unrecognised
labels become `NA` rather than being guessed at.

``` r
messy_readmit <- data.frame(
  patient_id        = paste0("pt_", 1:5),
  readmission_class = c("First", "READMIT", "Same_Episode", "late", "unmapped_value"),
  stringsAsFactors = FALSE
)

prep_classify_readmission(messy_readmit)
#> 
#>              index   late_readmission linked_readmission    new_readmission 
#>                  1                  1                  1                  1 
#>               <NA> 
#>                  1
#>   patient_id  readmission_class
#> 1       pt_1              index
#> 2       pt_2    new_readmission
#> 3       pt_3 linked_readmission
#> 4       pt_4   late_readmission
#> 5       pt_5               <NA>
```

### Pivot wide AST input to long format

``` r
wide_ast <- data.frame(
  patient_id = c("pt_010", "pt_011"),
  organism_name = c("Escherichia coli", "Klebsiella pneumoniae"),
  AMK = c("S", "R"),
  CIP = c("R", "S"),
  CTX = c("I", NA),
  stringsAsFactors = FALSE
)

ast_long <- prep_pivot_ast_wide_to_long(
  wide_ast,
  id_cols = c("patient_id", "organism_name"),
  antibiotic_cols = c("AMK", "CIP", "CTX"),
  remove_missing = TRUE
)

ast_long
#> # A tibble: 5 × 4
#>   patient_id organism_name         antibiotic_name antibiotic_value
#>   <chr>      <chr>                 <chr>           <chr>           
#> 1 pt_010     Escherichia coli      AMK             S               
#> 2 pt_010     Escherichia coli      CIP             R               
#> 3 pt_010     Escherichia coli      CTX             I               
#> 4 pt_011     Klebsiella pneumoniae AMK             R               
#> 5 pt_011     Klebsiella pneumoniae CIP             S
```

### Deduplicate AST records

[`prep_deduplicate_ast()`](https://saketlab.github.io/anumaan/reference/prep_deduplicate_ast.md)
runs in two sequential steps:

1.  **Exact-duplicate removal** — drops rows that are fully identical
    across patient, organism, antibiotic, date, *and* value. These are
    true redundant records (e.g. the same result entered twice) with no
    ambiguity.
2.  **Conflict resolution** — after exact duplicates are gone, detects
    groups where the same patient + organism + antibiotic + date has
    *different* values (e.g. one lab returns S, another returns R), then
    flags or resolves them via a strategy.

**Step 1 — exact duplicates removed automatically:**

``` r
ast_exact_dups <- data.frame(
  patient_id            = c("pt_030", "pt_030", "pt_031"),
  culture_date          = as.Date(c("2021-02-01", "2021-02-01", "2021-02-01")),
  organism_normalized   = c("escherichia coli", "escherichia coli", "klebsiella pneumoniae"),
  antibiotic_normalized = c("amikacin", "amikacin", "ciprofloxacin"),
  ast_value_harmonized  = c("R", "R", "S"),
  stringsAsFactors      = FALSE
)

# 3 rows in — rows 1 and 2 are identical on all 5 key columns
# step 1 removes the duplicate → 2 rows out, no conflicts flagged
out_exact <- prep_deduplicate_ast(ast_exact_dups, mode = "detect")
nrow(out_exact)
#> [1] 2
out_exact[, c("patient_id", "antibiotic_normalized", "ast_value_harmonized", "is_ast_duplicate")]
#>   patient_id antibiotic_normalized ast_value_harmonized is_ast_duplicate
#> 1     pt_030              amikacin                    R            FALSE
#> 2     pt_031         ciprofloxacin                    S            FALSE
```

**Step 2 — conflict detection after exact-dedup:**

``` r
ast_dups <- data.frame(
  patient_id            = c("pt_020", "pt_020", "pt_020"),
  culture_date          = as.Date(c("2021-01-15", "2021-01-15", "2021-01-15")),
  organism_normalized   = c("escherichia coli", "escherichia coli", "escherichia coli"),
  antibiotic_normalized = c("amikacin", "amikacin", "cefotaxime"),
  antibiotic_name       = c("Amikacin", "Amikacin", "Cefotaxime"),
  antibiotic_class      = c("Aminoglycosides", "Aminoglycosides", "Third-generation-cephalosporins"),
  antibiotic_value      = c("S", "R", "I"),
  ast_value_harmonized  = c("S", "R", "I"),
  stringsAsFactors      = FALSE
)

# detect mode: inspect which rows are flagged as conflicting
prep_deduplicate_ast(ast_dups, mode = "detect")[, c(
  "patient_id",
  "antibiotic_normalized",
  "ast_value_harmonized",
  "is_ast_duplicate"
)]
#> # A tibble: 1 × 6
#>   patient_id organism_normalized antibiotic_normalized culture_date
#>   <chr>      <chr>               <chr>                 <date>      
#> 1 pt_020     escherichia coli    amikacin              2021-01-15  
#> # ℹ 2 more variables: conflicting_values <chr>, n_rows <int>
#>   patient_id antibiotic_normalized ast_value_harmonized is_ast_duplicate
#> 1     pt_020              amikacin                    S             TRUE
#> 2     pt_020              amikacin                    R             TRUE
#> 3     pt_020            cefotaxime                    I            FALSE

# remove mode: resolve conflicts in one call — resistant_wins keeps "R"
prep_deduplicate_ast(ast_dups, mode = "remove", strategy = "resistant_wins")[, c(
  "patient_id",
  "antibiotic_normalized",
  "ast_value_harmonized"
)]
#> # A tibble: 1 × 6
#>   patient_id organism_normalized antibiotic_normalized culture_date
#>   <chr>      <chr>               <chr>                 <date>      
#> 1 pt_020     escherichia coli    amikacin              2021-01-15  
#> # ℹ 2 more variables: conflicting_values <chr>, n_rows <int>
#> # A tibble: 2 × 3
#>   patient_id antibiotic_normalized ast_value_harmonized
#>   <chr>      <chr>                 <chr>               
#> 1 pt_020     amikacin              R                   
#> 2 pt_020     cefotaxime            I
```

## MDR/XDR Classification

MDR (Multidrug Resistant) and XDR (Extensively Drug Resistant)
classification follows the Magiorakos 2012 criteria: MDR means resistant
to at least one agent in three or more antimicrobial categories; XDR
means susceptible to at most two categories. Both are computed per
`event_id`, so this must run *after* event creation.

### Collapse to class level

[`prep_classify_mdr_xdr()`](https://saketlab.github.io/anumaan/reference/prep_classify_mdr_xdr.md)
requires a `class_result_event` column, which comes from
[`prep_collapse_class_level()`](https://saketlab.github.io/anumaan/reference/prep_collapse_class_level.md)
(renaming its output `class_resistance` column). Pass `extra_cols` to
carry `organism_group` through the aggregation – it is dropped
otherwise, since the function only groups by
`event_col`/`organism_col`/`class_col`.

``` r
class_level <- prep_collapse_class_level(
  prepped,
  event_col          = "event_id",
  organism_col       = "organism_normalized",
  class_col          = "antibiotic_class",
  susceptibility_col = "antibiotic_value",
  extra_cols         = "organism_group"
)
#> # A tibble: 6 × 5
#>   antibiotic_class     R  `NA`     S     I
#>   <chr>            <int> <int> <int> <int>
#> 1 Aminoglycosides      1     0     0     0
#> 2 Colistin             0     1     0     0
#> 3 Fluoroquinolones     0     0     1     0
#> 4 Penicillins          0     0     0     1
#> 5 Vancomycin           0     0     1     0
#> 6 NA                   0     1     0     0

class_level <- dplyr::rename(class_level, class_result_event = class_resistance)
```

### Classify MDR/XDR

``` r
mdr_xdr <- prep_classify_mdr_xdr(class_level, organism_group_col = "organism_group")

dplyr::distinct(mdr_xdr, event_id, mdr, mdr_confidence, xdr, xdr_confidence)
#> # A tibble: 5 × 5
#>   event_id                             mdr   mdr_confidence xdr   xdr_confidence
#>   <chr>                                <lgl> <chr>          <lgl> <chr>         
#> 1 pt_001_blood_20210103_escherichia c… FALSE insufficient_… TRUE  high          
#> 2 pt_002_blood_20200103_staphylococcu… FALSE NA             TRUE  medium        
#> 3 pt_003_brain abscess_20210111_coagu… FALSE NA             TRUE  medium        
#> 4 pt_003_endotracheal aspirate (eta)_… FALSE NA             TRUE  medium        
#> 5 pt_004_superficial biopsy_20210115_… FALSE NA             TRUE  medium
```

With only 6 rows in this synthetic example, most events have too few
antibiotic classes tested for a confident call – `mdr_confidence`
reports `"insufficient_data"` rather than a false-confident answer. On a
real dataset with a fuller antibiotic panel per event, confidence rises
to `"medium"`/`"high"` as more categories are tested.

## Contaminants and Polymicrobial Episodes

### Flag likely contaminants

``` r
prepped <- prep_flag_contaminants(prepped, method = "auto")

prepped[, c("organism_normalized", "specimen_normalized", "is_contaminant", "contaminant_confidence", "contaminant_method")]
#>                    organism_normalized         specimen_normalized
#> 1                     escherichia coli                       Blood
#> 2                     escherichia coli                       Blood
#> 3                staphylococcus aureus                       Blood
#> 4     coagulase-negative staphylococci               Brain abscess
#> 5 non-fermenting gram-negative bacilli Endotracheal aspirate (ETA)
#> 6                                 <NA>          Superficial Biopsy
#>   is_contaminant contaminant_confidence contaminant_method
#> 1          FALSE                    low          heuristic
#> 2          FALSE                    low          heuristic
#> 3          FALSE                    low          heuristic
#> 4          FALSE                    low          heuristic
#> 5          FALSE                    low          heuristic
#> 6          FALSE                    low          heuristic
```

### Exclude fungal isolates when needed

``` r
prep_filter_fungal(prepped, group_col = "organism_group")[, c("patient_id", "organism_normalized", "organism_group")]
#>   patient_id                  organism_normalized        organism_group
#> 1     pt_001                     escherichia coli      Enterobacterales
#> 2     pt_001                     escherichia coli      Enterobacterales
#> 3     pt_002                staphylococcus aureus   Gram-positive cocci
#> 4     pt_003     coagulase-negative staphylococci   Gram-positive cocci
#> 5     pt_003 non-fermenting gram-negative bacilli Gram-negative bacilli
#> 6     pt_004                                 <NA>                  <NA>
```

### Flag and weight polymicrobial episodes

[`prep_flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/prep_flag_polymicrobial.md)
is easiest to interpret at the event level, where each row is already
attached to an episode identifier.

``` r
poly_input <- data.frame(
  patient_id = c("pt_poly_1", "pt_poly_1", "pt_poly_2"),
  event_id = c("ev_poly_1", "ev_poly_1", "ev_poly_2"),
  organism_normalized = c("escherichia coli", "klebsiella pneumoniae", "escherichia coli"),
  stringsAsFactors = FALSE
)

poly_flagged <- prep_flag_polymicrobial(
  poly_input,
  patient_col = "event_id",
  organism_col = "organism_normalized"
)
#> # A tibble: 2 × 2
#>   n_organisms n_groups
#>         <int>    <int>
#> 1           1        1
#> 2           2        1

poly_weighted <- prep_compute_poly_weights(
  poly_flagged,
  episode_col = "event_id",
  organism_col = "organism_normalized",
  polymicrobial_col = "is_polymicrobial",
  method = "equal"
)
#> # A tibble: 1 × 4
#>   mean_weight median_weight min_weight max_weight
#>         <dbl>         <dbl>      <dbl>      <dbl>
#> 1         0.5           0.5        0.5        0.5

poly_split <- prep_split_poly_episode(poly_weighted, strategy = "fractional")

poly_split
#> # A tibble: 3 × 7
#>   patient_id organism_normalized   n_organisms is_polymicrobial poly_weight
#>   <chr>      <chr>                       <int>            <int>       <dbl>
#> 1 pt_poly_1  escherichia coli                2                1         0.5
#> 2 pt_poly_1  klebsiella pneumoniae           2                1         0.5
#> 3 pt_poly_2  escherichia coli                1                0         1  
#> # ℹ 2 more variables: weight_method <chr>, weight_confidence <chr>
```

## Attrition and Analysis Readiness

### Track attrition across preprocessing stages

``` r
flow <- NULL
flow <- prep_attrition_flow(flow, raw_alias, "raw_input", "Unmapped raw extract", patient_col = "PID")
flow <- prep_attrition_flow(flow, prepped, "preprocessed", "After schema, cleaning, and derivation", patient_col = "patient_id", event_col = "event_id")

ready <- prep_filter_analysis_ready(
  prepped,
  patient_col = "patient_id",
  culture_date_col = "date_of_culture",
  organism_col = "organism_name",
  antibiotic_col = "antibiotic_name",
  ast_col = "ast_value_harmonized",
  contaminant_col = "is_contaminant"
)

flow <- prep_attrition_flow(flow, ready, "analysis_ready", "Rows retained for downstream analysis", patient_col = "patient_id", event_col = "event_id")
flow
#>            stage n_rows n_patients n_events n_removed
#> 1      raw_input      6          4       NA         0
#> 2   preprocessed      6          4        5         0
#> 3 analysis_ready      4          3        3         2
#>                                   reason
#> 1                   Unmapped raw extract
#> 2 After schema, cleaning, and derivation
#> 3  Rows retained for downstream analysis
```

### Final readiness checks

``` r
prep_missingness_report(
  ready,
  threshold = 20,
  cols = c("patient_id", "organism_name", "antibiotic_name", "ast_value_harmonized", "Age")
)
#>               col_name n_total n_missing pct_missing is_high_missing
#> 1           patient_id       4         0           0           FALSE
#> 2        organism_name       4         0           0           FALSE
#> 3      antibiotic_name       4         0           0           FALSE
#> 4 ast_value_harmonized       4         0           0           FALSE
#> 5                  Age       4         0           0           FALSE

prep_validate_analysis_ready(ready, min_rows = 1, stop_on_failure = FALSE)
```

## Summary

The preprocessing layer is designed to be used in two ways:

- **modular**, with explicit staged calls like the ones above
- **pipeline-first**, by wrapping the same steps inside
  [`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)

For development, the modular route is the more useful one because it
makes every transformation inspectable and keeps failures local to one
stage.

Once this workflow is stable for your dataset, the next step is usually
to:

1.  freeze the column map for that source
2.  run the staged workflow on a larger extract
3.  compare attrition before and after each filter
4.  wrap the validated sequence inside
    [`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)
