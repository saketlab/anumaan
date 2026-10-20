# AMR Data Preprocessing Pipeline

Master function that orchestrates the complete AMR data preprocessing
pipeline. Applies standardization, enrichment, and derivation in
sequence.

## Usage

``` r
run_preprocess(
  data,
  config = NULL,
  phases = "all",
  verbose = TRUE,
  validate = TRUE,
  generate_report = TRUE
)
```

## Arguments

- data:

  Data frame. Raw AMR dataset.

- config:

  Configuration object from amr_config(). If NULL, uses defaults.

- phases:

  Character vector. Which phases to run: "standardize", "enrich",
  "derive", or "all" (default). Allows partial pipeline execution.

- verbose:

  Logical. Print detailed progress messages. Default TRUE.

- validate:

  Logical. Run validation checks before and after. Default TRUE.

- generate_report:

  Logical. Generate preprocessing report. Default TRUE.

## Value

List with class "amr_result": - data: Preprocessed data frame - config:
Configuration used - log: List of processing logs and summaries -
report: Preprocessing report (if generate_report = TRUE) - metadata:
Pipeline execution metadata

## Details

Pipeline phases: 1. \*\*Standardization\*\*: Column mapping, value
normalization, date parsing 2. \*\*Enrichment\*\*: Derive missing
optional variables (Age, LOS, infection type) 3. \*\*Derivation\*\*:
Create analytical variables (event IDs, MDR/XDR, weights)

## Examples

``` r
raw <- data.frame(
  PID = c("pt_001", "pt_001", "pt_002", "pt_003"),
  Date.of.admission = c("2021-01-01", "2021-01-01", "2021-01-05", "2021-01-10"),
  CultureDate = c("2021-01-03", "2021-01-03", "2021-01-08", "2021-01-11"),
  Final.outcome = c("alive", "alive", "expired", "Discharged Alive"),
  Organism = c("E. coli", "E. coli", "Klebsiella pneumoniae", "Staphylococcus aureus"),
  Sample_type1_name = c("Blood", "Blood", "Urine", "Blood"),
  Antibiotic = c("Amikacin", "Ciprofloxacin", "Meropenem", "Vancomycin"),
  Result = c("Resistant", "Susceptible", "R", "S"),
  Gender = c("M", "M", "female", "F"),
  Age = c(34, 34, 61, 47)
)
result <- run_preprocess(raw, verbose = FALSE, generate_report = FALSE)
#> [v] Quality check passed: 4 rows x 10 cols, 100.0% complete
#> Auto fuzzy: 'Age' -> 'Age' (dist 0.00)
#> [prep_apply_column_map] Renamed 9 column(s).
#> [prep_apply_column_map] 1 column(s) not in map (kept as-is): PID
#> [pipeline] 'date_of_admission': 4 / 4 non-missing value(s) successfully parsed as Date (0 failed).
#> [pipeline] 'date_of_culture': 4 / 4 non-missing value(s) successfully parsed as Date (0 failed).
#> [prep_validate_date_logic] All date logic checks passed.
#> Standardized gender: M=2, F=2, NA=0
#> Standardized final outcome distribution:
#> 
#>     Died Survived 
#>        1        3 
#> Cleaning antibiotic values: 4 unique values found
#> Cleaned: 4 unique values -> 2 (S/I/R)
#> 
#> Value distribution:
#> 
#> R S 
#> 2 2 
#> No intermediate (I) values found. Returning data unchanged.
#> Normalized 4/4 organisms (100.0%)
#> Result: 3 unique organisms
#> is_MRSA: 0
#> is_MSSA: 0
#> is_MRCONS: 0
#> is_MSCONS: 0
#> Normalizing 4 unique antibiotic names against WHO reference (257 antibiotics)...
#> Normalized: 4 unique names -> 4
#> With antibiotic_class: 4 (100.0%)
#> With AWaRe category: 4 (100.0%)
#> Specimen normalization: 4/4 matched (100.0%)
#> Sterile classification: 3 sterile, 1 non-sterile, 0 others/ambiguous
#> 
#> Age enrichment summary:
#>   age_method age_confidence n
#> 1   provided           high 4
#> [prep_derive_los_from_dates] 0 filled; 4 remain missing.
#> 
#> LOS summary:
#>   mean_los median_los min_los max_los
#> 1       NA         NA      NA      NA
#> Inferring infection type using 2-day HAI cutoff...
#> Enriched infection_type: 4 rows filled
#> 
#> Infection type distribution:
#>   infection_type infection_type_method n
#> 1            HAI  inferred_2day_cutoff 3
#> 2            CAI  inferred_2day_cutoff 1
#> Inferring hospital department from contextual data...
#> Enriched department: 3 rows filled (LOW confidence - heuristic)
#> 
#> Department distribution:
#>   hospital_department department_confidence n
#> 1                 ICU                   low 3
#> 2                <NA>                  <NA> 1
#> Grooming 3 optional columns: infection_type, hospital_department, aware_category
#> 
#> Optional column completeness:
#>                column completeness_pct
#> 1      infection_type              100
#> 2      aware_category              100
#> 3 hospital_department               75
#> Assigned age bins: 4 binned, 0 unbinned
#> Extracted genus: 4 records, 3 unique genera
#> Extracted species: 4 records, 3 unique
#> Organism group distribution:
#> 
#>    Enterobacterales Gram-positive cocci 
#>                   3                   1 
#> Contaminants flagged: 0 total (high: 0, medium: 0)
#> Date columns not available for mortality classification
#> 
#> Mortality classification summary:
#>       Result
#> Method No
#>   <NA>  4
#> Mapped to RR pathogen: 4/4 (100.0%)
#> Mapped to RR drug: 4/4 (100.0%)
#> [prep_filter_minimally_usable] 4 rows retained; 0 removed (0.0%).
#> Warning: [prep_validate_analysis_ready] Validation issues found:
#>   - Only 4 rows present (minimum: 10).
#>  - Missing core fields: patient_id, culture_date, ast_value_harmonized.
#> [x] Quality issues detected:
#>   - 6 columns exceed 95% missing threshold: los_days, los_method, los_confidence, diagnosis_1, mortality_method
head(result$data)
#>      PID date_of_admission date_of_culture    final_outcome
#> 1 pt_001        2021-01-01      2021-01-03            alive
#> 2 pt_001        2021-01-01      2021-01-03            alive
#> 3 pt_002        2021-01-05      2021-01-08          expired
#> 4 pt_003        2021-01-10      2021-01-11 Discharged Alive
#>           organism_name specimen_type antibiotic_name antibiotic_value gender
#> 1               E. coli         Blood        Amikacin                R      M
#> 2               E. coli         Blood   Ciprofloxacin                S      M
#> 3 Klebsiella pneumoniae         Urine       Meropenem                R      F
#> 4 Staphylococcus aureus         Blood      Vancomycin                S      F
#>   Age outcome_std is_MRSA is_MSSA is_MRCONS is_MSCONS   organism_normalized
#> 1  34    Survived       0       0         0         0      escherichia coli
#> 2  34    Survived       0       0         0         0      escherichia coli
#> 3  61        Died       0       0         0         0 klebsiella pneumoniae
#> 4  47    Survived       0       0         0         0 staphylococcus aureus
#>        organism_group antibiotic_normalized antibiotic_class aware_category
#> 1    Enterobacterales              amikacin  Aminoglycosides         Access
#> 2    Enterobacterales         ciprofloxacin Fluoroquinolones          Watch
#> 3    Enterobacterales             meropenem      Carbapenems          Watch
#> 4 Gram-positive cocci         vancomycin_iv       Vancomycin          Watch
#>   specimen_normalized sample_category sterile_classification age_method
#> 1               Blood           Blood                Sterile   provided
#> 2               Blood           Blood                Sterile   provided
#> 3               Urine           Urine            Non-Sterile   provided
#> 4               Blood           Blood                Sterile   provided
#>   age_confidence los_days los_method los_confidence infection_type
#> 1           high       NA       <NA>           <NA>            HAI
#> 2           high       NA       <NA>           <NA>            HAI
#> 3           high       NA       <NA>           <NA>            HAI
#> 4           high       NA       <NA>           <NA>            CAI
#>   infection_type_method infection_type_confidence hospital_department
#> 1  inferred_2day_cutoff                    medium                 Icu
#> 2  inferred_2day_cutoff                    medium                 Icu
#> 3  inferred_2day_cutoff                    medium                <NA>
#> 4  inferred_2day_cutoff                    medium                 Icu
#>   diagnosis_1  department_method department_confidence Age_years Age_bin
#> 1        <NA> inferred_heuristic                   low        34   30-35
#> 2        <NA> inferred_heuristic                   low        34   30-35
#> 3        <NA>               <NA>                  <NA>        61   60-65
#> 4        <NA> inferred_heuristic                   low        47   45-50
#>   Age_resolved      org_genus org_species           org_group is_contaminant
#> 1           34    escherichia        coli    Enterobacterales          FALSE
#> 2           34    escherichia        coli    Enterobacterales          FALSE
#> 3           61     klebsiella  pneumoniae    Enterobacterales          FALSE
#> 4           47 staphylococcus      aureus Gram-positive cocci          FALSE
#>   contaminant_confidence contaminant_method mortality_infection
#> 1                    low          heuristic                  No
#> 2                    low          heuristic                  No
#> 3                    low          heuristic                  No
#> 4                    low          heuristic                  No
#>   mortality_method mortality_confidence           rr_pathogen          rr_drug
#> 1             <NA>                 <NA>      escherichia coli  Aminoglycosides
#> 2             <NA>                 <NA>      escherichia coli Fluoroquinolones
#> 3             <NA>                 <NA> klebsiella pneumoniae      Carbapenems
#> 4             <NA>                 <NA> staphylococcus aureus       Vancomycin
```
