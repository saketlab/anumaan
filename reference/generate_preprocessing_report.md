# Generate Preprocessing Report

Creates a comprehensive report documenting all transformations applied
during preprocessing, including column mappings, data quality metrics,
and transformation summaries.

## Usage

``` r
generate_preprocessing_report(
  raw_data,
  processed_data,
  preprocessing_log,
  config,
  output_file = NULL,
  include_plots = TRUE
)
```

## Arguments

- raw_data:

  Original input data frame (before preprocessing)

- processed_data:

  Final output data frame (after preprocessing)

- preprocessing_log:

  List containing logs from each preprocessing step

- config:

  AMR configuration object used

- output_file:

  Optional. Path to save report. If NULL, returns list. Supports .html,
  .pdf, .txt, .rds formats.

- include_plots:

  Logical. Include data quality plots. Default TRUE.

## Value

If output_file is NULL, returns a list with report components.
Otherwise, saves report to file and returns file path.

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
report <- generate_preprocessing_report(raw, result$data, result$log,
  result$config,
  include_plots = FALSE
)
report
#> ===================================================
#>    AMR Data Preprocessing Report
#> ===================================================
#> 
#> Report Generated: 2026-10-20 05:11:49 
#> Package Version: 0.1.0.9037 
#> R Version: R version 4.6.1 (2026-06-24) 
#> 
#> --- Raw Data Summary ---
#>   Rows: 4
#>   Columns: 10
#>   Memory: 3.9 Kb
#>   Missing Data: 0.0%
#> 
#> --- Column Mapping ---
#>   derived: 47 columns
#>   dropped: 8 columns
#>   exact_match: 8 columns
#>   fuzzy_match: 1 columns
#>   no_mapping: 1 columns
#> 
#> --- Processing Statistics ---
#>   Runtime: NA seconds
#>   Rows: 4 -> 4 (100.0% retained)
#>   Columns: 10 -> 49 (+39 derived)
#>   Variables Derived: 12
#> 
#> --- Data Quality ---
#>   Completeness: 100.0% -> 86.2%
#>   Duplicate Rows: 0 -> 0
#> 
#> ===================================================
#> Use summary(report) for detailed sections
#> Use export_report(report, 'file.html') to save
#> ===================================================
```
