# Filter Analysis-Ready Records

Keeps records that meet the full burden-analysis readiness criteria:

- Patient ID present

- Culture date present

- Organism name present and standardized

- Antibiotic name present and standardized

- Harmonized AST value in S/I/R

- Not flagged as a contaminant (if column present)

## Usage

``` r
prep_filter_analysis_ready(
  data,
  patient_col = "patient_id",
  culture_date_col = "culture_date",
  organism_col = "organism_name",
  antibiotic_col = "antibiotic_name",
  ast_col = "ast_value_harmonized",
  exclude_contaminants = TRUE,
  contaminant_col = "contaminant_flag"
)
```

## Arguments

- data:

  Data frame after preprocessing.

- patient_col:

  Character. Patient ID column. Default "patient_id".

- culture_date_col:

  Character. Culture date column. Default "culture_date".

- organism_col:

  Character. Organism column. Default "organism_name".

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_name".

- ast_col:

  Character. Harmonized AST column. Default "ast_value_harmonized".

- exclude_contaminants:

  Logical. Exclude rows flagged as contaminants. Default TRUE.

- contaminant_col:

  Character. Contaminant flag column. Default "contaminant_flag".

## Value

Filtered data frame.
