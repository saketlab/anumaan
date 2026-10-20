# Convert Wide Format to Long Format

Converts wide format data (where each antibiotic is a column) to long
format (one row per organism-antibiotic combination). This is the first
step before normalization and analysis.

## Usage

``` r
prep_pivot_ast_wide_to_long(
  data,
  antibiotic_cols = NULL,
  pattern = NULL,
  id_cols = c("patient_id", "event_id", "organism_name", "date_of_culture"),
  antibiotic_name_col = "antibiotic_name",
  antibiotic_value_col = "antibiotic_value",
  remove_missing = TRUE,
  create_event_id = FALSE
)
```

## Arguments

- data:

  Data frame in wide format (antibiotics as columns)

- antibiotic_cols:

  Character vector. Names of antibiotic columns to pivot. If NULL, will
  auto-detect based on pattern. Default NULL.

- pattern:

  Character. Regex pattern to identify antibiotic columns if
  antibiotic_cols not provided. Default NULL (no auto-detection).

- id_cols:

  Character vector. Columns to keep as identifiers (not pivoted).
  Default c("patient_id", "event_id", "organism_name",
  "date_of_culture").

- antibiotic_name_col:

  Character. Name for the new column containing antibiotic names.
  Default "antibiotic_name".

- antibiotic_value_col:

  Character. Name for the new column containing susceptibility results.
  Default "antibiotic_value".

- remove_missing:

  Logical. Remove rows where antibiotic_value is NA, empty, or "-".
  Default TRUE.

- create_event_id:

  Logical. Create event_id column if it doesn't exist (uses row
  numbers). Default FALSE.

## Value

Data frame in long format

## Examples

``` r
wide <- data.frame(
  patient_id      = c("P01", "P02", "P03"),
  organism_name   = c("Escherichia coli", "Klebsiella pneumoniae", "Escherichia coli"),
  date_of_culture = as.Date(c("2024-01-02", "2024-01-05", "2024-01-09")),
  Amikacin        = c("S", "R", "S"),
  Meropenem       = c("S", "R", NA),
  Ciprofloxacin   = c("R", "R", "S")
)
long <- prep_pivot_ast_wide_to_long(wide,
  antibiotic_cols = c("Amikacin", "Meropenem", "Ciprofloxacin"),
  id_cols = c("patient_id", "organism_name", "date_of_culture"),
  create_event_id = TRUE
)
#> Creating event_id column from row numbers...
#> Pivoting 3 antibiotic columns to long format...
#> Pivoted: 3 rows -> 9 rows (wide -> long)
#> Removed 1 rows with missing/empty values
#> Final long format: 8 rows x 6 columns (3 antibiotics, 3 events)
long
#> # A tibble: 8 × 6
#>   patient_id organism_name         date_of_culture event_id antibiotic_name
#>   <chr>      <chr>                 <date>             <int> <chr>          
#> 1 P01        Escherichia coli      2024-01-02             1 Amikacin       
#> 2 P01        Escherichia coli      2024-01-02             1 Meropenem      
#> 3 P01        Escherichia coli      2024-01-02             1 Ciprofloxacin  
#> 4 P02        Klebsiella pneumoniae 2024-01-05             2 Amikacin       
#> 5 P02        Klebsiella pneumoniae 2024-01-05             2 Meropenem      
#> 6 P02        Klebsiella pneumoniae 2024-01-05             2 Ciprofloxacin  
#> 7 P03        Escherichia coli      2024-01-09             3 Amikacin       
#> 8 P03        Escherichia coli      2024-01-09             3 Ciprofloxacin  
#> # ℹ 1 more variable: antibiotic_value <chr>
```
