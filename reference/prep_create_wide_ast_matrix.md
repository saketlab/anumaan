# Create Wide Format AST Matrix

Converts long format (one row per organism-antibiotic) to wide format
(one row per event with antibiotic columns). Useful for analysis and
machine learning applications.

## Usage

``` r
prep_create_wide_ast_matrix(
  data,
  event_col = "event_id",
  antibiotic_col = "antibiotic_normalized",
  susceptibility_col = "antibiotic_value",
  prefix = "abx_",
  keep_cols = c("patient_id", "organism_normalized", "date_of_culture")
)
```

## Arguments

- data:

  Data frame in long format

- event_col:

  Character. Event ID column. Default "event_id".

- antibiotic_col:

  Character. Antibiotic/class column to pivot. Default
  "antibiotic_normalized".

- susceptibility_col:

  Character. Susceptibility column. Default "antibiotic_value".

- prefix:

  Character. Prefix for pivoted columns. Default "abx\_".

- keep_cols:

  Character vector. Additional columns to keep from original data.
  Default c("patient_id", "organism_normalized", "date_of_culture").

## Value

Wide format data frame (one row per event)

## Examples

``` r
ast <- data.frame(
  event_id              = rep(c("E1", "E2"), each = 3),
  organism_normalized   = rep(c("Escherichia coli", "Klebsiella pneumoniae"), each = 3),
  antibiotic_normalized = rep(c("Meropenem", "Imipenem", "Ciprofloxacin"), 2),
  antibiotic_class      = rep(c("Carbapenems", "Carbapenems", "Fluoroquinolones"), 2),
  antibiotic_value      = c("S", "R", "R", "S", "S", "R")
)
prep_create_wide_ast_matrix(ast, keep_cols = "organism_normalized")
#> Creating wide format: pivoting 'antibiotic_normalized' column...
#> Created wide format: 2 events x 3 antibiotics
#> # A tibble: 2 × 5
#>   event_id organism_normalized   abx_Ciprofloxacin abx_Imipenem abx_Meropenem
#>   <chr>    <chr>                 <chr>             <chr>        <chr>        
#> 1 E1       Escherichia coli      R                 R            S            
#> 2 E2       Klebsiella pneumoniae R                 S            S            
```
