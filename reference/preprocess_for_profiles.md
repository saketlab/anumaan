# Preprocess AST Data for Resistance Profile Estimation (Pathway 1)

Orchestrates the full preprocessing pipeline required before marginal
resistance computation and convex optimisation. Calls existing `prep_*`
package functions in sequence and returns a wide tibble (one row per
isolate) with antibiotic-class columns ready for
[`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md).

## Usage

``` r
preprocess_for_profiles(
  data,
  col_map = list(isolate_col = "isolate_id", pathogen_col = "pathogen", ast_col =
    "ast_value", patient_col = "patient_id", date_col = "date_of_culture", geography_col
    = "state", specimen_col = "specimen_type", age_col = "age", dob_col = "dob",
    antibiotic_col = "antibiotic_name", class_col = "antibiotic_class", location_col =
    NULL, outcome_col = NULL),
  panel_map,
  stratify_by = NULL,
  outcome_col = NULL,
  who_table = NULL
)
```

## Arguments

- data:

  Data frame. Raw AST data in long or wide format.

- col_map:

  Named list mapping logical roles to actual column names. Same
  structure as in
  [`validate_profile_inputs()`](https://saketlab.github.io/anumaan/reference/validate_profile_inputs.md).
  Default names follow standard package conventions.

- panel_map:

  Named list. **Mandatory.** Maps each pathogen name to a character
  vector of antibiotic class names to include in profile estimation.
  Classes not in the panel are excluded and logged with reason code
  `"not_in_user_panel"`. Example:
  `list("Klebsiella pneumoniae" = c("3GC", "Carbapenems", "Fluoroquinolones"))`.

- stratify_by:

  Character vector or `NULL`. Stratification dimensions to carry through
  to output. Valid values: `"geography"`, `"year"`, `"outcome"`. Each
  requested dimension adds the corresponding column to the wide output.
  Default `NULL`.

- outcome_col:

  Character or `NULL`. Column name for patient outcome (e.g.
  `"final_outcome"`). Required when `"outcome"` is in `stratify_by`.
  Default `NULL`.

- who_table:

  Data frame or `NULL`. Custom WHO AWaRe classification table passed to
  [`prep_standardize_antibiotics()`](https://saketlab.github.io/anumaan/reference/prep_standardize_antibiotics.md).
  If `NULL`, uses the built-in `inst/extdata/WHO_aware_class.csv`.
  Default `NULL`.

## Value

A named list:

- `data_wide`:

  Tibble. One row per `isolate_id`. Columns: all metadata columns + one
  column per antibiotic class (values S / R / NA).

- `preprocessing_log`:

  Tibble. One row per pipeline step with columns `step`, `n_rows_in`,
  `n_rows_out`, `message`, `status`.

- `panel_exclusions`:

  Tibble. All organism-class combinations removed by the panel filter,
  with `reason_code` and `reason_text`.

- `col_map_resolved`:

  Named list. Updated `col_map` with any column names derived during
  preprocessing (e.g. the harmonised AST column name, the resolved class
  column name).

- `detected_format`:

  Character. `"long"` or `"wide"`.

## Details

**Pipeline steps (in order):**

1.  Input validation via
    [`validate_profile_inputs()`](https://saketlab.github.io/anumaan/reference/validate_profile_inputs.md).

2.  Wide -\> long conversion if wide format detected
    ([`prep_pivot_ast_wide_to_long()`](https://saketlab.github.io/anumaan/reference/prep_pivot_ast_wide_to_long.md)).

3.  AST value harmonisation – clean messy strings, resolve I values
    ([`prep_harmonize_ast()`](https://saketlab.github.io/anumaan/reference/prep_harmonize_ast.md)).

4.  Antibiotic name standardisation if no class column
    ([`prep_standardize_antibiotics()`](https://saketlab.github.io/anumaan/reference/prep_standardize_antibiotics.md)).

5.  Class assignment if still absent
    ([`prep_classify_antibiotic_class()`](https://saketlab.github.io/anumaan/reference/prep_classify_antibiotic_class.md)).

6.  Class-level collapse: any R in class =\> class = R
    ([`prep_collapse_class_level()`](https://saketlab.github.io/anumaan/reference/prep_collapse_class_level.md)).
    Counting unit is `isolate_id`.

7.  Year derivation from `date_col` when `"year"` is in `stratify_by`.

8.  Organism-specific panel filter: retain only user-specified classes
    per pathogen.

9.  Long -\> wide pivot with antibiotic-class columns
    ([`prep_create_wide_ast_matrix()`](https://saketlab.github.io/anumaan/reference/prep_create_wide_ast_matrix.md)).

**Counting unit:** All statistics are computed per unique `isolate_id`.
A patient with multiple isolates contributes one count per isolate, not
one count per patient.

## Examples

``` r
ast_long <- data.frame(
  isolate_id = rep(sprintf("ISO%02d", 1:6), each = 2),
  patient_id = rep(sprintf("P%02d", 1:6), each = 2),
  pathogen = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 6),
  antibiotic_name = rep(c("Meropenem", "Ceftriaxone"), 6),
  ast_value = rep(c("R", "S", "S", "R"), 3),
  date_of_culture = as.Date("2024-01-01") + rep(0:5, each = 2),
  specimen_type = "Blood",
  state = "Maharashtra",
  age = 40
)
res <- preprocess_for_profiles(ast_long,
  panel_map = list(
    "Klebsiella pneumoniae" = c("Carbapenems", "Third-generation-cephalosporins"),
    "Escherichia coli" = c("Carbapenems", "Third-generation-cephalosporins")
  )
)
#> 
#> [preprocess_for_profiles] Step 0: Validating inputs...
#> 
#> [validate_profile_inputs] 11 passed | 1 warnings | 0 failed
#>   [!] Antibiotic class column 'antibiotic_class' absent. Will derive classes from antibiotic names via WHO AWaRe reference.
#> 
#> [preprocess_for_profiles] Step 1: Format handling...
#> 
#> [preprocess_for_profiles] Step 2: Harmonising AST values...
#> Cleaning antibiotic values: 2 unique values found
#> Cleaned: 2 unique values -> 2 (S/I/R)
#> 
#> Value distribution:
#> 
#> R S 
#> 6 6 
#> No intermediate (I) values found. Returning data unchanged.
#> [prep_harmonize_ast] 12/12 rows harmonized; 0 intermediate values recoded.
#> 
#> [preprocess_for_profiles] Step 3: Antibiotic standardisation and class assignment...
#>   No antibiotic class column found. Deriving from antibiotic names via WHO AWaRe reference.
#> Normalizing 2 unique antibiotic names against WHO reference (257 antibiotics)...
#> Normalized: 2 unique names -> 2
#> With antibiotic_class: 12 (100.0%)
#> 
#> [preprocess_for_profiles] Step 4b: Deduplicating isolate x drug rows...
#> 
#> [preprocess_for_profiles] Step 4: Collapsing to class level (any R -> class R)...
#> Collapsing to antibiotic class level...
#> Collapsed: 12 rows -> 12 class-level rows
#> 
#> Class-level resistance distribution:
#> # A tibble: 2 × 3
#>   antibiotic_class                    R     S
#>   <chr>                           <int> <int>
#> 1 Carbapenems                         3     3
#> 2 Third-generation-cephalosporins     3     3
#> 
#> [preprocess_for_profiles] Step 5: Deriving year from date column...
#> 
#> [preprocess_for_profiles] Step 6: Applying organism-specific class panels...
#> 
#> [preprocess_for_profiles] Step 7: Pivoting to wide format (classes as columns)...
#> Creating wide format: pivoting 'antibiotic_class' column...
#> Created wide format: 6 events x 2 antibiotics
#> 
#> [preprocess_for_profiles] Complete. Output: 6 isolate(s) x 10 column(s). Steps logged: 10.
res$preprocessing_log
#> # A tibble: 10 × 5
#>    step                n_rows_in n_rows_out message                       status
#>    <chr>                   <int>      <int> <chr>                         <chr> 
#>  1 0_validate                 12         12 Validation passed. Format: l… ok    
#>  2 1_format                   12         12 Long format confirmed. No pi… ok    
#>  3 2_harmonise_ast            12         12 AST values harmonised into '… ok    
#>  4 3a_standardize_abx         12         12 Antibiotic names standardise… ok    
#>  5 3b_class_assignment        12         12 Antibiotic class column reso… ok    
#>  6 4b_dedup                   12         12 Dedup: 0 duplicate isolate-d… ok    
#>  7 4_collapse_class           12         12 Class-level collapse: 12 dru… ok    
#>  8 5_year_derivation          12         12 Year derived from 'date_of_c… ok    
#>  9 6_panel_filter             12         12 Panel filter: 12 -> 12 rows.… ok    
#> 10 7_wide_pivot               12          6 Long -> wide: 6 isolates x 2… ok    
```
