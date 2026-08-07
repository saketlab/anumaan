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
if (FALSE) { # \dontrun{
result <- preprocess_for_profiles(
  data = raw_ast,
  col_map = list(
    isolate_col = "isolate_id",
    pathogen_col = "organism",
    ast_col = "result",
    patient_col = "pid",
    date_col = "culture_date",
    geography_col = "state",
    specimen_col = "specimen",
    age_col = "age",
    antibiotic_col = "drug_name"
  ),
  panel_map = list(
    "Klebsiella pneumoniae" = c("3GC", "Carbapenems", "Fluoroquinolones"),
    "Escherichia coli"      = c("3GC", "Fluoroquinolones", "Aminoglycosides")
  ),
  stratify_by = c("geography", "year"),
  outcome_col = "final_outcome"
)

result$data_wide # ready for compute_marginals_from_data()
result$preprocessing_log # step-by-step row counts
result$panel_exclusions # classes dropped per organism
} # }
```
