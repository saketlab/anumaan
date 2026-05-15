# Handle Duplicate AST Results

Detects and optionally resolves conflicting AST records for the same
patient + organism + antibiotic + date combination.

## Usage

``` r
prep_deduplicate_ast(
  data,
  mode = c("detect", "remove"),
  strategy = c("resistant_wins", "susceptible_wins", "first"),
  patient_col = "patient_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized",
  date_col = "culture_date",
  ast_col = "ast_value_harmonized"
)
```

## Arguments

- data:

  Data frame with AST data in long format.

- mode:

  Character. `"detect"` (flag only) or `"remove"` (flag then resolve).
  Default `"detect"`.

- strategy:

  Character. Resolution strategy used only when `mode = "remove"`. One
  of `"resistant_wins"` (default), `"susceptible_wins"`, or `"first"`.

- patient_col:

  Character. Patient ID column. Default `"patient_id"`.

- organism_col:

  Character. Organism column. Default `"organism_normalized"`.

- antibiotic_col:

  Character. Antibiotic column. Default `"antibiotic_normalized"`.

- date_col:

  Character. Culture date column. Default `"culture_date"`.

- ast_col:

  Character. Harmonized AST value column. Default
  `"ast_value_harmonized"`.

## Value

- `mode = "detect"`: original data with `is_ast_duplicate` logical
  column added.

- `mode = "remove"`: data with conflicts resolved (one row per key) and
  no flag column.

## Details

- `"detect"`:

  Flags conflicting rows with `is_ast_duplicate = TRUE` and prints a QC
  summary of all conflict groups. Returns the full data frame with the
  flag column so you can inspect or filter before deciding how to
  resolve.

- `"remove"`:

  Runs the detect step first (flag + QC summary), then applies
  `strategy` to keep one row per key combination and drops the flag
  column from the returned data.
