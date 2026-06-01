# Handle Duplicate AST Results

Deduplicates AST records in two sequential steps:

1.  **Exact-duplicate removal**: drops rows that are fully identical
    across all five key columns (patient + organism + antibiotic +
    date + value). These are true redundant records with no ambiguity.

2.  **Conflict resolution**: handles the remaining cases where the same
    patient + organism + antibiotic + date group has *different* values.

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

  After exact-dedup, flags conflicting rows with
  `is_ast_duplicate = TRUE` and prints a QC summary. Returns the data
  with the flag column so you can inspect before deciding how to
  resolve.

- `"remove"`:

  Runs both steps: removes exact duplicates, then applies `strategy` to
  keep one row per key combination for any remaining conflicts.
