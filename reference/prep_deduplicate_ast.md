# Handle Duplicate AST Results

Deduplicates AST records in two sequential steps:

1.  **Exact-duplicate removal**: drops rows that are fully identical
    across every key column plus the value column (e.g. patient +
    organism + antibiotic + date + value). These are true redundant
    records with no ambiguity.

2.  **Conflict resolution**: handles the remaining cases where the same
    key group has *different* values – by construction, every group with
    more than one row left after step 1 is a genuine conflict (same key,
    different value), never an unresolved exact duplicate.

## Usage

``` r
prep_deduplicate_ast(
  data,
  mode = c("detect", "remove"),
  strategy = c("resistant_wins", "susceptible_wins", "first", "error", "drop_conflicts"),
  group_col = NULL,
  patient_col = "patient_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized",
  date_col = "culture_date",
  ast_col = "ast_value_harmonized",
  return_diagnostics = FALSE
)
```

## Arguments

- data:

  Data frame with AST data in long format.

- mode:

  Character. `"detect"` (flag only) or `"remove"` (flag then resolve).
  Default `"detect"`.

- strategy:

  Character. Resolution strategy used only when `mode = "remove"`.
  `"resistant_wins"` (default) and `"susceptible_wins"` pick one row per
  key group by clinical priority; `"first"` keeps whichever row came
  first (only appropriate once a caller has verified row order is
  meaningful or doesn't matter); `"error"` stops instead of resolving –
  a conflicting duplicate is never silently resolved by taking the first
  row or any priority order unless that policy was explicitly chosen;
  use this when a conflict should always be surfaced, not quietly picked
  for you. `"drop_conflicts"` removes every row belonging to a
  conflicting key group entirely (never guesses which value is right) –
  the surviving data has zero rows for that key, rather than one
  unverified row.

- group_col:

  Character or `NULL`. Additional stratum identifier (e.g. hospital)
  included in the key when supplied, so that a
  patient/organism/antibiotic/date combination that legitimately recurs
  in two different strata (e.g. a hospital-local patient id reused at a
  different hospital) is never treated as one conflicting record.
  Default `NULL` (key is patient/organism/antibiotic/date only, as
  before).

- patient_col:

  Character. Patient (or event) ID column. Default `"patient_id"`.

- organism_col:

  Character or `NULL`. Organism column; `NULL` drops it from the key
  (e.g. when the data has already been subset to one organism). Default
  `"organism_normalized"`.

- antibiotic_col:

  Character. Antibiotic (or antibiotic-class) column. Default
  `"antibiotic_normalized"`.

- date_col:

  Character or `NULL`. Culture date column; `NULL` drops it from the key
  (e.g. when date has already been folded into an event id upstream).
  Default `"culture_date"`.

- ast_col:

  Character. Harmonized AST value column. Default
  `"ast_value_harmonized"`.

- return_diagnostics:

  Logical. When `TRUE`, return a list with `data` (the usual return
  value), `exact_duplicates_removed` (the rows dropped in step 1, with
  their key columns, for a caller that needs its own per-key-group
  tally), and `conflict_groups` (the conflict summary data frame, one
  row per conflicting key group). Default `FALSE` (return just the data
  frame, as before).

## Value

- `mode = "detect"`: original data with `is_ast_duplicate` logical
  column added.

- `mode = "remove"`: data with conflicts resolved (one row per key, or
  zero rows for a key group dropped by `strategy = "drop_conflicts"`)
  and no flag column.

When `return_diagnostics = TRUE`, either of the above is wrapped as
`$data` in the list described above.

## Details

The key is
`c(group_col, patient_col, organism_col, antibiotic_col, date_col)`,
with `group_col`/`organism_col`/`date_col` dropped from the key when
`NULL` – lets this function key on exactly the columns a given caller's
data actually has (e.g. a hospital/stratum identifier via `group_col`
when patient/event ids are not already globally unique across strata, or
no `organism_col`/`date_col` at all once that information is already
folded into an event id upstream).

- `"detect"`:

  After exact-dedup, flags conflicting rows with
  `is_ast_duplicate = TRUE` and prints a QC summary. Returns the data
  with the flag column so you can inspect before deciding how to
  resolve.

- `"remove"`:

  Runs both steps: removes exact duplicates, then applies `strategy` to
  resolve any remaining conflicts.
