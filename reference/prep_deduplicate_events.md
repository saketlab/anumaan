# Deduplicate Events

Removes duplicate rows within groups. Two modes:

## Usage

``` r
prep_deduplicate_events(
  data,
  event_col = "event_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized",
  key_cols = NULL,
  keep = "first"
)
```

## Arguments

- data:

  Data frame.

- event_col:

  Character. Event ID column (event-aware mode). Default "event_id".

- organism_col:

  Character. Organism column (event-aware mode). Default
  "organism_normalized".

- antibiotic_col:

  Character. Antibiotic column (event-aware mode). Default
  "antibiotic_normalized".

- key_cols:

  Character vector. When supplied, switches to generic mode and uses
  these columns for duplicate detection. NULL uses event-aware mode.
  Default NULL.

- keep:

  Character. "first", "last", or "none" (drop all duplicates). Default
  "first".

## Value

Deduplicated data frame.

## Details

- Event-aware (default):

  Groups by `event_col` + `organism_col` + `antibiotic_col` and keeps
  first/last antibiotic test per group. Requires `event_col` to be
  present.

- Generic (key_cols supplied):

  Groups by `key_cols` and drops duplicates. When `keep = "none"`, both
  copies of every duplicate are removed. Replaces the former
  `remove_duplicate_rows()` helper.
