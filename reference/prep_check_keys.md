# Check join key quality

Reports missing, blank, duplicate, and placeholder keys in a column
before a merge. Run this on both sides of every planned join.

## Usage

``` r
prep_check_keys(
  data,
  key_col,
  table_label = "table",
  warn_missing_pct = 5,
  admission_col = "date_of_admission",
  culture_col = "date_of_culture"
)
```

## Arguments

- data:

  Data frame.

- key_col:

  Character. Name of the key column.

- table_label:

  Character. Label used in messages.

- warn_missing_pct:

  Numeric. Warn when proportion of missing keys exceeds this threshold
  (0-100). Default 5.

- admission_col:

  Character. Admission date column used for the optional paired
  missingness check. Default "date_of_admission".

- culture_col:

  Character. Culture date column used for the optional paired
  missingness check. Default "date_of_culture".

## Value

Invisibly returns a one-row summary tibble.
