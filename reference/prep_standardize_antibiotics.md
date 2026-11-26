# Normalize Antibiotic Names

Standardizes antibiotic names using fuzzy matching against WHO AWaRe
reference.

## Usage

``` r
prep_standardize_antibiotics(
  data,
  antibiotic_col = "antibiotic_name",
  who_table = NULL,
  add_class = TRUE,
  add_aware = TRUE
)
```

## Arguments

- data:

  Data frame with antibiotic data.

- antibiotic_col:

  Character. Column name with antibiotic names. Default
  "antibiotic_name".

- who_table:

  Data frame. WHO AWaRe classification table. If NULL, loads
  inst/extdata/WHO_aware_class.csv and appends
  inst/extdata/antibiotic_combination_supplement.csv (combinations found
  in surveillance data that the WHO table does not list).

- add_class:

  Logical. Add antibiotic_class column. Default TRUE.

- add_aware:

  Logical. Add aware_category column. Default TRUE.

## Value

Data frame with antibiotic_normalized, antibiotic_class, aware_category.

## Details

Combination drugs are matched as combinations: `-`, `+`, `/`, `&` and
"and" are treated as the same separator and word order is ignored, so
"piperacillin-tazobactam" matches "Piperacillin/Tazobactam" rather than
the single drug "Piperacillin". A single-drug reference only gets the
substring bonus when it covers every drug named in the input. A
combination with no reference row keeps its full name (components joined
by "/"), gets `NA` class and AWaRe category, and is reported in a
warning; it is never reduced to one of its components.
