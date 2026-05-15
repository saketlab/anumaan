# Normalize Organism Names

Standardizes organism names using organisms.csv reference file.
Automatically handles abbreviations (E. coli), case variations, and
typos.

## Usage

``` r
prep_standardize_organisms(
  data,
  organism_col = "organism_name",
  add_organism_group = TRUE,
  add_resistance_flags = TRUE
)
```

## Arguments

- data:

  Data frame with organism column.

- organism_col:

  Character. Organism column name. Default "organism_name".

- add_organism_group:

  Logical. Add organism_group column. Default TRUE.

- add_resistance_flags:

  Logical. Add MRSA/MRCONS flags. Default TRUE.

## Value

Data frame with organism_normalized, organism_group, and optionally
resistance flag columns.
