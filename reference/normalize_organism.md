# Normalize Organism Names

Normalizes organism names using organisms.csv reference file.
Automatically handles abbreviations (E. coli), case variations, and
typos.

## Usage

``` r
normalize_organism(
  data,
  organism_col = "organism_name",
  add_organism_group = TRUE,
  add_resistance_flags = TRUE
)
```

## Arguments

- data:

  Data frame with organism column

- organism_col:

  Character. Organism column name. Default "organism_name".

- add_organism_group:

  Logical. Add organism_group column from CSV. Default TRUE.

- add_resistance_flags:

  Logical. Add resistance flag columns (0/1). Default TRUE. Creates:
  is_MRSA, is_MSSA, is_MRCONS, is_MSCONS

## Value

Data frame with organism_normalized, organism_group, and optionally
resistance flag columns

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(organism = c("E. coli", "S. aureus", "MRSA"))
result <- normalize_organism(data, organism_col = "organism")
} # }
```
