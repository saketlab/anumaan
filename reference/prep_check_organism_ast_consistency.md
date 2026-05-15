# Check Organism-AST Consistency

Flags impossible/implausible AST results where a drug is inappropriate
for the organism type (e.g., a Gram-positive organism tested against a
Gram-negative-only drug).

## Usage

``` r
prep_check_organism_ast_consistency(
  data,
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized"
)
```

## Arguments

- data:

  Data frame with organism and antibiotic columns.

- organism_col:

  Character. Normalized organism column. Default "organism_normalized".

- antibiotic_col:

  Character. Normalized antibiotic column. Default
  "antibiotic_normalized".

## Value

Data frame with `is_ast_inconsistent` logical column added.

## Details

Uses `inst/extdata/Organism_AST_drugs.csv` which maps organism groups to
expected antibiotic panels.
