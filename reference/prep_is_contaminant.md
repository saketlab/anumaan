# Check if Organism is a Contaminant

Checks whether an organism name matches any known contaminant for a
specific syndrome or specimen type using flexible pattern matching.

## Usage

``` r
prep_is_contaminant(organism_name, syndrome = NULL, specimen_type = NULL)
```

## Arguments

- organism_name:

  Character vector of organism name(s) to check.

- syndrome:

  Character. Optional syndrome filter.

- specimen_type:

  Character. Optional specimen type filter.

## Value

Logical vector (same length as `organism_name`).
