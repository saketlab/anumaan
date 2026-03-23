# Classify MDR (Multidrug Resistant)

Classifies isolates as MDR using Magiorakos 2012 criteria.

## Usage

``` r
classify_mdr(data, definition = "Magiorakos", organism_group_col = "org_group")
```

## Arguments

- data:

  Data frame

- definition:

  Character. "Magiorakos" or "WHO". Default "Magiorakos".

- organism_group_col:

  Character. Organism group column for pathogen-specific thresholds.

## Value

Data frame with mdr, mdr_confidence, mdr_method, n_resistant_categories
columns

## References

Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
