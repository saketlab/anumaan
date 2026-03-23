# Classify XDR (Extensively Drug Resistant)

Classifies isolates as XDR using Magiorakos 2012 criteria.

## Usage

``` r
classify_xdr(data, definition = "Magiorakos", organism_group_col = "org_group")
```

## Arguments

- data:

  Data frame

- definition:

  Character. "Magiorakos" or "WHO". Default "Magiorakos".

- organism_group_col:

  Character. Organism group column.

## Value

Data frame with xdr, xdr_confidence, xdr_method columns

## References

Magiorakos AP et al. Clin Microbiol Infect. 2012;18(3):268-281.
