# Flag Organisms Unmatched in Reference

Flags organisms that could not be matched to the reference
organisms.csv. Adds a logical column `is_organism_unmatched`.

## Usage

``` r
prep_flag_organism_unmatched(data, organism_col = "organism_normalized")
```

## Arguments

- data:

  Data frame with organism_normalized column.

- organism_col:

  Character. Standard organism column. Default "organism_normalized".

## Value

Data frame with `is_organism_unmatched` column.
