# Classify Organism Group

Assigns organisms to taxonomic groups (Enterobacterales, Gram-positive,
etc.)

## Usage

``` r
classify_org_group(data, organism_col = "organism_normalized")
```

## Arguments

- data:

  Data frame

- organism_col:

  Character. Normalized organism column. Default "organism_normalized".

## Value

Data frame with org_group column added
