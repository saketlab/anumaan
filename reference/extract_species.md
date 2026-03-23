# Extract Species from Organism Name

Extracts bacterial species (second word) from normalized organism name.

## Usage

``` r
extract_species(data, organism_col = "organism_normalized")
```

## Arguments

- data:

  Data frame

- organism_col:

  Character. Organism column name. Default "organism_normalized".

## Value

Data frame with org_species column added
