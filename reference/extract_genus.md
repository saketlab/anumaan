# Extract Genus from Organism Name

Extracts bacterial genus (first word) from normalized organism name.

## Usage

``` r
extract_genus(data, organism_col = "organism_normalized")
```

## Arguments

- data:

  Data frame

- organism_col:

  Character. Organism column name. Default "organism_normalized".

## Value

Data frame with org_genus column added
