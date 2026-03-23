# Map Organism to RR Pathogen Category

Maps normalized organism names to RR (Relative Risk) pathogen categories
used in burden estimation.

## Usage

``` r
map_to_rr_pathogen(data, organism_col = "organism_normalized")
```

## Arguments

- data:

  Data frame

- organism_col:

  Character. Normalized organism column.

## Value

Data frame with rr_pathogen column added
