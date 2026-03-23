# Add RR Pathogen and Drug Mappings

Adds rr_pathogen and rr_drug columns to data by matching against the GBD
RR reference list.

## Usage

``` r
add_rr_mappings(
  data,
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_class"
)
```

## Arguments

- data:

  Data frame with organism and antibiotic columns

- organism_col:

  Character. Organism column. Default "organism_normalized".

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_class".

## Value

Data frame with rr_pathogen and rr_drug columns added

## Examples

``` r
if (FALSE) { # \dontrun{
data <- add_rr_mappings(data)
} # }
```
