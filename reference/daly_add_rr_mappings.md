# Add RR Pathogen and Drug Class Mappings

Maps organism names to GBD RR pathogen categories and antibiotic class
names to RR drug categories in a single pass. Either mapping is skipped
silently when the corresponding column is absent from `data`.

## Usage

``` r
daly_add_rr_mappings(
  data,
  organism_col = "organism_normalized",
  class_col = "antibiotic_class"
)
```

## Arguments

- data:

  Data frame with organism and/or antibiotic class columns.

- organism_col:

  Character. Column with organism names to map to `rr_pathogen`. Set to
  `NULL` to skip. Default `"organism_normalized"`.

- class_col:

  Character. Column with antibiotic class names to map to `rr_drug`. Set
  to `NULL` to skip. Default `"antibiotic_class"`.

## Value

Data frame with `rr_pathogen` and/or `rr_drug` columns appended.

## Details

Pathogen mapping uses the organism taxonomy
([`get_organism_taxonomy()`](https://saketlab.github.io/anumaan/reference/get_organism_taxonomy.md));
drug class mapping uses `inst/extdata/WHO_aware_class.csv`.
