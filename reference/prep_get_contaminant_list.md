# Get Contaminant List from Reference File

Loads contaminant organisms from `common_commensals.csv` and filters by
syndrome and/or specimen type.

## Usage

``` r
prep_get_contaminant_list(
  syndrome = NULL,
  specimen_type = NULL,
  return_all = FALSE
)
```

## Arguments

- syndrome:

  Character. Syndrome name to filter (e.g., "Bloodstream infections").

- specimen_type:

  Character. Specimen type to filter (e.g., "Blood culture").

- return_all:

  Logical. If TRUE, returns all contaminants across syndromes. Default
  FALSE.

## Value

A list with `names` (character vector) and `patterns` (list of pattern
information for flexible matching).
