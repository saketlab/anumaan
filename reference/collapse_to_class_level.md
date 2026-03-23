# Collapse to Class Level

Aggregates resistance at antibiotic class level instead of individual
drugs. Uses "any R in class -\> class R" logic.

## Usage

``` r
collapse_to_class_level(
  data,
  event_col = "event_id",
  organism_col = "organism_normalized",
  class_col = "antibiotic_class",
  susceptibility_col = "antibiotic_value",
  extra_cols = NULL
)
```

## Arguments

- data:

  Data frame with antibiotic class information

- event_col:

  Character. Event ID column. Default "event_id".

- organism_col:

  Character. Organism column. Default "organism_normalized".

- class_col:

  Character. Antibiotic class column. Default "antibiotic_class".

- susceptibility_col:

  Character. Susceptibility column. Default "antibiotic_value".

- extra_cols:

  Character vector or NULL. Additional columns to carry through the
  aggregation. Default NULL.

## Value

Aggregated data frame (one row per event-organism-class)

## Examples

``` r
if (FALSE) { # \dontrun{
class_level <- collapse_to_class_level(data)
} # }
```
