# Detect Schema Drift Across Centres

Compares column sets across a list of centre data frames. Flags columns
present in some centres but not others.

## Usage

``` r
prep_detect_schema_drift(data_list, reference_centre = NULL)
```

## Arguments

- data_list:

  Named list of data frames (one per centre).

- reference_centre:

  Character or NULL. Name of the centre to treat as the reference
  schema. If NULL, the union of all columns is used.

## Value

Tibble: column \| present_in (comma-separated centres) \| missing_from
\| n_centres_present
