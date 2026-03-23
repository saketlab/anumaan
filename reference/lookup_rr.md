# Lookup Relative Risk Values

Looks up RR values from table. Implements 3GC/4GC proxy logic.

## Usage

``` r
lookup_rr(data, rr_table = NULL)
```

## Arguments

- data:

  Data frame with rr_pathogen and rr_drug columns

- rr_table:

  Data frame with RR values. If NULL, uses built-in.

## Value

Data frame with rr_value column added
