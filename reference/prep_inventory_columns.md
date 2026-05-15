# Inventory Columns

Returns a tibble describing every column: class, distinct values, and
missingness. Used to detect completely empty columns and schema drift.

## Usage

``` r
prep_inventory_columns(data)
```

## Arguments

- data:

  Data frame.

## Value

Tibble: col_name \| class \| n_distinct \| pct_missing \| has_any_value
