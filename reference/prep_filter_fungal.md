# Filter Fungal Organisms

Removes rows where the organism belongs to a fungal group. Logs the
count of rows removed.

## Usage

``` r
prep_filter_fungal(
  data,
  group_col = "organism_group",
  organism_col = "organism_normalized"
)
```

## Arguments

- data:

  Data frame.

- group_col:

  Character. Organism group column. Default "organism_group".

- organism_col:

  Character. Normalized organism column used as fallback when
  `group_col` is absent. Default "organism_normalized".

## Value

Data frame with fungal rows removed.
