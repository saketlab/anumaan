# Collapse to Class Level

Aggregates resistance at antibiotic class level instead of individual
drugs. Uses "any R in class -\> class R" logic.

## Usage

``` r
prep_collapse_class_level(
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
ast <- data.frame(
  event_id              = rep(c("E1", "E2"), each = 3),
  organism_normalized   = rep(c("Escherichia coli", "Klebsiella pneumoniae"), each = 3),
  antibiotic_normalized = rep(c("Meropenem", "Imipenem", "Ciprofloxacin"), 2),
  antibiotic_class      = rep(c("Carbapenems", "Carbapenems", "Fluoroquinolones"), 2),
  antibiotic_value      = c("S", "R", "R", "S", "S", "R")
)
prep_collapse_class_level(ast)
#> Collapsing to antibiotic class level...
#> Collapsed: 6 rows -> 4 class-level rows
#> 
#> Class-level resistance distribution:
#> # A tibble: 2 × 3
#>   antibiotic_class     R     S
#>   <chr>            <int> <int>
#> 1 Carbapenems          1     1
#> 2 Fluoroquinolones     2     0
#> # A tibble: 4 × 9
#>   event_id organism_normalized   antibiotic_class class_resistance
#>   <chr>    <chr>                 <chr>            <chr>           
#> 1 E1       Escherichia coli      Carbapenems      R               
#> 2 E1       Escherichia coli      Fluoroquinolones R               
#> 3 E2       Klebsiella pneumoniae Carbapenems      S               
#> 4 E2       Klebsiella pneumoniae Fluoroquinolones R               
#> # ℹ 5 more variables: n_drugs_in_class <int>, n_resistant <int>,
#> #   pct_resistant_in_class <dbl>, drugs_tested <chr>, collapse_method <chr>
```
