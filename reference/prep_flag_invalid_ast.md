# Flag Invalid AST Values

Adds a logical column `is_ast_invalid` that is TRUE when the harmonized
AST value is not in `c("S", "I", "R")` and not NA.

## Usage

``` r
prep_flag_invalid_ast(data, col = "ast_value_harmonized")
```

## Arguments

- data:

  Data frame.

- col:

  Character. Harmonized AST column. Default "ast_value_harmonized".

## Value

Data frame with `is_ast_invalid` column added.
