# Harmonize AST Values

Wrapper that applies the full AST cleaning pipeline in sequence:

1.  [`prep_clean_ast_values()`](https://saketlab.github.io/anumaan/reference/prep_clean_ast_values.md) -
    extract clean S/I/R from raw strings (handles exact matches, regex
    patterns, and keyword scanning in one pass)

2.  [`prep_recode_intermediate_ast()`](https://saketlab.github.io/anumaan/reference/prep_recode_intermediate_ast.md) -
    resolve I values

Retains the original `ast_col` column as `ast_value_raw` and adds
`ast_value_harmonized` and `intermediate_recoded`.

## Usage

``` r
prep_harmonize_ast(
  data,
  antibiotic_col = "antibiotic_name_std",
  ast_col = "ast_value_raw",
  colistin_to_s = TRUE,
  others_to_r = TRUE
)
```

## Arguments

- data:

  Data frame with AST data.

- antibiotic_col:

  Character. Standardized antibiotic name column. Default
  "antibiotic_name_std".

- ast_col:

  Character. Raw AST value column. Default "ast_value_raw".

- colistin_to_s:

  Logical. Recode Colistin I -\> S. Default TRUE.

- others_to_r:

  Logical. Recode non-Colistin I -\> R. Default TRUE.

## Value

Data frame with `ast_value_harmonized` and `intermediate_recoded`
columns added.
