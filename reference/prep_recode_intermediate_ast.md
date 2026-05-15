# Recode Intermediate (I) Susceptibility Values

Converts "I" (Intermediate) values to either "S" or "R" based on
antibiotic. For Colistin: I -\> S (clinical guideline). For all other
antibiotics: I -\> R (conservative surveillance approach).

## Usage

``` r
prep_recode_intermediate_ast(
  data,
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  colistin_to_s = TRUE,
  others_to_r = TRUE
)
```

## Arguments

- data:

  Data frame with antibiotic susceptibility data.

- antibiotic_col:

  Character. Column with antibiotic names. Default "antibiotic_name".

- value_col:

  Character. Column with S/I/R values. Default "antibiotic_value".

- colistin_to_s:

  Logical. Convert Colistin I to S. Default TRUE.

- others_to_r:

  Logical. Convert other antibiotics' I to R. Default TRUE.

## Value

Data frame with recoded intermediate values.
