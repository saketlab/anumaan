# Recode Intermediate (I) Susceptibility Values

Converts "I" (Intermediate) values to either "S" or "R" based on
antibiotic type. For Colistin: I -\> S (following clinical guidelines)
For all other antibiotics: I -\> R (conservative approach for
surveillance)

## Usage

``` r
recode_intermediate(
  data,
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  colistin_to_s = TRUE,
  others_to_r = TRUE
)
```

## Arguments

- data:

  Data frame with antibiotic susceptibility data

- antibiotic_col:

  Character. Column name with antibiotic names. Default
  "antibiotic_name".

- value_col:

  Character. Column name with S/I/R values. Default "antibiotic_value".

- colistin_to_s:

  Logical. Convert Colistin I to S. Default TRUE.

- others_to_r:

  Logical. Convert other antibiotics' I to R. Default TRUE.

## Value

Data frame with recoded intermediate values

## Examples

``` r
if (FALSE) { # \dontrun{
# Recode all I values according to standard rules
data <- recode_intermediate(data)

# Only recode Colistin
data <- recode_intermediate(data, others_to_r = FALSE)
} # }
```
