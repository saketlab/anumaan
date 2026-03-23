# Clean Antibiotic Susceptibility Values

Extracts clean S/I/R values from messy antibiotic result columns.
Handles common patterns like "S (HIGH LEVEL)", "R Escherichia coli",
etc.

## Usage

``` r
clean_antibiotic_values(data, value_col = "antibiotic_value", strict = FALSE)
```

## Arguments

- data:

  Data frame with antibiotic susceptibility data

- value_col:

  Character. Column name containing susceptibility values. Default
  "antibiotic_value".

- strict:

  Logical. If TRUE, only accept S/I/R values. If FALSE, attempt to parse
  from messy strings. Default FALSE.

## Value

Data frame with cleaned antibiotic_value column

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  antibiotic_value = c("S", "R   E. coli", "S (HIGH LEVEL)", "I")
)
clean_antibiotic_values(data)
} # }
```
