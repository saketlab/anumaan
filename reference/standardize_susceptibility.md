# Standardize Susceptibility Values

Maps various susceptibility/resistance representations to standard "S"
(Susceptible), "R" (Resistant), "I" (Intermediate) values.

## Usage

``` r
standardize_susceptibility(data, col = "antibiotic_value")
```

## Arguments

- data:

  Data frame containing susceptibility column

- col:

  Character. Name of susceptibility column. Default "antibiotic_value".

## Value

Data frame with standardized susceptibility values and added column
`antibiotic_value_std`
