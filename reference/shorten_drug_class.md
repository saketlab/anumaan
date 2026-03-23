# Shorten Antibiotic Class Names

Maps long antibiotic class names to common abbreviations used in
GBD-style figures.

## Usage

``` r
shorten_drug_class(x)
```

## Arguments

- x:

  Character vector of antibiotic class names.

## Value

Character vector of shortened names.

## Examples

``` r
shorten_drug_class(c("Carbapenems", "Third-generation-cephalosporins"))
#> [1] "Carbapenems" "3GC"        
```
