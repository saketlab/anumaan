# Create Resistance Profile

Generates a resistance profile string for each event summarizing
resistance patterns. Useful for identifying common resistance
phenotypes.

## Usage

``` r
prep_create_resistance_profile(
  data,
  event_col = "event_id",
  antibiotic_col = "antibiotic_normalized",
  susceptibility_col = "antibiotic_value",
  format = "resistant_list",
  class_col = "antibiotic_class"
)
```

## Arguments

- data:

  Data frame with resistance data

- event_col:

  Character. Event ID column. Default "event_id".

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_normalized".

- susceptibility_col:

  Character. Susceptibility column. Default "antibiotic_value".

- format:

  Character. Output format: - "resistant_list": List resistant drugs
  only (default) - "full_pattern": Full S/R pattern string -
  "class_summary": Resistant classes only

- class_col:

  Character. Class column (required for format = "class_summary").
  Default "antibiotic_class".

## Value

Data frame with resistance_profile column added

## Examples

``` r
ast <- data.frame(
  event_id              = rep(c("E1", "E2"), each = 3),
  organism_normalized   = rep(c("Escherichia coli", "Klebsiella pneumoniae"), each = 3),
  antibiotic_normalized = rep(c("Meropenem", "Imipenem", "Ciprofloxacin"), 2),
  antibiotic_class      = rep(c("Carbapenems", "Carbapenems", "Fluoroquinolones"), 2),
  antibiotic_value      = c("S", "R", "R", "S", "S", "R")
)
prep_create_resistance_profile(ast)
#> Creating resistance profiles (format: resistant_list)...
#> 
#> Top 10 resistance profiles:
#>        resistance_profile n
#> 1           Ciprofloxacin 1
#> 2 Ciprofloxacin; Imipenem 1
#>   event_id   organism_normalized antibiotic_normalized antibiotic_class
#> 1       E1      Escherichia coli             Meropenem      Carbapenems
#> 2       E1      Escherichia coli              Imipenem      Carbapenems
#> 3       E1      Escherichia coli         Ciprofloxacin Fluoroquinolones
#> 4       E2 Klebsiella pneumoniae             Meropenem      Carbapenems
#> 5       E2 Klebsiella pneumoniae              Imipenem      Carbapenems
#> 6       E2 Klebsiella pneumoniae         Ciprofloxacin Fluoroquinolones
#>   antibiotic_value      resistance_profile n_resistant profile_format
#> 1                S Ciprofloxacin; Imipenem           2 resistant_list
#> 2                R Ciprofloxacin; Imipenem           2 resistant_list
#> 3                R Ciprofloxacin; Imipenem           2 resistant_list
#> 4                S           Ciprofloxacin           1 resistant_list
#> 5                S           Ciprofloxacin           1 resistant_list
#> 6                R           Ciprofloxacin           1 resistant_list
prep_create_resistance_profile(ast, format = "full_pattern")
#> Creating resistance profiles (format: full_pattern)...
#> 
#> Top 10 resistance profiles:
#>                         resistance_profile n
#> 1 Ciprofloxacin:R; Imipenem:R; Meropenem:S 1
#> 2 Ciprofloxacin:R; Imipenem:S; Meropenem:S 1
#>   event_id   organism_normalized antibiotic_normalized antibiotic_class
#> 1       E1      Escherichia coli             Meropenem      Carbapenems
#> 2       E1      Escherichia coli              Imipenem      Carbapenems
#> 3       E1      Escherichia coli         Ciprofloxacin Fluoroquinolones
#> 4       E2 Klebsiella pneumoniae             Meropenem      Carbapenems
#> 5       E2 Klebsiella pneumoniae              Imipenem      Carbapenems
#> 6       E2 Klebsiella pneumoniae         Ciprofloxacin Fluoroquinolones
#>   antibiotic_value                       resistance_profile n_resistant
#> 1                S Ciprofloxacin:R; Imipenem:R; Meropenem:S           2
#> 2                R Ciprofloxacin:R; Imipenem:R; Meropenem:S           2
#> 3                R Ciprofloxacin:R; Imipenem:R; Meropenem:S           2
#> 4                S Ciprofloxacin:R; Imipenem:S; Meropenem:S           1
#> 5                S Ciprofloxacin:R; Imipenem:S; Meropenem:S           1
#> 6                R Ciprofloxacin:R; Imipenem:S; Meropenem:S           1
#>   n_tested profile_format
#> 1        3   full_pattern
#> 2        3   full_pattern
#> 3        3   full_pattern
#> 4        3   full_pattern
#> 5        3   full_pattern
#> 6        3   full_pattern
```
