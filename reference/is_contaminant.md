# Check if Organism is a Contaminant

Checks if a given organism name matches any known contaminant for a
specific syndrome or specimen type using flexible pattern matching.

## Usage

``` r
is_contaminant(organism_name, syndrome = NULL, specimen_type = NULL)
```

## Arguments

- organism_name:

  Character vector. Organism name(s) to check.

- syndrome:

  Character. Optional syndrome name to filter contaminants.

- specimen_type:

  Character. Optional specimen type to filter contaminants.

## Value

Logical vector indicating if each organism is a contaminant

## Examples

``` r
# Check if organism is a blood culture contaminant
is_contaminant("Staph epidermidis", syndrome = "Bloodstream infections")
#> [1] TRUE
is_contaminant("E. coli", syndrome = "Bloodstream infections")
#> [1] FALSE

# Check multiple organisms
is_contaminant(c("Staph", "Klebsiella"), specimen_type = "Blood culture")
#> [1] FALSE FALSE
```
