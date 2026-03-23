# Load RR (Relative Risk) Reference Data

Loads the rr_list_gbd.xlsx file containing pathogen-drug combinations
with relative risk values for burden estimation.

## Usage

``` r
load_rr_reference()
```

## Value

Data frame with Pathogen, Drug, and RR columns

## Examples

``` r
if (FALSE) { # \dontrun{
rr_data <- load_rr_reference()
head(rr_data)
} # }
```
