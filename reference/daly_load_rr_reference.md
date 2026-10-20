# Load RR (Relative Risk) Reference Data

Loads the bundled rr_list_gbd.csv file containing pathogen-drug
combinations with relative risk values for burden estimation.

## Usage

``` r
daly_load_rr_reference()
```

## Value

Data frame with Pathogen, Drug, and RR columns

## Examples

``` r
rr_data <- daly_load_rr_reference()
#> [v] Loaded 86 pathogen-drug RR combinations
head(rr_data)
#>                  pathogen                                                  drug
#> 1 acinetobacter baumannii Anti-pseudomonal penicillin/Beta-lactamase inhibitors
#> 2 acinetobacter baumannii                 Beta-lactam/Beta-lactamase inhibitors
#> 3 acinetobacter baumannii                                            Carbapenem
#> 4 acinetobacter baumannii                      Fourth-generation cephalosporins
#> 5 acinetobacter baumannii                       Third-generation cephalosporins
#> 6 acinetobacter baumannii                                       Aminoglycosides
#>   sample_size mean_rr lower_bound upper_bound
#> 1         948    1.31        1.12        1.52
#> 2        1555    1.27        1.11        1.44
#> 3        3232    1.42        1.27        1.58
#> 4        1439    1.31        1.14        1.51
#> 5        2055    1.35        1.13        1.62
#> 6        2066    1.10        0.97        1.25
```
