# Calculate Deaths by Underlying Cause (D_J)

Computes D_J, the number of deaths for each underlying cause J from
population-level vital registration or mortality data. Optionally
stratified by grouping variables such as age group, sex, year, or
location.

## Usage

``` r
daly_calc_deaths_by_cause(
  pop_data,
  cause_col = "cause_of_death",
  deaths_col = NULL,
  groupby_cols = NULL
)
```

## Arguments

- pop_data:

  Data frame. Population-level vital registration or mortality data.
  Required – this function does not accept facility-level data.

- cause_col:

  Character. Column containing the underlying cause of death (cause J),
  e.g. an ICD-10 code or cause name. Default `"cause_of_death"`.

- deaths_col:

  Character. Column with pre-aggregated death counts. Set to `NULL` if
  each row represents one individual death record. Default `NULL`.

- groupby_cols:

  Character vector. Additional stratification columns (e.g.,
  `c("Age_bin", "gender", "year")`). Default `NULL`.

## Value

Data frame with columns: `cause_col`, any `groupby_cols`, `D_J` (death
count), `D_J_method` (`"population"`), `D_J_confidence` (`"high"`).

## Details

This function requires population-level data. If only facility-level
data are available, use
[`daly_calc_deaths_by_syndrome()`](https://saketlab.github.io/anumaan/reference/daly_calc_deaths_by_syndrome.md)
with `facility_data` instead, which directly counts deaths by syndrome.

## References

Antimicrobial Resistance Collaborators. Global burden of bacterial
antimicrobial resistance in 2019. Lancet. 2022.

## Examples

``` r
vital_reg <- data.frame(
  cause_of_death     = c("A41.9", "A41.9", "J18.9", "J18.9", "I21.9", "N39.0"),
  is_infection_death = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  syndrome           = c("BSI", "BSI", "LRI", "LRI", NA, "UTI")
)
d_j <- daly_calc_deaths_by_cause(vital_reg)
#> Computing D_J from population-level data...
#> D_J computed: 6 total deaths across 4 cause(s).
d_j
#> # A tibble: 4 × 4
#>   cause_of_death   D_J D_J_method D_J_confidence
#>   <chr>          <int> <chr>      <chr>         
#> 1 A41.9              2 population high          
#> 2 I21.9              1 population high          
#> 3 J18.9              2 population high          
#> 4 N39.0              1 population high          
```
