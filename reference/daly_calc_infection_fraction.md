# Calculate Infection Fraction of Deaths by Cause (S_J)

Computes S_J, the fraction of deaths for underlying cause J that were
related to infection. Values range from 0 (no infection involvement) to
1 (all deaths for cause J involved infection).

## Usage

``` r
daly_calc_infection_fraction(
  pop_data,
  cause_col = "cause_of_death",
  infection_flag_col = "is_infection_death",
  groupby_cols = NULL
)
```

## Arguments

- pop_data:

  Data frame. Population-level mortality data. Required.

- cause_col:

  Character. Underlying cause of death column (cause J). Default
  `"cause_of_death"`.

- infection_flag_col:

  Character. Binary column (TRUE/FALSE or 1/0) in `pop_data` indicating
  whether the death involved infection. Default `"is_infection_death"`.

- groupby_cols:

  Character vector. Stratification columns present in `pop_data`.
  Default `NULL`.

## Value

Data frame with columns: `cause_col`, any `groupby_cols`, `D_J` (total
deaths), `infection_deaths` (infection-related death count), `S_J`
(infection fraction 0-1), `S_J_method`, `S_J_confidence`.

## Details

This function requires population-level data where each death record
carries a binary flag indicating whether infection was involved. If only
facility data are available, use
[`daly_calc_deaths_by_syndrome()`](https://saketlab.github.io/anumaan/reference/daly_calc_deaths_by_syndrome.md)
with `facility_data` instead.

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
s_j <- daly_calc_infection_fraction(vital_reg)
#> Computing S_J from population-level data...
#> S_J computed: mean infection fraction = 0.750 across 4 cause(s).
s_j
#> # A tibble: 4 × 6
#>   cause_of_death   D_J infection_deaths   S_J S_J_method S_J_confidence
#>   <chr>          <int>            <dbl> <dbl> <chr>      <chr>         
#> 1 A41.9              2                2     1 population high          
#> 2 I21.9              1                0     0 population high          
#> 3 J18.9              2                2     1 population high          
#> 4 N39.0              1                1     1 population high          
```
