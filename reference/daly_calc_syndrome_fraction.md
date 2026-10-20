# Calculate Infectious Syndrome Fraction (M_LJ)

Computes M_LJ, the fraction of infection-related deaths for underlying
cause J that are attributed to infectious syndrome L. This distributes
infection deaths across clinical syndrome categories (e.g., bloodstream
infection, pneumonia, urinary tract infection).

## Usage

``` r
daly_calc_syndrome_fraction(
  pop_data,
  cause_col = "cause_of_death",
  syndrome_col = "syndrome",
  infection_flag_col = "is_infection_death",
  groupby_cols = NULL
)
```

## Arguments

- pop_data:

  Data frame. Population-level mortality data with cause and syndrome
  columns. Required.

- cause_col:

  Character. Underlying cause of death column (cause J). Default
  `"cause_of_death"`.

- syndrome_col:

  Character. Infectious syndrome column (syndrome L), e.g.,
  `"infectious_syndrome"`. Default `"syndrome"`.

- infection_flag_col:

  Character. Binary column (TRUE/FALSE or 1/0) indicating infection
  involvement in `pop_data`. Default `"is_infection_death"`.

- groupby_cols:

  Character vector. Additional stratification columns. Default `NULL`.

## Value

Data frame with columns: `cause_col`, `syndrome_col`, any
`groupby_cols`, `infection_deaths_LJ` (deaths for syndrome L given cause
J), `infection_deaths_J` (total infection deaths for cause J), `M_LJ`
(syndrome fraction 0-1), `M_LJ_method`, `M_LJ_confidence`.

## Details

This function requires population-level data. If only facility data are
available, use
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
m_lj <- daly_calc_syndrome_fraction(vital_reg)
#> Computing M_LJ from population-level data...
#> M_LJ computed: 3 cause-syndrome combination(s) across 3 cause(s).
m_lj
#> # A tibble: 3 × 7
#>   cause_of_death syndrome infection_deaths_LJ infection_deaths_J  M_LJ
#>   <chr>          <chr>                  <int>              <int> <dbl>
#> 1 A41.9          BSI                        2                  2     1
#> 2 J18.9          LRI                        2                  2     1
#> 3 N39.0          UTI                        1                  1     1
#> # ℹ 2 more variables: M_LJ_method <chr>, M_LJ_confidence <chr>
```
