# Calculate Infectious Syndrome Fraction (M_LJ)

Computes M_LJ, the fraction of infection-related deaths for underlying
cause J that are attributed to infectious syndrome L. This distributes
infection deaths across clinical syndrome categories (e.g., bloodstream
infection, pneumonia, urinary tract infection).

## Usage

``` r
calculate_syndrome_fraction(
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
[`calculate_syndrome_deaths()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_deaths.md)
with `facility_data` instead.

## References

Antimicrobial Resistance Collaborators. Global burden of bacterial
antimicrobial resistance in 2019. Lancet. 2022.

## Examples

``` r
if (FALSE) { # \dontrun{
m_lj <- calculate_syndrome_fraction(
  pop_data           = vital_reg,
  cause_col          = "icd10_cause",
  syndrome_col       = "infectious_syndrome",
  infection_flag_col = "is_infectious",
  groupby_cols       = c("Age_bin")
)
} # }
```
