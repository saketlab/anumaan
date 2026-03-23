# Calculate Infection Fraction of Deaths by Cause (S_J)

Computes S_J, the fraction of deaths for underlying cause J that were
related to infection. Values range from 0 (no infection involvement) to
1 (all deaths for cause J involved infection).

## Usage

``` r
calculate_infection_fraction(
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
[`calculate_syndrome_deaths()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_deaths.md)
with `facility_data` instead.

## References

Antimicrobial Resistance Collaborators. Global burden of bacterial
antimicrobial resistance in 2019. Lancet. 2022.

## Examples

``` r
if (FALSE) { # \dontrun{
s_j <- calculate_infection_fraction(
  pop_data           = vital_reg,
  cause_col          = "icd10_cause",
  infection_flag_col = "is_infectious",
  groupby_cols       = c("Age_bin", "gender")
)
} # }
```
