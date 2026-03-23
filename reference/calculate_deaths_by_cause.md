# Calculate Deaths by Underlying Cause (D_J)

Computes D_J, the number of deaths for each underlying cause J from
population-level vital registration or mortality data. Optionally
stratified by grouping variables such as age group, sex, year, or
location.

## Usage

``` r
calculate_deaths_by_cause(
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
[`calculate_syndrome_deaths()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_deaths.md)
with `facility_data` instead, which directly counts deaths by syndrome.

## References

Antimicrobial Resistance Collaborators. Global burden of bacterial
antimicrobial resistance in 2019. Lancet. 2022.

## Examples

``` r
if (FALSE) { # \dontrun{
# One row per death record
d_j <- calculate_deaths_by_cause(
  pop_data  = vital_reg,
  cause_col = "icd10_cause"
)

# Pre-aggregated counts, stratified by age and sex
d_j <- calculate_deaths_by_cause(
  pop_data     = vital_reg,
  cause_col    = "icd10_cause",
  deaths_col   = "n_deaths",
  groupby_cols = c("Age_bin", "gender")
)
} # }
```
