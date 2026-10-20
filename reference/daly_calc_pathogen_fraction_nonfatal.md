# Calculate non-fatal pathogen distribution (P'\_{Lk})

Computes the non-fatal pathogen distribution for an infectious syndrome
(L): the fractional contribution of each pathogen (k) among patients who
were discharged alive. Its `P_Lk_prime` element is the input expected by
[`daly_calc_incidence_from_cfr()`](https://saketlab.github.io/anumaan/reference/daly_calc_incidence_from_cfr.md)
and `daly_calc_yld()`.

## Usage

``` r
daly_calc_pathogen_fraction_nonfatal(
  data,
  syndrome_col,
  syndrome_name,
  specimen_col = NULL,
  specimen_name = NULL,
  polymicrobial_col,
  patient_col,
  pathogen_col,
  outcome_col,
  discharged_value = "Discharged",
  glass_ref = NULL,
  facility_col = NULL,
  facility_name = NULL,
  pathogen_name = NULL
)
```

## Arguments

- data:

  Data frame of facility-level microbiology records.

- syndrome_col:

  Character. Column containing infectious syndrome labels (L).

- syndrome_name:

  Character. Syndrome to analyse.

- specimen_col:

  Character or NULL. Column containing specimen type. If `NULL`, no
  specimen filter is applied.

- specimen_name:

  Character or NULL. Specimen to restrict to (e.g., `"Blood"`). Required
  when `specimen_col` is provided.

- polymicrobial_col:

  Character. Column flagging polymicrobial patients (1 = polymicrobial,
  0 = monomicrobial).

- patient_col:

  Character. Unique patient identifier column.

- pathogen_col:

  Character. Pathogen (organism) column (k).

- outcome_col:

  Character. Final patient outcome column.

- discharged_value:

  Character. Value of `outcome_col` indicating a non-fatal (discharged)
  outcome. Default `"Discharged"`.

- glass_ref:

  Character vector of valid pathogen names, or a data frame with columns
  `specimen` and `pathogen`. Applied to polymicrobial patients only.
  When `specimen_name` is `NULL` and `glass_ref` is a data frame, all
  pathogens in the reference are used regardless of specimen. `NULL`
  skips GLASS filtering.

- facility_col:

  Character or NULL. Facility identifier column. When provided without
  `facility_name`, returns both facility-level and pooled P_LK.

- facility_name:

  Character or NULL. Restricts to a single facility.

- pathogen_name:

  Character vector or NULL. Filter to specific pathogen(s).

## Value

A list with `P_Lk_prime` (pooled distribution, with column `P_Lk_prime`)
and, when `facility_col` is supplied without `facility_name`,
`facility_level`.

## Details

Weighting mirrors
[`daly_calc_pathogen_fraction_fatal()`](https://saketlab.github.io/anumaan/reference/daly_calc_pathogen_fraction_fatal.md):
each patient contributes total weight 1, split equally across their
valid pathogens, and GLASS filtering applies to polymicrobial patients
only.

## Examples

``` r
fac <- data.frame(
  patient_id = sprintf("P%02d", 1:6), syndrome = "BSI", poly = 0,
  organism = c("E. coli", "E. coli", "E. coli", "K. pneumoniae", "K. pneumoniae", "S. aureus"),
  outcome = c("Discharged", "Discharged", "Death", "Discharged", "Discharged", "Death")
)
daly_calc_pathogen_fraction_nonfatal(fac, "syndrome", "BSI",
  polymicrobial_col = "poly", patient_col = "patient_id",
  pathogen_col = "organism", outcome_col = "outcome"
)
#> P'LK computed: 2 pathogen(s), N^NF_L = 4 patient(s).
#> $P_Lk_prime
#> # A tibble: 2 × 4
#>   organism      N_NF_LK N_NF_L P_Lk_prime
#>   <chr>           <dbl>  <int>      <dbl>
#> 1 E. coli             2      4        0.5
#> 2 K. pneumoniae       2      4        0.5
#> 
```
