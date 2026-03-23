# Calculate Deaths by Infectious Syndrome (D_L)

Computes D_L, the number of deaths due to each infectious syndrome L.

## Usage

``` r
calculate_syndrome_deaths(
  d_j = NULL,
  s_j = NULL,
  m_lj = NULL,
  pop_data = NULL,
  facility_data = NULL,
  cause_col = "cause_of_death",
  syndrome_col = "syndrome",
  syndrome = NULL,
  deaths_col = NULL,
  infection_flag_col = "is_infection_death",
  outcome_col = "final_outcome",
  death_value = "Died",
  patient_col = "patient_id",
  facility_col = NULL,
  facility_name = NULL,
  groupby_cols = NULL
)
```

## Arguments

- d_j:

  Data frame. Output of
  [`calculate_deaths_by_cause()`](https://saketlab.github.io/anumaan/reference/calculate_deaths_by_cause.md).
  Supply together with `s_j` and `m_lj` for Mode 1. Default `NULL`.

- s_j:

  Data frame. Output of
  [`calculate_infection_fraction()`](https://saketlab.github.io/anumaan/reference/calculate_infection_fraction.md).
  Default `NULL`.

- m_lj:

  Data frame. Output of
  [`calculate_syndrome_fraction()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_fraction.md).
  Default `NULL`.

- pop_data:

  Data frame. Raw population-level data for Mode 2. Default `NULL`.

- facility_data:

  Data frame. Facility-level data for Mode 3 fallback. Default `NULL`.

- cause_col:

  Character. Underlying cause column used in Modes 1 and 2. Default
  `"cause_of_death"`.

- syndrome_col:

  Character. Column containing syndrome labels in both population and
  facility data. Default `"syndrome"`.

- syndrome:

  Character. Optional. Name of a specific syndrome to filter results to
  (e.g., `"Bloodstream infection"`). `NULL` returns all syndromes.
  Default `NULL`.

- deaths_col:

  Character. Pre-aggregated deaths column in `pop_data` (Mode 2 only).
  `NULL` means each row is one death record. Default `NULL`.

- infection_flag_col:

  Character. Binary infection flag column in `pop_data` (Modes 1 and 2).
  Default `"is_infection_death"`.

- outcome_col:

  Character. Outcome column in `facility_data` (Mode 3). Default
  `"final_outcome"`.

- death_value:

  Character. Value in `outcome_col` that represents death. Accepts any
  string, e.g., `"Died"` or `"Death"`. Default `"Died"`.

- patient_col:

  Character. Unique patient identifier column in `facility_data`. Deaths
  are counted as distinct patients, not rows. Default `"patient_id"`.

- facility_col:

  Character. Column containing facility names in `facility_data`.
  Required when `facility_name` is specified. Default `NULL`.

- facility_name:

  Character. Name of a specific facility to restrict the analysis to.
  `NULL` uses all facilities combined. Default `NULL`.

- groupby_cols:

  Character vector. Additional stratification columns (e.g.,
  `c("Age_bin", "gender")`). Default `NULL`.

## Value

Data frame with columns: `syndrome_col`, any `groupby_cols`, `D_L`
(deaths by syndrome), `D_L_method`, `D_L_confidence`.

## Details

Three operating modes, evaluated in priority order:

1.  **Pre-computed components** (HIGH confidence): pass the outputs of
    [`calculate_deaths_by_cause()`](https://saketlab.github.io/anumaan/reference/calculate_deaths_by_cause.md),
    [`calculate_infection_fraction()`](https://saketlab.github.io/anumaan/reference/calculate_infection_fraction.md),
    and
    [`calculate_syndrome_fraction()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_fraction.md)
    via `d_j`, `s_j`, and `m_lj`. Computes: D_L = sum_J(D_J \* S_J \*
    M_LJ).

2.  **Raw population data** (HIGH confidence): pass `pop_data` and all
    three components are computed internally before combining. Computes:
    D_L = sum_J(D_J \* S_J \* M_LJ).

3.  **Facility fallback** (LOW confidence): when no population inputs
    are available, counts the number of **unique patients** who died for
    each syndrome directly from `facility_data`.

Use `syndrome` to restrict results to a single syndrome of interest. Use
`facility_col` and `facility_name` to restrict the facility fallback to
a specific site.

## References

Antimicrobial Resistance Collaborators. Global burden of bacterial
antimicrobial resistance in 2019. Lancet. 2022.

## Examples

``` r
if (FALSE) { # \dontrun{
# Mode 1: pass pre-computed components, all syndromes
d_j <- calculate_deaths_by_cause(pop_data = vr, cause_col = "icd10")
s_j <- calculate_infection_fraction(pop_data = vr, cause_col = "icd10")
m_lj <- calculate_syndrome_fraction(
  pop_data = vr, cause_col = "icd10", syndrome_col = "syndrome"
)
d_l <- calculate_syndrome_deaths(
  d_j = d_j, s_j = s_j, m_lj = m_lj,
  cause_col = "icd10", syndrome_col = "syndrome"
)

# Mode 2: raw population data, filter to one syndrome
d_l <- calculate_syndrome_deaths(
  pop_data     = vital_reg,
  cause_col    = "icd10_cause",
  syndrome_col = "infectious_syndrome",
  syndrome     = "Bloodstream infection"
)

# Mode 3: facility fallback, one facility, all syndromes
d_l <- calculate_syndrome_deaths(
  facility_data = amr_data,
  syndrome_col  = "specimen_normalized",
  outcome_col   = "final_outcome",
  death_value   = "Died",
  patient_col   = "patient_id",
  facility_col  = "location",
  facility_name = "Mumbai"
)

# Mode 3: facility fallback, all facilities, one syndrome, stratified by age
d_l <- calculate_syndrome_deaths(
  facility_data = amr_data,
  syndrome_col  = "specimen_normalized",
  syndrome      = "Blood",
  death_value   = "Died",
  groupby_cols  = c("Age_bin")
)
} # }
```
