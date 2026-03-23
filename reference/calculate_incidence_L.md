# Calculate syndrome incidence from deaths, CFR, and CR_L (formula-based)

Estimates the number of incident cases of syndrome L using:

## Usage

``` r
calculate_incidence_L(
  deaths_L,
  cfr_lk_tbl,
  P_Lk_prime_tbl,
  CR_L = 1,
  pathogen_col = "pathogen",
  cfr_col = "CFR_LK",
  plk_col = "P_Lk_prime",
  facility_col = NULL,
  deaths_col = "deaths"
)
```

## Arguments

- deaths_L:

  Numeric scalar (pooled mode) or data frame with a `facility_col`
  column and a `deaths_col` column (facility-level mode).

- cfr_lk_tbl:

  Data frame with at minimum columns `pathogen_col` and `cfr_col`.
  Typically the `cfr_table` element from
  [`calculate_cfr_lk()`](https://saketlab.github.io/anumaan/reference/calculate_cfr_lk.md).

- P_Lk_prime_tbl:

  Data frame with `pathogen_col` and `plk_col`. Use the `P_Lk_prime`
  (pooled) element from `calculate_P_Lk_prime_BSI()` or
  [`calculate_P_Lk_prime()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk_prime.md).

- CR_L:

  Numeric scalar. CFR adjustment factor from
  [`calculate_CR_L()`](https://saketlab.github.io/anumaan/reference/calculate_CR_L.md).
  Default `1` (no adjustment).

- pathogen_col:

  Character. Pathogen column in both tables. Default `"pathogen"`.

- cfr_col:

  Character. CFR column in `cfr_lk_tbl`. Default `"CFR_LK"`.

- plk_col:

  Character. P'LK column in `P_Lk_prime_tbl`. Default `"P_Lk_prime"`.

- facility_col:

  Character or NULL. Facility identifier. When provided, `cfr_lk_tbl`
  and `P_Lk_prime_tbl` must each contain `facility_col`, and `deaths_L`
  must be a data frame with `facility_col` + `deaths_col`.

- deaths_col:

  Character. Column in `deaths_L` data frame containing death counts.
  Default `"deaths"`. Ignored when `deaths_L` is a scalar.

## Value

Data frame with columns `deaths`, `CFR_L`, `CR_L`, `I_L` (incident
cases), and `facility_col` when applicable.

## Details

\$\$I_L = \frac{D_L}{\text{CFR}\_L \times \text{CR}\_L}\$\$

where the syndrome-level CFR is the pathogen-weighted average:

\$\$\text{CFR}\_L = \sum_k P'\_{Lk} \times \text{CFR}\_{Lk}\$\$

Use this when you have population- or facility-level death counts and
want to back-calculate incidence. For a direct patient count from
facility data, use
[`count_incident_cases()`](https://saketlab.github.io/anumaan/reference/count_incident_cases.md)
instead.
