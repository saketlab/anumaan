# Compute YLD Attributable to Resistance (Direct, Profile-Specific)

Primary, direct implementation of the profile-specific attributable-YLD
equation:

## Usage

``` r
daly_calc_yld_attributable(
  profiles_with_los,
  P_Lk_prime_tbl,
  yld_ref = NULL,
  DW_sepsis = NULL,
  pathogen_col = "pathogen",
  plk_col = "P_Lk_prime",
  incidence_col = "N_NF_L",
  facility_col = NULL,
  facility_name = NULL,
  facility_state_map = NULL,
  state_col = "state",
  state_name = NULL,
  probability_col = "probability",
  dominant_class_col = "dominant_class",
  losr_col = "LOSR_years",
  loss_col = "LOSS_years"
)
```

## Arguments

- profiles_with_los:

  Named list (one entry per pathogen) from
  `daly_assign_rr_to_profiles(..., los_r_col = ..., los_s_col = ...)`.
  Each entry must have `probability_col`, `dominant_class_col`, and
  `losr_col` columns.

- P_Lk_prime_tbl:

  Data frame: the `P_Lk_prime` (pooled) or `facility_level` element from
  `daly_calc_pathogen_fraction_nonfatal()`. Must contain `pathogen_col`,
  `plk_col`, and `incidence_col`.

- yld_ref:

  Data frame with columns `location_name` and `DW_sepsis`. Ignored when
  `DW_sepsis` is supplied.

- DW_sepsis:

  Numeric scalar or `NULL`. Disability weight used directly for every
  row when supplied (recommended – see "Disability weight caveat").

- pathogen_col:

  Character. Default `"pathogen"`.

- plk_col:

  Character. P'\_Lk column in `P_Lk_prime_tbl`. Default `"P_Lk_prime"`.

- incidence_col:

  Character. I_L column in `P_Lk_prime_tbl`. Default `"N_NF_L"` (the
  column already produced by `daly_calc_pathogen_fraction_nonfatal()`).

- facility_col:

  Character or `NULL`. Facility identifier column in `P_Lk_prime_tbl`.
  Default `NULL`.

- facility_name:

  Character or `NULL`. Restrict to a single facility. Default `NULL`.

- facility_state_map:

  Data frame with `facility_col` and `state_col`, mapping each facility
  to a state. Used only when `DW_sepsis` is not supplied and
  `facility_col` is present.

- state_col:

  Character. Default `"state"`.

- state_name:

  Character or `NULL`. State for the pooled/no-facility DW lookup.
  `NULL` uses the "India" row. Ignored when `DW_sepsis` is supplied.

- probability_col:

  Character. R'\_k_delta column in `profiles_with_los`. Default
  `"probability"`.

- dominant_class_col:

  Character. d\*(delta) column. Default `"dominant_class"`.

- losr_col:

  Character. LOSR_k,d\*(delta) column (years). Default `"LOSR_years"`.

- loss_col:

  Character. LOSS_k,d\*(delta) column (years) in `profiles_with_los`.
  Default `"LOSS_years"`.

## Value

`P_Lk_prime_tbl` augmented with `I_L`, `DW_L`, `profile_term` (=
sum_delta R'\_k_delta \* (LOSR_k_delta - LOSS_k_delta)), and
`YLD_attributable`.

## Details

YLD_attrib_k = I_L \* P'\_Lk \* DW_L \* sum\_{delta in D_k}
\[R'\_k_delta \* (LOSR_k,d\*(delta) - LOSS_k,d\*(delta))\]

The excess `(LOSR - LOSS)` term makes this a genuine counterfactual: how
much MORE disability burden exists because these infections were
resistant, versus if they had been susceptible to the same (dominant)
class. D_k excludes the all-susceptible profile from the sum entirely
(rather than including a zero `LOSR - LOSS` term for it) – see
[`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md)
for why.

Obtains I_L, P'\_Lk, DW_L, R'\_k_delta, d\*(delta), LOSR, and LOSS the
same way as
[`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md)
– see that function's documentation for full component-by-component
sourcing and the disability-weight caveat. Does **not** require a
baseline YLD or a PAF_LOS computed first; unlike the pre-existing
multiplicative pathway
([`daly_calc_paf_los()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_los.md)
multiplied against a pooled baseline YLD), this is the direct equation.
[`daly_calc_paf_los()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_los.md)
remains available standalone for diagnostics / comparison against this
direct calculation.
