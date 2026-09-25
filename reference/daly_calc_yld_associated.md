# Compute YLD Associated with Resistance (Direct, Profile-Specific)

Primary, direct implementation of the profile-specific associated-YLD
equation:

## Usage

``` r
daly_calc_yld_associated(
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
  losr_col = "LOSR_years"
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

## Value

`P_Lk_prime_tbl` augmented with `I_L`, `DW_L`, `profile_term` (=
sum_delta R'\_k_delta \* LOSR_k_delta), and `YLD_associated`.

## Details

YLD_assoc_k = I_L \* P'\_Lk \* DW_L \* sum\_{delta in D_k} \[R'\_k_delta
\* LOSR_k,d\*(delta)\]

where D_k is the set of RESISTANT resistance profiles for pathogen k
(the all-susceptible profile has no dominant class and is excluded – see
"Details"). d\*(delta) is the dominant drug class of profile delta (GBD
max rule, assigned by
[`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)).

This does **not** require computing a baseline/pooled YLD first: unlike
the pre-existing multiplicative pathway (a pooled baseline YLD
multiplied by an associated-burden fraction), this is the direct
equation, evaluated per pathogen (and per facility, if `P_Lk_prime_tbl`
is facility-level) from its own inputs.

**Component sourcing**:

- **I_L** and **P'\_Lk**: read directly from `P_Lk_prime_tbl` – the
  `P_Lk_prime` (pooled) or `facility_level` element already returned by
  `daly_calc_pathogen_fraction_nonfatal()`, which computes both the
  non-fatal incidence count and the pathogen fraction from the same
  observed cohort. Incidence is never back-calculated from deaths, CFR,
  or LOS.

- **DW_L**: resolved exactly as in the package's pre-existing YLD weight
  mechanism – either the `DW_sepsis` scalar, or a lookup in `yld_ref`
  (loaded from `inst/extdata/Proxy_YLD_per_case.csv`) by state/facility.
  See "Disability weight caveat" below.

- **R'\_k_delta**, **d\*(delta)**, and **LOSR**: read from
  `profiles_with_los`, the list returned by
  [`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
  when called with `los_r_col`/`los_s_col` set.

**Why D_k excludes the all-susceptible profile**: this matches the
pre-existing convention elsewhere in this package that "associated"
burden is a partition over resistant profiles only (the removed
`daly_calc_fraction_associated_yld()`'s `Fraction_k` summed over
resistant profiles only;
[`daly_calc_paf_los()`](https://saketlab.github.io/anumaan/reference/daly_calc_paf_los.md)'s
all-susceptible term already contributes exactly zero by construction,
since RR=1 there). `YLD_associated` therefore represents the burden
occurring specifically among pathogen-k patients with a resistant
profile, not the total burden across all pathogen-k patients.

**Disability weight caveat**: `Proxy_YLD_per_case.csv`'s `DW_sepsis`
column is empirically a duration-inclusive "YLD per case" proxy
(`proxy_yld_per_case == yld_days_per_case / 365`), not a pure disability
weight. Because this function already multiplies by profile-specific
`LOSR_years` separately, sourcing `DW_L` from `yld_ref` risks
double-counting duration; a warning is emitted when `DW_sepsis` is not
supplied directly. Prefer supplying a pure disability weight via the
`DW_sepsis` scalar for this pathway.
