# Calculate YLD per pathogen

Computes Years Lived with Disability (YLD) attributable to each pathogen
K within syndrome L:

## Usage

``` r
calculate_YLD(
  incidence_data,
  P_Lk_prime_tbl,
  yld_ref = NULL,
  DW_sepsis = NULL,
  avg_los_years = NULL,
  state_name = NULL,
  facility_col = NULL,
  facility_name = NULL,
  facility_state_map = NULL,
  state_col = "state",
  pathogen_col = "pathogen",
  pathogen_name = NULL,
  plk_col = "P_Lk_prime",
  incidence_col = "n_cases"
)
```

## Arguments

- incidence_data:

  Numeric scalar (pooled) or data frame with `facility_col` +
  `incidence_col` columns (facility-level).

- P_Lk_prime_tbl:

  Data frame from
  [`calculate_P_Lk_prime()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk_prime.md).
  Use the `P_Lk_prime` element for pooled mode or the `facility_level`
  element for facility-level mode.

- yld_ref:

  Data frame with columns `location_name` and `DW_sepsis`. Loaded from
  `inst/extdata/Proxy_YLD_per_case.xlsx`. Optional if `DW_sepsis` is
  provided.

- DW_sepsis:

  Numeric scalar or NULL. If provided, this value is used directly for
  all rows and `yld_ref` is ignored.

- avg_los_years:

  Numeric scalar or NULL. Overall average length of stay in years across
  all patients (resistant and susceptible combined), used to convert DW
  to a duration-weighted value:
  `effective_DW = DW_sepsis * avg_los_years`. If NULL, DW_sepsis is used
  as-is (caller is responsible for duration weighting).

- state_name:

  Character or NULL. State for YLD weight in pooled mode. NULL uses the
  India row. Ignored in facility-level mode when `DW_sepsis` is
  provided.

- facility_col:

  Character or NULL. Facility identifier column.

- facility_name:

  Character or NULL. Restrict to one facility only.

- facility_state_map:

  Data frame with `facility_col` and `state_col` columns, mapping each
  facility to its state. Required in facility-level mode only when
  `DW_sepsis` is not provided.

- state_col:

  Character. Column in `facility_state_map` containing state names that
  match `yld_ref$location_name`. Default `"state"`.

- pathogen_col:

  Character. Pathogen column in `P_Lk_prime_tbl`. Default `"pathogen"`.

- pathogen_name:

  Character vector or NULL. If provided, restricts output to those
  pathogen(s) only.

- plk_col:

  Character. P'LK column in `P_Lk_prime_tbl`. Default `"P_Lk_prime"`.

- incidence_col:

  Character. Incidence column when `incidence_data` is a data frame.
  Default `"n_cases"`.

## Value

Data frame with columns: `pathogen_col`, `P_Lk_prime`, `incidence_L`,
`DW_sepsis`, `YLD`, and `facility_col` / `state_col` when in
facility-level mode.

## Details

\$\$YLD_K = \text{Incidence}\_L \times P'\_{LK} \times
\text{DW\\sepsis}\$\$

The YLD weight per incident case is drawn from a GBD-derived reference
table stratified by Indian state (`yld_ref`) when `DW_sepsis` is not
provided. If `DW_sepsis` is provided directly, that scalar is used for
all rows and `yld_ref` is ignored.

**Facility-level mode** (when `facility_col` is supplied and
`facility_name` is NULL): `P_Lk_prime_tbl` must contain a `facility_col`
column (i.e., the `facility_level` element from
[`calculate_P_Lk_prime()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk_prime.md)).
`incidence_data` must be a data frame with `facility_col` and an
incidence count column. A `facility_state_map` data frame linking each
facility to its state is required only when `DW_sepsis` is not supplied.

**Pooled / no-facility mode**: `incidence_data` is a single numeric
scalar and `P_Lk_prime_tbl` has no facility column. `state_name` selects
the YLD weight; defaults to "India" if NULL.
