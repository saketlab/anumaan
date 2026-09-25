# Assign Per-Class LOS RR to Resistance Profiles (Max Rule)

For each resistance profile delta (from compute_resistance_profiles()),
determines the profile-level RR_kd_LOS using the GBD max rule: RR_kd_LOS
= max over c in C_R(d) of RR_kc_LOS \[if C_R(d) non-empty\] = 1 \[if d =
all-susceptible\] where C_R(d) = {c : d_c = 1}. The CI reported for each
profile is that of its dominant (max-RR) class.

## Usage

``` r
daly_assign_rr_to_profiles(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1,
  los_r_col = NULL,
  los_s_col = NULL,
  los_unit = c("years", "days")
)
```

## Arguments

- profiles_output:

  Named list from compute_resistance_profiles().

- rr_table:

  Data frame from daly_fit_los_rr() or daly_fit_los_rr_distribution().
  Must have columns pathogen_col, class_col, rr_col, and optionally
  CI_lower / CI_upper.

- pathogen_col:

  Character. Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- rr_col:

  Character. Default `"RR_LOS"`.

- fallback_rr:

  Numeric. RR for resistant classes with no match. Default `1` (no
  attributable effect).

- los_r_col:

  Character or `NULL`. Column in `rr_table` with absolute mean LOS for
  resistant patients (e.g. `"mean_LOS_R"`). When supplied together with
  `los_s_col`, `LOSR_years` / `LOSS_years` are added to the output.
  Default `NULL` (no absolute LOS propagated; behaviour identical to
  before this parameter existed).

- los_s_col:

  Character or `NULL`. Column in `rr_table` with absolute mean LOS for
  susceptible patients (e.g. `"mean_LOS_S"`). Required when `los_r_col`
  is supplied.

- los_unit:

  Character. Unit of `los_r_col` / `los_s_col` in `rr_table`: `"years"`
  (used as-is) or `"days"` (divided by 365, matching the day-to-year
  conversion already used in
  [`daly_fit_los_rr_distribution()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr_distribution.md)).
  Default `"years"`. Ignored when `los_r_col` is `NULL`.

## Value

Named list (one entry per pathogen): original profiles data frame
augmented with RR_LOS_profile, dominant_class, (if available)
CI_lower_profile / CI_upper_profile, and (when `los_r_col` / `los_s_col`
are supplied) LOSR_years / LOSS_years.

## Details

If `rr_table` was fit with a `syndrome_name` filter, its RR values are
syndrome-specific but get applied here to all profiles of a pathogen –
i.e. it assumes syndrome-invariant LOS prolongation. Refit with
`syndrome_name = NULL` for a RR pooled across syndromes.

Optionally, when `los_r_col` and `los_s_col` are supplied, also
propagates the **absolute** mean LOS for the profile's dominant class
(resistant and susceptible) as `LOSR_years` / `LOSS_years`. This does
not change the RR assignment above; it is an additive lookup against the
dominant class already selected by the max rule, using whichever
absolute-LOS columns are present in `rr_table` (e.g. `mean_LOS_R` /
`mean_LOS_S` from
[`daly_fit_los_rr_distribution()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr_distribution.md),
or `mean_los_resistant` / `mean_los_susceptible` from
[`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md)).
`LOSR_years`/`LOSS_years` are set to `NA` for the all-susceptible
profile, which has no dominant class.
