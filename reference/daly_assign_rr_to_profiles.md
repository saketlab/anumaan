# Assign Per-Class LOS RR to Resistance Profiles (Max Rule)

For each resistance profile delta, determines the profile-level relative
LOS using the GBD max rule: RR_kd_LOS = max over c in C_R(d) of
RR_kc_LOS \[if C_R(d) non-empty\] = 1 \[if d = all-susceptible\] where
C_R(d) is the set of classes the profile is resistant to. The class
attaining the maximum is the profile's dominant class d\*. The CI
reported for each profile is that of its dominant class.

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
  los_unit = c("years", "days"),
  profile_col = "profile_delta",
  profile_probability_col = "R_mean",
  profile_hospital_col = "group",
  profile_classes_col = "classes",
  profile_class_order = NULL
)
```

## Arguments

- profiles_output:

  Legacy list from
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
  or a data frame of profile strings (see "Profile input formats").

- rr_table:

  The list returned by
  [`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md),
  or a data frame with columns `pathogen_col`, `class_col`, `rr_col` and
  optionally `CI_lower` / `CI_upper` (e.g. from
  [`daly_fit_los_rr_distribution()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr_distribution.md)).

- pathogen_col:

  Character. Pathogen column in `rr_table` (and in string profiles).
  Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- rr_col:

  Character. Default `"RR_LOS"`.

- fallback_rr:

  Numeric. RR for resistant classes with no match. Default `1` (no
  attributable effect).

- los_r_col, los_s_col:

  Character or `NULL`. Columns in a data-frame `rr_table` with absolute
  LOS for resistant / susceptible patients (e.g. `"median_LOS_R"` /
  `"median_LOS_S"`). Ignored when `rr_table` is a
  [`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md)
  list.

- los_unit:

  Character. Unit of `los_r_col` / `los_s_col`: `"years"` or `"days"`.
  Default `"years"`.

- profile_col:

  Character. Profile-string column. Default `"profile_delta"`.

- profile_probability_col:

  Character. Profile probability column in string profiles. Default
  `"R_mean"`. Copied to `probability`.

- profile_hospital_col:

  Character or `NULL`. Hospital column in string profiles. Default
  `"group"`; `NULL` or absent means pooled profiles. Copied to
  `hospital`.

- profile_classes_col:

  Character. Column of string profiles holding the class order as one
  `;`-separated string. Default `"classes"`.

- profile_class_order:

  `NULL` (read the order from `profile_classes_col`), a character vector
  (the class of each letter, in order) or a named list of such vectors,
  one per pathogen.

## Value

Named list (one entry per pathogen): the profiles data frame augmented
with `RR_LOS_profile`, `dominant_class`, (if available)
`CI_lower_profile` / `CI_upper_profile`, `LOSR_years` / `LOSS_years`,
and for string profiles `probability`, `hospital` and `los_level`.

## Details

**Ties**: when several resistant classes share the maximum RR, the first
of them in the class order wins (`profile_class_order` for string
profiles, the `classes` element for legacy profiles). A message reports
how many profiles had a tie.

If `rr_table` was fit with a `syndrome_name` filter, its RR values are
syndrome-specific but get applied here to all profiles of a pathogen –
i.e. it assumes syndrome-invariant LOS prolongation.

## Profile input formats

- Legacy:

  The named list from
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
  (`$profiles` with 0/1 class columns and `$classes`).

- Profile strings:

  A data frame with one row per (hospital x) pathogen x profile, where
  the profile is a string of R/S letters such as `"RRSR"` – e.g. the
  output of
  [`profile_aggregate_for_daly()`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)
  (`group`, `pathogen`, `profile_delta`, `R_mean`), or a burden table
  with `hospital`, `pathogen`, `profile`, `mean` and `classes`. When the
  default column names are absent, those of the burden table are used.

## Class of each letter

The class order is read from `profile_classes_col` (e.g.
`"Fluoroquinolones;Aminoglycosides;Carbapenems"`), which must hold a
single value per pathogen. `profile_class_order` can be given instead,
or to override it (a warning is raised if they differ). Class names are
matched to the LOS results ignoring case and punctuation, so
`"Third_generation_cephalosporins"` matches
`"Third-generation-cephalosporins"`.

## LOS for each profile

When `rr_table` is the list returned by
[`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md),
the profile takes the adjusted LOS of its dominant class: from
`rr_table$hospital` for the profile's own hospital when the profiles
carry a hospital column (`los_level = "hospital"`), otherwise the
all-hospital values in `rr_table$rr` (`los_level = "network"`). Days are
converted to years. With a data-frame `rr_table`, absolute LOS is taken
from `los_r_col` / `los_s_col` as before. `LOSR_years` / `LOSS_years`
are `NA` for the all-susceptible profile, which has no dominant class.

## Two options for profile-level LOS

This function is **option 2**: per-class LOS estimates mapped to
profiles by the max rule. **Option 1** is
[`daly_fit_los_distn_by_profile()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_distn_by_profile.md)
(Section 5.3), which estimates LOS directly for each patient-level
resistance profile against the fully susceptible profile.
