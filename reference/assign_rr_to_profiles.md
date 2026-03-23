# Assign Per-Class LOS RR to Resistance Profiles (Max Rule)

For each resistance profile delta (from compute_resistance_profiles()),
determines the profile-level RR_kd_LOS using the GBD max rule: RR_kd_LOS
= max over c in C_R(d) of RR_kc_LOS \[if C_R(d) non-empty\] = 1 \[if d =
all-susceptible\] where C_R(d) = {c : d_c = 1}. The CI reported for each
profile is that of its dominant (max-RR) class.

For each resistance profile delta, determines the profile-level RR using
the GBD max rule: RR_kd = max over resistant classes of RR_kc. The CI
reported for each profile is that of its dominant (max-RR) class. Also
used for mortality ORs when called with `rr_col = "OR_death"`.

## Usage

``` r
assign_rr_to_profiles(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1
)

assign_rr_to_profiles(
  profiles_output,
  rr_table,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rr_col = "RR_LOS",
  fallback_rr = 1
)
```

## Arguments

- profiles_output:

  Named list from
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).

- rr_table:

  Data frame from
  [`fit_mortality_rr_logistic()`](https://saketlab.github.io/anumaan/reference/fit_mortality_rr_logistic.md)
  or a LOS RR table. Must have `pathogen_col`, `class_col`, `rr_col`,
  and optionally `CI_lower` / `CI_upper`.

- pathogen_col:

  Character. Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- rr_col:

  Character. Default `"RR_LOS"`. Use `"OR_death"` for mortality.

- fallback_rr:

  Numeric. RR for classes with no match. Default `1`.

## Value

Named list (one entry per pathogen): original profiles data frame
augmented with RR_LOS_profile, dominant_class, and (if available)
CI_lower_profile / CI_upper_profile.

Named list (one per pathogen): profiles data frame augmented with
`RR_LOS_profile`, `dominant_class`, and CI columns.
