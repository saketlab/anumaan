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
  fallback_rr = 1
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

## Value

Named list (one entry per pathogen): original profiles data frame
augmented with RR_LOS_profile, dominant_class, and (if available)
CI_lower_profile / CI_upper_profile.

## Details

If `rr_table` was fit with a `syndrome_name` filter, its RR values are
syndrome-specific but get applied here to all profiles of a pathogen –
i.e. it assumes syndrome-invariant LOS prolongation. Refit with
`syndrome_name = NULL` for a RR pooled across syndromes.
