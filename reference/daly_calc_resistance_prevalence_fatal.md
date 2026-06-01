# Calculate Fatal Resistance Prevalence (R_k)

Computes the profile-weighted expected proportion of deaths attributable
to resistance (R_k) and the expected odds ratio of death (E\[OR_death\])
for each pathogen.

## Usage

``` r
daly_calc_resistance_prevalence_fatal(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
)
```

## Arguments

- profiles_with_rr:

  Named list returned by
  [`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
  or
  [`daly_filter_profiles_to_rr_classes()`](https://saketlab.github.io/anumaan/reference/daly_filter_profiles_to_rr_classes.md).

- probability_col:

  Character. Profile probability column. Default `"probability"`.

- rr_profile_col:

  Character. Profile-level OR column. Default `"RR_LOS_profile"`.

## Value

Named list (one entry per pathogen) each containing `per_profile` data
frame, `R_k`, and `E_OR_k`.
