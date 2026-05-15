# Compute PAF for length of stay per resistance profile

Computes the population attributable fraction (PAF) for length of stay
from a named list of profile data frames, each containing
resistance-profile probabilities and profile-level relative risk of LOS.

## Usage

``` r
daly_calc_paf_los(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
)
```

## Arguments

- profiles_with_rr:

  Named list returned by `assign_rr_to_profiles()` or
  `filter_profiles_to_rr_classes()`. Each entry is a data frame of
  resistance profiles for one pathogen.

- probability_col:

  Character. Column name for profile probability (must sum to 1 within
  each pathogen). Default `"probability"`.

- rr_profile_col:

  Character. Column name for profile-level LOS relative risk. Default
  `"RR_LOS_profile"`.

- profile_col:

  Character. Column name for the profile identifier. Default
  `"profile"`.

## Value

Named list (one entry per pathogen) with columns `profile`,
`probability`, `rr_profile_col`, `numerator`, `PAF_LOS`, and
`denominator` added to the input profile data frame.

## Details

For each pathogen the formula is: PAF = sum_d R'\_kd \* (RR_kd - 1) /
(1 + sum_d R'\_kd \* (RR_kd - 1))

where R'\_kd is the probability of resistance profile d and RR_kd is the
corresponding relative LOS multiplier.
