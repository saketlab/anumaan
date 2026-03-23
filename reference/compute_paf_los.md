# Compute LOS Population Attributable Fraction per Resistance Profile

Computes PAF_kd_LOS for each pathogen k and profile d using the GBD
multi-exposure Levin formula:

## Usage

``` r
compute_paf_los(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
)
```

## Arguments

- profiles_with_rr:

  Named list from assign_rr_to_profiles().

- probability_col:

  Character. Profile probability column. Default `"probability"`.

- rr_profile_col:

  Character. Profile-level RR column. Default `"RR_LOS_profile"`.

- profile_col:

  Character. Profile label column. Default `"profile"`.

## Value

Named list (one entry per pathogen): per_profile: profiles data frame
augmented with numerator, PAF_LOS, denominator. PAF_k: overall PAF for
pathogen k. denominator: 1 + sum_d R'\_kd(RR_kd-1).

## Details

PAF_kd_LOS = R'\_kd \* (RR_kd - 1) / \[1 + sum_d R'\_kd \* (RR_kd - 1)\]

where R'\_kd is the profile probability from
compute_resistance_profiles() and RR_kd is from assign_rr_to_profiles().
The denominator equals E_d\[RR_kd\] (expected RR over all profiles) and
is shared across profiles for pathogen k. The all-susceptible profile
contributes 0 (RR = 1).

Overall PAF for pathogen k: PAF_k = sum_d PAF_kd_LOS = \[sum_d
R'\_kd(RR_kd-1)\] / \[1 + sum_d R'\_kd(RR_kd-1)\]
