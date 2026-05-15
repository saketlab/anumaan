# Compute Mortality Population Attributable Fraction per Resistance Profile

Computes PAF_kd_mortality for each pathogen k and resistance profile
delta using the GBD multi-exposure Levin formula, substituting the
mortality odds ratio (OR_death) from `fit_mortality_rr_logistic()` in
place of the LOS relative risk used by `compute_paf_los()`:

## Usage

``` r
daly_calc_paf_mortality(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile"
)
```

## Arguments

- profiles_with_rr:

  Named list from `assign_rr_to_profiles()` called with
  `rr_col = "OR_death"`. Each element is a profile data frame for one
  pathogen.

- probability_col:

  Character. Profile probability column name. Default `"probability"`.

- rr_profile_col:

  Character. Profile-level OR column as produced by
  `assign_rr_to_profiles()`. Default `"RR_LOS_profile"`.

- profile_col:

  Character. Profile label column. Default `"profile"`.

## Value

Named list (one entry per pathogen) containing:

- `per_profile`: profile data frame augmented with `numerator_mort` (=
  \\R'\_{K\delta}(\text{OR}\_{K\delta}-1)\\), `PAF_mortality` (=
  numerator / denominator), and `denominator_mort` (= \\1 +
  \sum\_\delta\\ numerator).

- `PAF_k_mort`: overall mortality PAF for pathogen k.

- `denominator_mort`: shared denominator \\E\[\text{OR}\_k\]\\.

## Details

\$\$\text{PAF}\_{kd,\text{mort}} =
\frac{R'\_{K\delta}\\(\text{OR}\_{K\delta} - 1)} {1 + \sum\_\delta
R'\_{K\delta}\\(\text{OR}\_{K\delta} - 1)}\$\$

The all-susceptible profile carries OR = 1 and contributes 0. The
denominator equals \\E\[\text{OR}\_k\] = \sum\_\delta R'\_{K\delta}
\cdot \text{OR}\_{K\delta}\\, numerically identical to the denominator
produced by `compute_paf_los()`.

Overall mortality PAF for pathogen k:

\$\$\text{PAF}\_{k,\text{mort}} = \sum\_\delta
\text{PAF}\_{kd,\text{mort}} = \frac{\sum\_\delta
R'\_{K\delta}(\text{OR}\_{K\delta}-1)} {1 + \sum\_\delta
R'\_{K\delta}(\text{OR}\_{K\delta}-1)}\$\$

**Usage pipeline:**

      # 1. Fit mortality OR per class
      mort_or <- fit_mortality_rr_logistic(data, ...)

      # 2. Assign OR to profiles via max rule (rr_col = "OR_death")
      profiles_with_or <- assign_rr_to_profiles(
          profiles_output,
          rr_table = mort_or,
          rr_col   = "OR_death"
      )

      # 3. Compute per-profile and overall mortality PAF
      paf_mort <- daly_calc_paf_mortality(profiles_with_or)
