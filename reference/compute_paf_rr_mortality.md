# Compute Mortality Population Attributable Fraction per Resistance Profile

Computes PAF_kd_mortality for each pathogen k and resistance profile
delta using the GBD multi-exposure Levin formula, substituting the
mortality odds ratio (OR_death) from
[`fit_mortality_rr_logistic()`](https://saketlab.github.io/anumaan/reference/fit_mortality_rr_logistic.md)
in place of the LOS relative risk used by
[`compute_paf_los()`](https://saketlab.github.io/anumaan/reference/compute_paf_los.md):

Computes PAF_kd_mortality for each pathogen x resistance profile using
the GBD multi-exposure Levin formula with mortality OR substituted for
LOS RR:

## Usage

``` r
compute_paf_rr_mortality(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile",
  facility_col = NULL,
  facility_name = NULL
)

compute_paf_rr_mortality(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  profile_col = "profile",
  facility_col = NULL,
  facility_name = NULL
)
```

## Arguments

- profiles_with_rr:

  Named list from
  [`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
  called with `rr_col = "OR_death"`.

- probability_col:

  Character. Default `"probability"`.

- rr_profile_col:

  Character. Default `"RR_LOS_profile"`.

- profile_col:

  Character. Default `"profile"`.

- facility_col:

  Character or `NULL`. Facility identifier column. Default `NULL`.

- facility_name:

  Character or `NULL`. If provided, stored in the output for provenance
  tracking. Default `NULL`.

## Value

Named list (one entry per pathogen) containing:

- `per_profile`: profile data frame augmented with `numerator_mort` (=
  \\R'\_{K\delta}(\text{OR}\_{K\delta}-1)\\), `PAF_mortality` (=
  numerator / denominator), and `denominator_mort` (= \\1 +
  \sum\_\delta\\ numerator).

- `PAF_k_mort`: overall mortality PAF for pathogen k.

- `denominator_mort`: shared denominator \\E\[\text{OR}\_k\]\\.

Named list (one per pathogen) with `per_profile`, `PAF_k_mort`, and
`denominator_mort`.

## Details

\$\$\text{PAF}\_{kd,\text{mort}} =
\frac{R'\_{K\delta}\\(\text{OR}\_{K\delta} - 1)} {1 + \sum\_\delta
R'\_{K\delta}\\(\text{OR}\_{K\delta} - 1)}\$\$

The all-susceptible profile carries OR = 1 and contributes 0. The
denominator equals \\E\[\text{OR}\_k\] = \sum\_\delta R'\_{K\delta}
\cdot \text{OR}\_{K\delta}\\, numerically identical to the denominator
produced by
[`compute_paf_los()`](https://saketlab.github.io/anumaan/reference/compute_paf_los.md).

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
      paf_mort <- compute_paf_rr_mortality(profiles_with_or)

\$\$\text{PAF}\_{kd} = \frac{R'\_{K\delta}(OR\_{K\delta} - 1)} {1 +
\sum\_\delta R'\_{K\delta}(OR\_{K\delta} - 1)}\$\$

**Pipeline:**

      or_table          <- fit_mortality_rr_logistic(data, ...)
      p0_res            <- compute_p0(data, ...)
      rr_table          <- convert_or_to_rr(or_table, p0 = p0_res$p0)
      profiles_with_or  <- assign_rr_to_profiles(profiles_output,
                               rr_table = or_table, rr_col = "OR_death")
      paf_mort          <- compute_paf_rr_mortality(profiles_with_or)
