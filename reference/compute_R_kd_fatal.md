# Compute Fatal Prevalence of Resistance (R_kd)

Computes the fatal prevalence of resistance R_kd for each pathogen k and
resistance profile delta, using per-profile mortality odds ratios
(OR_death) from
[`fit_mortality_rr_logistic()`](https://saketlab.github.io/anumaan/reference/fit_mortality_rr_logistic.md)
in place of LOS relative risks. The formula is identical to
[`compute_fraction_associated()`](https://saketlab.github.io/anumaan/reference/compute_fraction_associated.md)
– only the RR source changes.

## Usage

``` r
compute_R_kd_fatal(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
)
```

## Arguments

- profiles_with_rr:

  Named list from
  [`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
  called with `rr_col = "OR_death"`.

- probability_col:

  Character. Profile probability column. Default `"probability"`.

- rr_profile_col:

  Character. Profile-level OR column as produced by
  [`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md).
  Default `"RR_LOS_profile"`.

## Value

Named list (one entry per pathogen) containing:

- `per_profile`: profile data frame augmented with `numerator_R_kd` (=
  \\R'\_{K\delta} \cdot \text{OR}\_{K\delta}\\) and `R_kd` (= numerator
  / E\[OR_k\]) per profile.

- `R_k`: overall fatal prevalence of resistance for pathogen k (sum of
  `R_kd` over resistant profiles only).

- `E_OR_k`: expected OR = \\\sum\_\delta R'\_{K\delta} \cdot
  \text{OR}\_{K\delta}\\.

## Details

\$\$E\[\text{OR}\_k\] = \sum\_\delta R'\_{K\delta} \cdot
\text{OR}\_{K\delta}\$\$

\$\$R\_{kd} = \frac{R'\_{K\delta} \cdot
\text{OR}\_{K\delta}}{E\[\text{OR}\_k\]}\$\$

\$\$R\_{k} = \sum\_{\delta \neq \delta_0} R\_{kd}\$\$

where \\R'\_{K\delta}\\ is the profile probability from
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
and \\\text{OR}\_{K\delta}\\ is assigned by
[`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
with `rr_col = "OR_death"` (the all-susceptible profile carries OR = 1).

**Usage pipeline:**

      # 1. Fit mortality OR per class
      mort_or <- fit_mortality_rr_logistic(data, ...)

      # 2. Assign OR to profiles via max rule
      profiles_with_or <- assign_rr_to_profiles(
          profiles_output,
          rr_table = mort_or,
          rr_col   = "OR_death"
      )

      # 3. Compute R_kd (fatal prevalence of resistance)
      rkd <- compute_R_kd_fatal(profiles_with_or)
