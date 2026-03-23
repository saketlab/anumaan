# Compute Per-Profile Fatal Prevalence of Resistance R_Kdelta\* (Eq. 13)

For each \*\*resistant\*\* profile delta of pathogen K, computes the
fatal prevalence of resistance:

## Usage

``` r
compute_fatal_resistance_prevalence(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile",
  dominant_class_col = "dominant_class",
  facility_col = NULL,
  facility_name = NULL
)
```

## Arguments

- profiles_with_rr:

  Named list from `assign_rr_to_profiles(rr_col = "RR_death")`, one
  element per pathogen. Each element must contain `probability_col`,
  `rr_profile_col`, and `dominant_class_col`.

- probability_col:

  Character. Profile prevalence column R'\_Kdelta. Default
  `"probability"`.

- rr_profile_col:

  Character. Dominant-class converted RR. Default `"RR_LOS_profile"`.

- dominant_class_col:

  Character. Column identifying the dominant class (or
  `"all_susceptible"`). Default `"dominant_class"`.

- facility_col:

  Character or `NULL`. Facility identifier column. Default `NULL`.

- facility_name:

  Character or `NULL`. If provided, stored in the output for provenance
  tracking. Default `NULL`.

## Value

Named list (one per pathogen). Each element:

- `per_profile`:

  Data frame of resistant profiles augmented with `R_star` (R\*\_Kdelta
  per profile) and `numerator_delta`.

- `R_K_star`:

  Scalar: \\\sum\_\delta R^\*\_{K\delta}\\ = total fatal resistance
  prevalence for pathogen K.

- `sum_r_prime`:

  Scalar: \\\sum\_\delta R'\_{K\delta}\\ (resistant profiles only).

- `susceptible_fraction`:

  Scalar: \\1 - \sum\_\delta R'\_{K\delta}\\.

- `denominator`:

  Scalar: shared denominator for all profiles.

- `n_resistant_profiles`:

  Integer: number of resistant profiles.

## Details

\$\$ R^\*\_{K\delta} = \frac{R'\_{K\delta} \cdot RR\_{Kd^\*}} {\left(1 -
\textstyle\sum\_\delta R'\_{K\delta}\right) + \textstyle\sum\_\delta
R'\_{K\delta} \cdot RR\_{Kd^\*}} \$\$

where the sums \\\sum\_\delta\\ run over **resistant profiles only**
(profiles with at least one resistant antibiotic class, i.e.
`dominant_class != "all_susceptible"`). The term \\(1 - \sum\_\delta
R'\_{K\delta})\\ is the susceptible fraction and must be \> 0 for the
formula to give R\*\_Kdelta \< 1.

**Important:** pass the output of
[`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
*directly* – do *not* call
[`filter_profiles_to_rr_classes()`](https://saketlab.github.io/anumaan/reference/filter_profiles_to_rr_classes.md)
first, because renormalisation would set \\\sum R'\_{K\delta} = 1\\ and
collapse the denominator, producing R\*\_Kdelta = 1 for every pathogen.

## References

Bhaswati Ganguli. DALY Methodology for AMR (YLD notes). March 2026. Eq.
13.
