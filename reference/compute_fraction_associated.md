# Compute Associated-Burden Fractions per Resistance Profile

Computes the fraction of total expected YLD burden that \*occurs in\*
infections with each resistance profile delta (YLDs associated with
resistance). This is a \*\*burden partition\*\*, not a counterfactual.

## Usage

``` r
compute_fraction_associated(
  profiles_with_rr,
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
)
```

## Arguments

- profiles_with_rr:

  Named list from assign_rr_to_profiles() or
  filter_profiles_to_rr_classes(). Each entry is a profile data frame.

- probability_col:

  Character. Profile probability column. Default `"probability"`.

- rr_profile_col:

  Character. Profile-level RR column. Default `"RR_LOS_profile"`.

## Value

Named list (one entry per pathogen) containing:

- `per_profile`: profile data frame augmented with `numerator_assoc` (=
  p \* rr) and `fraction_assoc` (= p \* rr / E_RR_k).

- `Fraction_k`: overall associated fraction for the pathogen (sum of
  `fraction_assoc` over all resistant profiles).

- `E_RR_k`: expected RR = sum_delta R'\_K_delta \* RR_K_delta.

## Details

For a single drug-class d the formula simplifies to:

Fraction_assoc_Kd = R'\_Kd \* RR_Kd / \[(1 - R'\_Kd) + R'\_Kd \* RR_Kd\]

With resistance profiles the denominator becomes the expected RR across
ALL profiles (including the all-susceptible profile with RR = 1):

E_RR_k = sum_delta R'\_K_delta \* RR_K_delta = 1 + sum_delta R'\_K_delta
\* (RR_K_delta - 1) \[equivalent\]

Per-profile associated fraction: fraction_K_delta = R'\_K_delta \*
RR_K_delta / E_RR_k

Overall associated fraction (all resistant profiles combined):
Fraction_k = sum\_{delta != 0} fraction_K_delta

where delta != 0 denotes profiles with at least one resistant class.

Note: E_RR_k is numerically identical to the \`denominator\` produced by
compute_paf_los() – both equal 1 + sum_d R'\_kd\*(RR_kd - 1).
