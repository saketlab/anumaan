# Formally Check Resistance Profile Probability Constraints

Validates that the profile probability distributions returned by
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
satisfy all mathematical requirements: non-negativity, sum-to-one, and
reconstruction of the input marginal and pairwise resistance rates
within a specified tolerance.

## Usage

``` r
check_profile_constraints(
  profiles_output,
  marginals = NULL,
  pairwise = NULL,
  tolerance = 1e-06,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rate_col = "marginal_resistance",
  class1_col = "antibiotic_class_1",
  class2_col = "antibiotic_class_2",
  pairwise_rate_col = "pairwise_resistance_prevalence"
)
```

## Arguments

- profiles_output:

  Named list returned by
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).
  Each element must contain `$profiles` (data frame with `probability`
  and binary class columns) and optionally `$constraint_residuals`.

- marginals:

  Data frame or `NULL`. Original marginal resistance rates used as QP
  constraints. When provided, reconstructed marginals are compared
  against these values. Must have columns `pathogen_col`, `class_col`,
  and `rate_col`. Default `NULL`.

- pairwise:

  Data frame or `NULL`. Original pairwise co-resistance rates. When
  provided, reconstructed pairwise values are also checked. Must have
  columns `pathogen_col`, `class1_col`, `class2_col`, and
  `pairwise_rate_col`. Default `NULL`.

- tolerance:

  Numeric. Maximum absolute residual considered a pass. Default `1e-6`.

- pathogen_col:

  Character. Pathogen column in `marginals` / `pairwise`. Default
  `"pathogen"`.

- class_col:

  Character. Class column in `marginals`. Default `"antibiotic_class"`.

- rate_col:

  Character. Marginal rate column in `marginals`. Default
  `"marginal_resistance"`.

- class1_col:

  Character. First class column in `pairwise`. Default
  `"antibiotic_class_1"`.

- class2_col:

  Character. Second class column in `pairwise`. Default
  `"antibiotic_class_2"`.

- pairwise_rate_col:

  Character. Pairwise rate column. Default
  `"pairwise_resistance_prevalence"`.

## Value

A tibble with one row per constraint per pathogen:

- `pathogen`:

  Pathogen name.

- `constraint_type`:

  One of `"nonneg"`, `"sum_to_one"`, `"marginal"`, `"pairwise"`.

- `constraint_name`:

  Specific label (e.g. `"marg_3GC"`, `"pair_3GC_Carbapenems"`).

- `target`:

  The constraint target value (`NA` for `nonneg`).

- `reconstructed`:

  The value implied by \\\hat{p}\\.

- `abs_residual`:

  Absolute difference between target and reconstructed.

- `pass`:

  Logical. `TRUE` if `abs_residual < tolerance` or the structural check
  passed.

- `source`:

  Either `"stored_residuals"` (from solve-time values) or `"recomputed"`
  (from direct reconstruction).

Also prints a per-pathogen summary of passed and failed checks.

## Details

The function uses two complementary sources of evidence:

1.  **Constraint residuals** already stored in
    `profiles_output[[pathogen]]$constraint_residuals` – these are \\M
    \hat{p} - v\\ computed at solve time and are available without
    re-supplying the original rates.

2.  **Direct reconstruction** from the binary class-indicator columns in
    the profiles data frame when `marginals` and/or `pairwise` are
    supplied – independently verifies the residuals and catches any
    post-solve modifications to the probability vector.

## Examples

``` r
if (FALSE) { # \dontrun{
marg <- compute_marginal_resistance(amr_clean)
co_res <- compute_pairwise_coresistance(marg)
rp <- compute_resistance_profiles(marg, co_res)

# Using stored constraint residuals only
checks <- check_profile_constraints(rp)

# Full re-verification with original rates
checks <- check_profile_constraints(
  rp,
  marginals = marg$marginal,
  pairwise  = NULL,
  tolerance = 1e-5
)

# Inspect failures
checks[!checks$pass, ]
} # }
```
