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
if (requireNamespace("quadprog", quietly = TRUE)) {
  set.seed(1)
  ast_class <- data.frame(
    isolate_id       = rep(sprintf("ISO%02d", 1:80), each = 3),
    organism_name    = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 120),
    org_group        = "Enterobacterales",
    antibiotic_class = rep(c("Carbapenems", "3GC", "Fluoroquinolones"), times = 80),
    antibiotic_value = sample(c("S", "R"), 240, replace = TRUE)
  )
  marg <- compute_marginal_resistance(ast_class)
  marg$marginal
  co_res <- compute_pairwise_coresistance(marg)
  co_res[["Klebsiella pneumoniae"]]$prevalence
  rp <- compute_resistance_profiles(marg, co_res, exclude_near_zero = FALSE)
  rp[["Klebsiella pneumoniae"]]$profiles
  checks <- check_profile_constraints(rp,
    marginals = marg$marginal,
    pathogen_col = "organism_name", class_col = "antibiotic_class",
    rate_col = "marginal_resistance"
  )
  checks[!checks$pass, ]
}
#> Collapsed to class level: 80 isolates | 2 pathogens | 3 classes.
#> No classes with marginal_resistance <= 0.
#> min_n_tested = 30: all 6 pathogen-class combination(s) passed the threshold.
#> class_long after min_n_tested filter: 240 isolate-class row(s) retained.
#> 'Escherichia coli': 3 classes, 40 isolates -- co-resistance matrix built.
#> 'Klebsiella pneumoniae': 3 classes, 40 isolates -- co-resistance matrix built.
#> 
#> Processing 2 pathogen(s) with 1 core(s)...
#> 
#> [1/2] Processing 'Escherichia coli'...
#> 'Escherichia coli': n=3 classes -> 8 profiles. Max |residual| = 0.00000.
#> 
#> [2/2] Processing 'Klebsiella pneumoniae'...
#> 'Klebsiella pneumoniae': n=3 classes -> 8 profiles. Max |residual| = 0.00000.
#> 
#> [check_profile_constraints] Results (tolerance = 1e-06):
#>   [OK] Escherichia coli                          11/11 passed | max |residual| = 8.93e-10
#>   [OK] Klebsiella pneumoniae                     11/11 passed | max |residual| = 5.14e-10
#> # A tibble: 0 × 8
#> # ℹ 8 variables: pathogen <chr>, constraint_type <chr>, constraint_name <chr>,
#> #   target <dbl>, reconstructed <dbl>, abs_residual <dbl>, pass <lgl>,
#> #   source <chr>
```
