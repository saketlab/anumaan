# Estimate Resistance Profile Probabilities via Convex Optimisation

The public-facing solver for Pathway 1. Accepts pre-computed marginal
and pairwise resistance rates (from
[`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md)
/
[`compute_pairwise_from_data()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_from_data.md),
or from
[`validate_aggregate_inputs()`](https://saketlab.github.io/anumaan/reference/validate_aggregate_inputs.md)
for GBD-style pre-computed inputs), enumerates all \\2^n\\ binary
resistance profiles, builds the constraint matrix, and solves the
simplex-constrained weighted least-squares QP:

## Usage

``` r
estimate_profiles_convex(
  marginals,
  pairwise = NULL,
  panel_map = NULL,
  exclude_near_zero = TRUE,
  zero_threshold = 0,
  top_n_classes = NULL,
  lambda = 1e-08,
  sigma_sq = 1,
  solver = c("osqp", "quadprog"),
  n_cores = 1L,
  pathogen_col = "pathogen",
  class_col = "antibiotic_class",
  rate_col = "marginal_resistance",
  n_tested_col = "n_tested"
)
```

## Arguments

- marginals:

  Tibble. Output of
  [`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md)
  or a validated aggregate table. Must have at minimum `pathogen_col`,
  `class_col`, and `rate_col`.

- pairwise:

  Tibble or `NULL`. Output of
  [`compute_pairwise_from_data()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_from_data.md).
  When `NULL`, independence \\P(A \cap B) = P(A)P(B)\\ is assumed for
  all class pairs. Default `NULL`.

- panel_map:

  Named list or `NULL`. If supplied, classes for each pathogen are
  restricted to the panel. When `NULL`, all classes present in
  `marginals` are used. Default `NULL`.

- exclude_near_zero:

  Logical. Drop classes with `marginal_resistance <= zero_threshold`.
  Default `TRUE`.

- zero_threshold:

  Numeric. Classes at or below this value are considered near-zero.
  Default `0`.

- top_n_classes:

  Integer or `NULL`. Cap classes per pathogen to top N by `n_tested`.
  Limits 2^n explosion. Default `NULL`.

- lambda:

  Numeric. Ridge regularisation added to QP Hessian diagonal. Default
  `1e-8`.

- sigma_sq:

  Numeric. Uniform constraint variance. Default `1`.

- solver:

  Character. `"osqp"` (default, preferred) or `"quadprog"` (fallback).

- n_cores:

  Integer. Parallel cores. Default `1L`.

- pathogen_col:

  Character. Default `"pathogen"`.

- class_col:

  Character. Default `"antibiotic_class"`.

- rate_col:

  Character. Default `"marginal_resistance"`.

- n_tested_col:

  Character. Default `"n_tested"`.

## Value

Tibble with one row per (stratum x) pathogen x profile:

- Stratum columns:

  Any geography / year / outcome columns from `marginals`.

- `pathogen`:

  Pathogen name.

- `profile_set_type`:

  `"aggregate_convex"`.

- `profile_class_set`:

  Ordered class names joined by `"|"`.

- `profile_delta`:

  Binary profile label, e.g. `"RSS"`.

- `profile_probability`:

  Estimated probability \\\hat{p}\_\delta\\.

- `estimator`:

  `"convex"`.

- `convergence_flag`:

  Logical; `FALSE` if QP failed and uniform distribution was returned.

- `identifiability_flag`:

  Logical; `TRUE` when the constraint matrix is rank-deficient (system
  is underdetermined).

- `max_abs_residual`:

  Maximum absolute constraint residual.

- `notes`:

  Semicolon-separated notes: capped pairs, fallback pairs.

## Details

\$\$\hat{p} = \arg\min\_{p \in \Delta} \\Mp - v\\^2 + \lambda
\\p\\^2\$\$

Stratification columns (geography, year, outcome) present in `marginals`
are detected automatically and carried through to the output.

## Examples

``` r
if (requireNamespace("quadprog", quietly = TRUE)) {
  marginals <- data.frame(
    pathogen = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 3),
    antibiotic_class = rep(c("Carbapenems", "3GC", "Fluoroquinolones"), 2),
    n_tested = c(420L, 460L, 390L, 280L, 310L, 265L),
    n_resistant = c(126L, 299L, 195L, 45L, 155L, 133L)
  )
  marginals$marginal_resistance <- marginals$n_resistant / marginals$n_tested
  out <- estimate_profiles_convex(marginals,
    panel_map = list(
      "Klebsiella pneumoniae" = c("Carbapenems", "3GC"),
      "Escherichia coli" = c("Carbapenems", "3GC")
    ),
    pathogen_col = "pathogen", class_col = "antibiotic_class",
    rate_col = "marginal_resistance", n_tested_col = "n_tested"
  )
  out[, c("pathogen", "profile_delta", "profile_probability")]
}
#> [estimate_profiles_convex] 2 stratum x pathogen combination(s) | solver: osqp | n_cores: 1
#>   [1/2] 'Klebsiella pneumoniae': n=2 classes -> 4 profiles | max|resid|=1.05e-09 [UNDERDETERMINED]
#>   [2/2] 'Escherichia coli': n=2 classes -> 4 profiles | max|resid|=1.70e-09 [UNDERDETERMINED]
#> 
#> [estimate_profiles_convex] Done. 8 row(s) | 2 pathogen(s) | 1 stratum combination(s).
#> # A tibble: 8 × 3
#>   pathogen              profile_delta profile_probability
#>   <chr>                 <chr>                       <dbl>
#> 1 Klebsiella pneumoniae SS                         0.245 
#> 2 Klebsiella pneumoniae RS                         0.455 
#> 3 Klebsiella pneumoniae SR                         0.105 
#> 4 Klebsiella pneumoniae RR                         0.195 
#> 5 Escherichia coli      SS                         0.420 
#> 6 Escherichia coli      RS                         0.420 
#> 7 Escherichia coli      SR                         0.0804
#> 8 Escherichia coli      RR                         0.0804
```
