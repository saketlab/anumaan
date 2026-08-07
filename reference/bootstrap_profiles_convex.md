# Bootstrap Uncertainty Intervals for Resistance Profile Probabilities

Quantifies uncertainty in resistance-profile probability estimates by
repeatedly resampling isolate counts from a binomial distribution
implied by the observed `n_tested` and `n_resistant` values, refitting
the convex optimisation QP for each replicate, and returning percentile
confidence intervals across replicates.

## Usage

``` r
bootstrap_profiles_convex(
  marginals,
  coresistance_output = NULL,
  B = 500L,
  seed = 123L,
  alpha = 0.05,
  n_cores = 1L,
  exclude_near_zero = TRUE,
  top_n_classes = NULL,
  sigma_sq = 1,
  ridge = 1e-08,
  pathogen_col = "organism_name",
  class_col = "antibiotic_class",
  n_tested_col = "n_tested",
  n_resistant_col = "n_resistant",
  org_group_col = "org_group"
)
```

## Arguments

- marginals:

  Data frame. Marginal resistance rates with at minimum columns
  `pathogen_col`, `class_col`, `n_tested_col`, `n_resistant_col`. This
  is the `$marginal` element from
  [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md),
  or a flat tibble with the same structure from
  [`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md).

- coresistance_output:

  Named list or `NULL`. Output of
  [`compute_pairwise_coresistance()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_coresistance.md)
  (one entry per pathogen with `$T_matrix`, `$R_matrix`, `$prevalence`).
  When supplied, pairwise co-resistance is resampled per replicate. When
  `NULL`, independence fallback \\P(A \cap B) = P(A) P(B)\\ is used in
  every replicate. Default `NULL`.

- B:

  Integer. Number of bootstrap replicates. Default `500`.

- seed:

  Integer. Random seed for reproducibility. Default `123`.

- alpha:

  Numeric. Two-sided interval width. Default `0.05` (i.e., 95%
  intervals).

- n_cores:

  Integer. Parallel cores via
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html).
  Default `1L`.

- exclude_near_zero:

  Logical. Passed to
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md).
  Default `TRUE`.

- top_n_classes:

  Integer or `NULL`. Cap on number of classes per pathogen. Default
  `NULL`.

- sigma_sq:

  Numeric. QP constraint variance. Default `1`.

- ridge:

  Numeric. QP Hessian ridge term. Default `1e-8`.

- pathogen_col:

  Character. Pathogen column in `marginals`. Default `"organism_name"`.

- class_col:

  Character. Class column in `marginals`. Default `"antibiotic_class"`.

- n_tested_col:

  Character. Tested-count column. Default `"n_tested"`.

- n_resistant_col:

  Character. Resistant-count column. Default `"n_resistant"`.

- org_group_col:

  Character. Organism-group column (passed through to QP engine).
  Default `"org_group"`.

## Value

A named list, one entry per pathogen. Each entry is a tibble with
columns:

- `profile`:

  Profile label (e.g. `"RSS"`).

- `probability_mean`:

  Mean profile probability across B replicates.

- `probability_median`:

  Median across replicates.

- `lower`:

  Lower percentile (`alpha / 2`).

- `upper`:

  Upper percentile (`1 - alpha / 2`).

- `n_replicates_converged`:

  Number of replicates for which the QP converged (non-uniform
  solution).

- `convergence_rate`:

  Proportion of replicates that converged.

The point-estimate `profiles` tibble (from a single run on the original
marginals) is stored as the attribute `"point_estimate"`.

## Details

**Resampling mechanism:** For each bootstrap replicate \\b\\ and each
(pathogen, class) cell: \$\$n\_{\text{resistant}}^{(b)} \sim
\text{Binomial}(n\_{\text{tested}},\\ \hat{r}\_{kd})\$\$ The bootstrap
marginal is then \\\hat{r}\_{kd}^{(b)} = n\_{\text{resistant}}^{(b)} /
n\_{\text{tested}}\\. Pairwise co-resistance rates (if supplied) are
resampled analogously. The QP is re-solved on each bootstrap marginal
using the existing
[`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
engine.

**Performance:** For \\n \leq 12\\ classes the QP solves in
milliseconds; B = 500 replicates will complete in seconds. For \\n \geq
14\\ classes use `n_cores > 1` or reduce B.

## Examples

``` r
if (FALSE) { # \dontrun{
marg <- compute_marginal_resistance(amr_clean)
co_res <- compute_pairwise_coresistance(marg)

boot <- bootstrap_profiles_convex(
  marginals = marg$marginal,
  coresistance_output = co_res,
  B = 500,
  seed = 42
)

# 95% intervals for K. pneumoniae
boot[["Klebsiella pneumoniae"]]

# Point estimate
attr(boot[["Klebsiella pneumoniae"]], "point_estimate")
} # }
```
