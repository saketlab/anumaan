# Classify why a fit's sampler diagnostics look the way they do

Produces a ranked primary/secondary/not-indicated recommendation set
instead of one generic line, so a reviewer is told the SPECIFIC failure
mode (divergence vs. energy vs. Omega geometry vs. fixed-effect sparsity
vs. treedepth) rather than a boilerplate "increase adapt_delta or
simplify the model" that is only actually appropriate for divergences.

## Usage

``` r
.probit_classify_fitting_issues(
  diag,
  grouped_diag,
  degeneracy_stats,
  residual_structure,
  beta_summary = NULL,
  thresholds = .probit_diagnostic_thresholds()
)
```

## Arguments

- diag:

  `fit$diagnostics` (canonical scalar diagnostic fields).

- grouped_diag:

  `fit$diagnostics_detail$grouped` – canonical per-parameter_group
  Rhat/ESS table (rows: beta, Omega, lp\_\_, z_free, any declared
  random-effect blocks). May be `NULL`.

- degeneracy_stats:

  Result of
  [`.omega_degeneracy_stats()`](https://saketlab.github.io/anumaan/reference/dot-omega_degeneracy_stats.md),
  or `NULL` for identity-residual fits.

- residual_structure:

  `"identity"` or `"correlated"`.

- beta_summary:

  Per-covariate-family beta Rhat/ESS breakdown from
  `plot_probit_beta_family_diagnostics()$summary`, or `NULL`.

- thresholds:

  List from
  [`.probit_diagnostic_thresholds()`](https://saketlab.github.io/anumaan/reference/dot-probit_diagnostic_thresholds.md).

## Value

Named list: `primary`, `secondary`, `not_indicated` (character vectors,
possibly empty) and `bucket` (single string identifying which top-level
failure mode was selected as primary, for testing).
