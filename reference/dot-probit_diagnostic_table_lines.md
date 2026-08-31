# Compact diagnostic-summary table for the interpretation page

Compact diagnostic-summary table for the interpretation page

## Usage

``` r
.probit_diagnostic_table_lines(
  diag,
  degeneracy_stats,
  residual_structure,
  thresholds = .probit_diagnostic_thresholds()
)
```

## Arguments

- diag:

  `fit$diagnostics` (canonical scalar diagnostic fields).

- degeneracy_stats:

  Result of
  [`.omega_degeneracy_stats()`](https://saketlab.github.io/anumaan/reference/dot-omega_degeneracy_stats.md),
  or `NULL` for identity-residual fits.

- residual_structure:

  `"identity"` or `"correlated"`.

- thresholds:

  List from
  [`.probit_diagnostic_thresholds()`](https://saketlab.github.io/anumaan/reference/dot-probit_diagnostic_thresholds.md).

## Value

Character vector of pre-formatted monospace table lines (header,
separator, one row per diagnostic).
