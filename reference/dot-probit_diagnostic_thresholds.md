# Central thresholds for diagnostic-report interpretation

Single place these thresholds are defined, so the interpretation logic
and the report table stay in sync instead of drifting independently.

## Usage

``` r
.probit_diagnostic_thresholds(
  ebfmi_min = 0.3,
  rhat_max = 1.01,
  rhat_severe = 1.05,
  treedepth_max = 0L,
  divergent_max = 0L,
  omega_small_eigen_threshold = 0.05,
  omega_small_eigen_fraction_warning = 0.5,
  omega_condition_number_warning = 50
)
```

## Arguments

- ebfmi_min:

  Warn below this minimum chain E-BFMI. Default `0.3` (Stan's own
  guidance threshold).

- rhat_max:

  Warn above this Rhat. Default `1.01`.

- rhat_severe:

  Structural Rhat above this is treated as severe rather than
  borderline. Default `1.05`.

- treedepth_max:

  Treedepth-saturation count above this triggers the treedepth failure
  mode. Default `0L` (any saturation).

- divergent_max:

  Divergent-transition count above this triggers the divergence failure
  mode. Default `0L` (any divergence).

- omega_small_eigen_threshold:

  Per-draw smallest-eigenvalue-of-Omega value below which a draw is
  counted as near-singular. Default `0.05`.

- omega_small_eigen_fraction_warning:

  Fraction of posterior draws flagged near-singular (via
  `omega_small_eigen_threshold`) above which Omega geometry is reported
  as a primary concern. Default `0.50`.

- omega_condition_number_warning:

  Median Omega condition number above which it is called out explicitly.
  Default `50`.

## Value

Named list of thresholds.
