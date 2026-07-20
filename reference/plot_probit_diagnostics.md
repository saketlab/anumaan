# Generate diagnostic plots for a fitted Bayesian multivariate probit model

Generate diagnostic plots for a fitted Bayesian multivariate probit
model

## Usage

``` r
plot_probit_diagnostics(
  fit_obj,
  output_dir,
  experiment_id,
  pathogen,
  max_params = 30L
)
```

## Arguments

- fit_obj:

  Named list returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- output_dir:

  Directory where the PDF will be saved.

- experiment_id:

  Character. Experiment identifier (used in filename and titles).

- pathogen:

  Character. Pathogen name (used in plot titles).

- max_params:

  Integer. Maximum number of parameters shown in trace/rank plots.
  Defaults to 30. Includes `lp__`, all `tau_*`, and controlled subsets
  of `beta`, random-effect terms, random-effect correlations, and
  residual-correlation terms when available.

## Value

Invisibly returns the path to the saved PDF, or `NULL` if skipped.
