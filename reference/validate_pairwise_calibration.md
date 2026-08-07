# Observed-versus-Model Pairwise Co-resistance Validation

For every adequately co-tested hospital x pathogen x class-pair
combination, compares the observed pairwise resistance rate
(\\n\_{RR}/n\_{\text{cotested}}\\) against the model's implied joint
probability, computed differently by residual structure:

- Identity:

  \\P(Y\_{d_1}=1, Y\_{d_2}=1 \mid \theta) =
  \Phi(\mu\_{ed_1})\\\Phi(\mu\_{ed_2})\\ – the two classes are combined
  only through the shared fixed and random effects in \\\mu\\, not
  through any residual correlation (the identity model has none). **This
  checks whether the identity-residual model reproduces observed
  pairwise resistance through shared covariates and random effects; it
  is not evidence that a residual-independence assumption is correct.**
  Persistent pairwise miscalibration despite good marginal calibration
  suggests unmodelled class-to-class dependence that a
  correlated-residual fit might capture.

- Correlated:

  \\P(Y\_{d_1}=1, Y\_{d_2}=1 \mid \theta) = \Phi_2(\mu\_{ed_1},
  \mu\_{ed_2}; \rho\_{d_1 d_2})\\ – the standard bivariate-probit
  identity \\P(Z_1\>0,Z_2\>0) = \Phi_2(\mu_1,\mu_2;\rho)\\ for
  \\(Z_1,Z_2)\sim N((\mu_1,\mu_2), \[\[1,\rho\],\[\rho,1\]\])\\, with
  \\\rho\_{d_1 d_2} = \Omega\_{d_1 d_2}\\ from the fitted residual
  correlation matrix (computed via pbivnorm). This actually uses the
  fitted correlation structure, unlike the identity-model product
  formula – checking whether \\\Omega\\ reproduces observed
  co-resistance.

## Usage

``` r
validate_pairwise_calibration(
  fitted_model,
  n_posterior_draws_for_validation = 2000L,
  seed = 123L,
  ci_level = 0.95,
  min_cotested = NULL
)
```

## Arguments

- fitted_model:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- n_posterior_draws_for_validation:

  Integer. Posterior draws used. Default `2000L`.

- seed:

  Integer. Random seed (draw subsampling only). Default `123L`.

- ci_level:

  Numeric. Credible interval coverage. Default `0.95`.

- min_cotested:

  Integer or `NULL`. Eligibility threshold. When `NULL` (default),
  eligibility is taken from `fitted_model$eligibility_report$pairwise`.

## Value

Tibble, one row per eligible hospital x pathogen x class-pair:
`n_cotested`, `observed_RR/RS/SR/SS`, `observed_pairwise_resistance`,
`model_pairwise_mean/lower/upper`, `absolute_error`,
`interval_contains_observed`.
