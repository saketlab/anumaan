# Summarize a posterior D x D correlation matrix from stored fit draws

Extracts posterior mean/median/CI for each off-diagonal cell of a named
correlation-matrix generated quantity (e.g. "Omega", "R_hospital",
"R_patient", "R_admission"). These are three distinct quantities and
must not be combined into one table: Omega is the event-level residual
cross-class correlation (only estimated when residual_structure ==
"correlated"); R_hospital/R_patient/R_admission are the correlation,
across antibiotic classes, of that random effect's own tendencies
(estimated regardless of residual_structure, since the
hospital/patient/admission random effects always use a
diag_pre_multiply(tau, L_corr) parameterisation). Returns NULL if the
requested matrix_var is not present in fit\$draws (e.g. Omega on an
identity-residual fit, or R_admission on a 1- or 2-RE fit).

## Usage

``` r
summarize_fit_correlation_matrix(
  fit,
  matrix_var,
  class_cols,
  ci_level = 0.95,
  block_index = NULL
)
```

## Arguments

- fit:

  A fitted_model object (or its lightweight saved form) with a `$draws`
  posterior::draws_array retaining the requested matrix_var.

- matrix_var:

  Character. One of "Omega", "R_hospital", "R_patient", "R_admission".

- class_cols:

  Character vector of antibiotic class column names, in the same order
  used to fit the model (dimension order of matrix_var).

- ci_level:

  Numeric credible-interval level. Default 0.95.

- block_index:

  Integer or NULL. When the fit used the generic random-effect
  architecture (Stage 1), per-block correlation matrices are emitted as
  `R_block[r,i,j]` (block index r as the FIRST subscript), not
  `R_hospital[i,j]`/`R_patient[i,j]`/etc. Pass `matrix_var = "R_block"`
  and the 1-based block index here to look up a specific block's
  correlation matrix; leave NULL for 2-index variables (`Omega`, or
  legacy `R_hospital` etc. from a fit made with the old hardcoded
  1re/2re/3re Stan models).

## Value

A tibble with one row per off-diagonal class pair, or NULL if matrix_var
was not estimated for this fit.
