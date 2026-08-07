# Summarise Validation Checks into a Profile-Validation Status

Combines the outputs of
[`validate_marginal_calibration()`](https://saketlab.github.io/anumaan/reference/validate_marginal_calibration.md),
[`validate_pairwise_calibration()`](https://saketlab.github.io/anumaan/reference/validate_pairwise_calibration.md),
[`validate_complete_profile_calibration()`](https://saketlab.github.io/anumaan/reference/validate_complete_profile_calibration.md),
and
[`mask_and_validate_ast()`](https://saketlab.github.io/anumaan/reference/mask_and_validate_ast.md)
into a single `profile_validation_status` field. **This is deliberately
separate from `diagnostic_status`** (see
[`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)):
`diagnostic_status` describes HMC sampler/model-fitting health (R-hat,
ESS, divergences, E-BFMI, tree-depth); `profile_validation_status`
describes predictive agreement with the observed data and data support
for that agreement. A model can converge cleanly
(`diagnostic_status = "pass"`) while its profile predictions are poorly
calibrated, or vice versa, and the two should never be merged into one
flag.

## Usage

``` r
compute_profile_validation_status(
  marginal_tbl = NULL,
  pairwise_tbl = NULL,
  complete_profile_tbl = NULL,
  masked_summary = NULL,
  thresholds = list()
)
```

## Arguments

- marginal_tbl:

  Tibble from
  [`validate_marginal_calibration()`](https://saketlab.github.io/anumaan/reference/validate_marginal_calibration.md),
  or `NULL` if not run.

- pairwise_tbl:

  Tibble from
  [`validate_pairwise_calibration()`](https://saketlab.github.io/anumaan/reference/validate_pairwise_calibration.md),
  or `NULL` if not run.

- complete_profile_tbl:

  Tibble from
  [`validate_complete_profile_calibration()`](https://saketlab.github.io/anumaan/reference/validate_complete_profile_calibration.md),
  or `NULL` if not run.

- masked_summary:

  Tibble (`$summary` element) from
  [`mask_and_validate_ast()`](https://saketlab.github.io/anumaan/reference/mask_and_validate_ast.md),
  or `NULL` if not run.

- thresholds:

  Named list overriding any of: `max_mean_abs_error_marginal` (default
  `0.10`), `min_marginal_coverage` (default `0.80`),
  `max_mean_abs_error_pairwise` (default `0.10`),
  `min_pairwise_coverage` (default `0.80`),
  `min_complete_profile_panels_evaluated` (default `1L`),
  `max_masked_brier` (default `0.25`), `min_masked_auroc` (default
  `0.6`).

## Value

Named list: `status` (one of `"pass"`, `"in_sample_validation_only"`,
`"warning_marginal_miscalibration"`,
`"warning_pairwise_miscalibration"`,
`"warning_sparse_complete_profiles"`, `"fail_masked_ast_validation"`, or
`"not_evaluated"` when no inputs were supplied), `reasons` (character
vector of every triggered condition – including
`"masked_validation_in_sample_only"` as an informational reason even
when a more severe status wins – the most severe status is returned but
all are listed), and `thresholds_used` (the resolved threshold list, for
audit trail).

## Details

All thresholds are explicit parameters with documented defaults – none
are silently hard-coded elsewhere. Pass `thresholds` to override any
subset; unspecified entries keep their default.

**In-sample masked validation can never produce a bare `"pass"` on its
own.** When `masked_summary$validation_mode == "in_sample_no_refit"`
(see
[`mask_and_validate_ast()`](https://saketlab.github.io/anumaan/reference/mask_and_validate_ast.md)),
the masked cells were already seen during fitting, so predicting them
well is not genuine validation. If every threshold is otherwise cleared,
the status becomes `"in_sample_validation_only"` rather than `"pass"` –
a real `"pass"` requires either `validation_mode == "refit_masked"` (a
genuine holdout) or no masked-AST check at all (the other
observed-versus-model checks being used purely descriptively).
