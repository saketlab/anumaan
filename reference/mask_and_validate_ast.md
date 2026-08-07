# Masked-AST Holdout Validation

Masks a reproducible, stratified subset of genuinely **observed** AST
cells, predicts those cells' resistance probability from a fitted model,
and scores the predictions against their true values.

## Usage

``` r
mask_and_validate_ast(
  fitted_model,
  fixed_effects = NULL,
  event_id_col = "event_id",
  outcome_col = NULL,
  fraction_to_mask = 0.05,
  max_mask_per_cell = NULL,
  seed = 123L,
  refit = FALSE,
  n_posterior_draws_for_validation = 2000L,
  min_tested_after_mask = 30L,
  min_resistant_after_mask = 5L,
  min_susceptible_after_mask = 5L,
  panel_eligibility = list(),
  prior_config = NULL,
  sampler_config = NULL,
  n_gibbs_burnin = 10L,
  n_gibbs_kept = 20L
)
```

## Arguments

- fitted_model:

  List returned by
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- fixed_effects:

  Character vector. Required when `refit = TRUE` (not recoverable from
  `fitted_model` due to dummy-coding of categorical covariates); must
  match the original fitting call.

- event_id_col:

  Character. Column in `fitted_model$event_metadata` holding the
  original event id; only used when `refit = TRUE`. Default
  `"event_id"`.

- outcome_col:

  Character or `NULL`. Only used when `refit = TRUE`.

- fraction_to_mask:

  Numeric in (0, 1\]. Target fraction of observed cells to mask per
  hospital x pathogen x class cell, before the eligibility floor and
  `max_mask_per_cell` caps are applied. Default `0.05`.

- max_mask_per_cell:

  Integer or `NULL`. Optional additional cap on cells masked per
  hospital x pathogen x class cell.

- seed:

  Integer. Random seed for mask selection (and posterior draw
  subsampling). Default `123L`.

- refit:

  Logical. See modes above. Default `FALSE`.

- n_posterior_draws_for_validation:

  Integer. Default `2000L`.

- min_tested_after_mask, min_resistant_after_mask,
  min_susceptible_after_mask:

  Integer. Panel-eligibility floors enforced during masking. Default
  `30L`, `5L`, `5L`.

- panel_eligibility:

  Named list. Forwarded to
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
  when `refit = TRUE`.

- prior_config, sampler_config:

  Named list or `NULL`. Forwarded when `refit = TRUE`; default to the
  original fit's `prior_config_used`/`sampler_config_used`.

- n_gibbs_burnin, n_gibbs_kept:

  Integer. Gibbs burn-in/kept-iteration counts forwarded to
  [`.gibbs_conditional_profile_probs()`](https://saketlab.github.io/anumaan/reference/dot-gibbs_conditional_profile_probs.md)
  for correlated-residual fits' masked-cell prediction. Ignored for
  identity-residual fits. Default `10L`/`20L`.

## Value

Named list: `predictions` (one row per masked cell: hospital, pathogen,
event id, class, true value, predicted probability, log score, Brier
score, predicted class at 0.5, calibration group, validation mode),
`summary` (one-row aggregate: counts, mean log/Brier score, accuracy,
AUROC \\only when both classes are present in the masked set\\, masking
configuration), and `refit_model` (the refit `fitted_model`, or `NULL`
when `refit = FALSE`).

## Details

**Two validation modes:**

- `refit = FALSE` (default):

  Predicts masked cells using the *already-fitted* `fitted_model`.
  Because \\\mu\_{ed}\\ never depends on class \\d\\'s own observed
  value (only on event-level covariates and random-effect memberships),
  this is honestly an **in-sample conditional-prediction diagnostic**
  (`validation_mode = "in_sample_no_refit"`), not a genuine held-out
  test – the masked cells still contributed to estimating
  `beta`/`hospital_effect` etc. It is useful for spotting miscalibrated
  subgroups quickly, but is optimistic relative to a real holdout.

- `refit = TRUE` (preferred, more expensive):

  Sets the selected cells to `NA` in a copy of
  `fitted_model$event_metadata`, re-runs
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
  on the masked data, and predicts the held-out cells from that refit
  (`validation_mode = "refit_masked"`). This is a genuine holdout and is
  the recommended final procedure once computationally affordable.

**Masking procedure:** reproducible via `seed`; stratified by hospital x
pathogen x class x observed R/S state; for every hospital x pathogen x
class cell, at most enough cells are masked to keep
`n_tested/n_resistant/n_susceptible` at or above
`min_tested_after_mask`/`min_resistant_after_mask`/
`min_susceptible_after_mask` (panel eligibility is not destroyed); only
observed (non-`NA`) cells are ever candidates; and no event is ever left
with zero observed classes across *all* `class_cols` (which would
otherwise silently drop that event from a refit).

**Prediction**: for `fitted_model$residual_structure == "identity"`,
\\\Phi(\mu\_{ed})\\ is the exact model-implied probability for a masked
cell regardless of what else is known about that event (classes are
conditionally independent given \\\mu\\). For `"correlated"`, using
\\\Phi(\mu\_{ed})\\ alone would ignore every OTHER still-observed class
for that event – wrong whenever classes are correlated. Instead each
masked event's full `class_cols` panel is rebuilt with the masked
cell(s) set to `NA` and every other class at its true observed value,
and
[`.gibbs_conditional_profile_probs()`](https://saketlab.github.io/anumaan/reference/dot-gibbs_conditional_profile_probs.md)
is used to get \\P(\text{masked class}=R \mid \text{that event's other
observed classes}, \theta)\\ – the same conditional-imputation machinery
used for DALY profile generation, not a marginal prediction.
