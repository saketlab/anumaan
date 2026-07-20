# Fit Bayesian Hierarchical Multivariate Probit Model for Resistance Profiles

Accepts wide-format event-level AST data (one row per organism-event,
one column per antibiotic class with 0 / 1 / `NA` values) and fits a
hierarchical Bayesian multivariate probit model via cmdstanr.

## Usage

``` r
fit_bayesian_multivariate_probit(
  event_class_data,
  class_cols,
  fixed_effects,
  random_effects,
  pathogen = NULL,
  pathogen_col = "pathogen",
  event_id_col = "event_id",
  eligible_pairs = NULL,
  outcome_col = NULL,
  reserve_drug_cols = NULL,
  panel_eligibility = list(),
  residual_structure = c("identity", "correlated"),
  estimand = "observed_stewardship_event_mix",
  prior_config = list(),
  sampler_config = list(),
  show_messages = TRUE,
  ...
)
```

## Arguments

- event_class_data:

  Data frame. One row per organism-event. Antibiotic class columns hold
  0 (susceptible), 1 (resistant), or `NA` (not tested). Metadata columns
  hold covariates and grouping variables.

- class_cols:

  Character vector. Names of the antibiotic class columns. Required.

- fixed_effects:

  Character vector. Event-level covariate column names. Required.

- random_effects:

  Character vector of length 1, 2, or 3. Grouping column names. First
  element: hospital (upper); second (optional): patient; third
  (optional): admission. Required.

- pathogen:

  Character or `NULL`. When supplied, filters `event_class_data` to rows
  where `pathogen_col` equals this value before fitting. Recommended:
  fit one pathogen at a time.

- pathogen_col:

  Character. Column identifying the pathogen. Default `"pathogen"`.

- event_id_col:

  Character. Column in `event_class_data` holding unique event
  identifiers. Default `"event_id"`.

- eligible_pairs:

  Tibble or `NULL`. Hospital x pathogen pairs to include. `NULL` uses
  all pairs present in the data.

- outcome_col:

  Character or `NULL`. Patient outcome column. Only used downstream to
  split R_ALL and R_NF cohorts – does not enter the probit likelihood.
  Default `NULL`.

- reserve_drug_cols:

  Character vector or `NULL`. Class columns to exclude from the main
  model.

- panel_eligibility:

  Named list. Eligibility thresholds per hospital x class cell:
  `min_tested` (default 30), `min_resistant` (default 5),
  `min_susceptible` (default 5). Cells not meeting thresholds are
  reported but fitting still proceeds.

- residual_structure:

  Character. Controls the residual covariance structure. `"identity"`
  (default): classes are conditionally independent given fixed and
  random effects – residual covariance = I_D. `"correlated"`: estimates
  a full residual correlation matrix \\\Omega\\ via LKJCholesky prior.
  Use `"correlated"` only when panel co-testing overlap is adequate
  (check `$eligibility_report$pairwise`); otherwise \\\Omega\\ is driven
  mainly by the LKJ prior and adds identifiability risk.

- estimand:

  Character. Identifies the target quantity. Only
  `"observed_stewardship_event_mix"` is currently supported.

- prior_config:

  Named list. Any subset of `beta_sd` (default 1.5), `tau_sd` (default
  1.0), `lkj_eta` (default 2.0).

- sampler_config:

  Named list. Sampler settings forwarded to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html):
  `chains` (4), `iter_warmup` (1000), `iter_sampling` (1000),
  `adapt_delta` (NULL, uses Stan default), `max_treedepth` (NULL),
  `seed` (123), `parallel_chains` (NULL), `max_param_count` (NULL – set
  to a positive integer to stop if approximate parameter count exceeds
  the threshold). Any additional entries are forwarded via `...`.

- show_messages:

  Logical. Print sampling progress. Default `TRUE`.

- ...:

  Additional arguments forwarded to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

Named list with elements: `draws`, `diagnostics`, `diagnostics_detail`,
`fit`, `data_long`, `index_maps`, `X_design`, `class_cols`,
`event_metadata`, `n_re_levels`, `upper_re_col`, `middle_re_col`,
`lower_re_col`, `pathogen_col`, `pathogen_fitted`, `estimand`,
`prior_config_used`, `sampler_config_used`, `eligibility_report`.

`diagnostics` is a one-row tibble reported in **two scopes**, because
the model has `N_events * D` latent `z_free` nuisance parameters from
probit data augmentation. These latent variables are excluded from
`draws`, `draws_summary.csv`, and
[`plot_probit_diagnostics()`](https://saketlab.github.io/anumaan/reference/plot_probit_diagnostics.md),
but are still part of the Stan fit:

- `max_rhat_structural`, `min_ess_bulk_structural`,
  `min_ess_tail_structural`, `converged_structural` – computed only over
  the retained structural parameters such as `beta`, random-effect
  terms, `tau_*`, `R_*`, `Omega` when present, and `lp__`. This is the
  scope that matches `draws` and diagnostic plots, and is the
  recommended pass/fail signal for the resistance-profile model.

- `max_rhat_full`, `min_ess_bulk_full`, `min_ess_tail_full`,
  `converged_full` – the same computation over every Stan parameter,
  including `z_free`. A small number of `z_free` entries crossing the
  Rhat 1.01 or ESS 100 thresholds can occur even when the structural
  parameters are well behaved. Therefore, `converged_full = FALSE` while
  `converged_structural = TRUE` should not by itself be read as model
  failure.

- `latent_diagnostic_warning` – `TRUE` when `converged_structural` is
  `TRUE` but `converged_full` is `FALSE`. This is informational
  metadata, not a failure signal.

- `n_divergent`, `n_treedepth_sat`, `ebfmi_min`, `ebfmi_mean`, and
  `ebfmi_by_chain` – sampler-health diagnostics, not parameter-scope
  dependent.

`diagnostics_detail` contains structural-parameter, full-parameter, and
chain-level diagnostic tables.

## Details

**Model (per pathogen):** \$\$Y\_{ed} = \mathbf{1}(Z\_{ed} \> 0)\$\$
\$\$Z\_{ed} = \mathbf{x}\_e^\top\beta_d +
\text{hospital\\effect}\_{d,h(e)} \\\[+\\
\text{patient\\effect}\_{d,p(e)}\] \\\[+\\
\text{admission\\effect}\_{d,a(e)}\] + (L\_\Omega\\\eta_e)\_d\$\$
\$\$\eta_e \sim N(0, I_D),\quad L\_\Omega \sim
\text{LKJCholesky}(\eta)\$\$

The per-event latent noise \\\eta_e\\ makes \\L\_\Omega\\ (and therefore
\\\Omega = L\_\Omega L\_\Omega^\top\\) identifiable from the observed
binary outcomes. Hospital (and lower-level) effects use a
\\\text{diag}(\tau) L\_{\text{corr}} z\_{\text{raw}}\\ parameterisation
so the D-dim random effect vectors are correlated across antibiotic
classes.

**Single-pathogen design:** pass `pathogen` to restrict the fit. Run
once per pathogen and orchestrate in the analysis repository.

**Random effects** (`random_effects`): 1, 2, or 3 grouping column names
defining the clustering hierarchy (e.g. hospital; hospital + patient;
hospital + patient + admission). The first element is the upper-most
level; subsequent elements are nested within it. Nested levels receive
globally unique composite keys built internally. Any hierarchical
grouping variable can occupy any slot – the labels
hospital/patient/admission are semantic conventions, not constraints.
Isolate- or sample-event-level effects can be passed as the second or
third element, but note the returned object uses generic level names
(`upper_re_col`, `middle_re_col`, `lower_re_col`).

**Missing-AST handling:** Observed AST cells impose sign constraints on
the latent variables; untested cells impose no sign constraint and are
represented as unconstrained latent values. The model assumes testing is
conditionally ignorable given included covariates and random effects. It
does NOT correct for selective testing bias. See `residual_structure`
for control over residual covariance complexity.

**Fixed-effect missingness:** the function warns but does NOT silently
impute covariates. Imputation is the caller's responsibility; columns
with any remaining NA values after the call will cause Stan to fail.
