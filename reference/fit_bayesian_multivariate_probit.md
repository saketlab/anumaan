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
  random_effects = list(),
  profile_group_col = NULL,
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
  compute = list(),
  show_messages = TRUE,
  save_full_latent_diagnostics = FALSE,
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

  Character vector or named list of random-intercept blocks. Use
  [`list()`](https://rdrr.io/r/base/list.html) for no random-effect
  blocks. When non-empty, the legacy character-vector form names
  grouping columns; the list form uses `name`, `group_col`, and optional
  `terms = "intercept"`.

- profile_group_col:

  Character scalar or `NULL`. Column used for downstream profile
  aggregation, eligibility summaries, and validation. It is independent
  of the random-effect specification. Defaults to the first
  random-effect grouping column when random effects exist; it is
  required for fixed-effects-only fits.

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

- compute:

  Named list controlling the Stan execution backend. Supported fields
  are `backend` (`"cpu"` default or `"opencl"`), `opencl_platform_id`,
  `opencl_device_id`, and `allow_cpu_fallback` (`FALSE` by default).
  OpenCL changes the compilation/sampling backend only; it does not
  change the statistical model, priors, diagnostics, or downstream
  estimands.

- show_messages:

  Logical. Print sampling progress. Default `TRUE`.

- save_full_latent_diagnostics:

  Logical. When `FALSE` (default), `diagnostics_detail$all_parameters`
  (per-parameter Rhat/ESS including the `N_events * D` latent `z_free`
  nuisance parameters – tens of thousands of rows for realistic data) is
  returned with 0 rows to keep the fitted object small; the full-scope
  summary scalars (`max_rhat_full`, `min_ess_bulk_full`,
  `min_ess_tail_full`, `n_params_full`, `pct_rhat_full_above_1_01`,
  `pct_ess_bulk_full_below_100`) are always populated in `diagnostics`
  regardless. Set `TRUE` to retain the full per-parameter table for
  debugging.

- ...:

  Additional arguments forwarded to
  [`cmdstanr::sample()`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html).

## Value

Named list with elements: `draws`, `diagnostics`, `diagnostics_detail`,
`fit`, `data_long`, `index_maps`, `X_event`, `event_re_idx`,
`class_cols`, `event_metadata`, `n_re_levels`, `upper_re_col`,
`middle_re_col`, `lower_re_col`, `patient_key_col`, `admission_key_col`,
`pathogen_col`, `pathogen_fitted`, `residual_structure`, `estimand`,
`prior_config_used`, `sampler_config_used`, `compute_config_used`, and
`eligibility_report`.

`diagnostics` is a one-row monitored summary. The main diagnostic fields
`max_rhat`, `min_ess_bulk`, and `min_ess_tail` are computed over the
monitored parameters retained in `draws`. These exclude the latent
`z_free` nuisance parameters used for probit data augmentation.
`converged_structural`/`converged_full` incorporate tree-depth
saturation as well as R-hat/ESS/divergences/E-BFMI (a run that
repeatedly saturates `max_treedepth` is not reported as converged).
`diagnostic_status` is a canonical multi-level status – one of `"pass"`,
`"warning_rhat"`, `"warning_low_bulk_ess"`, `"warning_low_tail_ess"`,
`"warning_treedepth"`, `"fail_rhat"`, `"fail_energy"`, or
`"fail_divergent"` (most severe first) – intended as the single
authoritative status field for callers, replacing ad hoc re-derivations
downstream.

Full-fit diagnostic fields, including `max_rhat_all`,
`min_ess_bulk_all`, and `min_ess_tail_all`, are reported separately for
visibility and include all Stan parameters, including `z_free`.

`diagnostics_detail` contains `monitored_parameters`, `all_parameters`
(per-parameter Rhat/ESS with `parameter_group`), `grouped`
(per-`parameter_group` Rhat/ESS quantiles and percentage-over-threshold,
not only the single worst parameter), and `chains` (per-chain
`n_sampling`, `n_divergent`, `n_treedepth_sat`, `ebfmi`,
`mean_accept_stat`, `mean_stepsize`, `mean_treedepth`, `max_treedepth`)
– the canonical diagnostic tables; callers should not recompute these
independently from `fit$fit$sampler_diagnostics()`.

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

**Random effects** (`random_effects`): an optional character vector or
list of grouping blocks defining the clustering hierarchy. Use
[`list()`](https://rdrr.io/r/base/list.html) for a fixed-effects-only
model. When blocks are present, the first element is the upper-most
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
