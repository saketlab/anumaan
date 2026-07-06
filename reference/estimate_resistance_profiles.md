# Estimate Resistance Profiles: Pathway 1 (Convex) or Pathway 2 (Bayesian)

Top-level dispatcher. Runs either the convex optimisation pathway
(Pathway 1, aggregate surveillance data) or the Bayesian hierarchical
multivariate probit pathway (Pathway 2, facility-level AST data).

## Usage

``` r
estimate_resistance_profiles(
  data,
  method = c("convex", "bayesian"),
  pairwise = NULL,
  panel_map = NULL,
  class_cols = NULL,
  fixed_effects = NULL,
  random_effects = NULL,
  pathogen = NULL,
  pathogen_col = "pathogen",
  eligible_pairs = NULL,
  outcome_col = NULL,
  nonfatal_values = c("Discharged", "Survived", "Alive", "discharged", "survived",
    "alive"),
  panel_eligibility = list(),
  residual_structure = c("identity", "correlated"),
  estimand = "observed_stewardship_event_mix",
  prior_config = list(),
  sampler_config = list(),
  n_posterior_draws_for_profiles = 2000L
)
```

## Arguments

- data:

  For Pathway 1: marginals tibble from
  [`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md).
  For Pathway 2: wide-format event-level tibble (one row per event;
  class columns 0/1/NA).

- method:

  Character. `"convex"` or `"bayesian"`.

- pairwise:

  Tibble or `NULL`. Pathway 1 only.

- panel_map:

  Named list or `NULL`. Pathway 1 only.

- class_cols:

  Character vector. Pathway 2 only. Required.

- fixed_effects:

  Character vector. Pathway 2 only. Required.

- random_effects:

  Character vector (1-3 elements). Pathway 2 only. Required. Elements:
  hospital; \[+patient\]; \[+admission\].

- pathogen:

  Character or `NULL`. Pathway 2 only. When supplied, filters data to a
  single pathogen before fitting.

- pathogen_col:

  Character. Default `"pathogen"`.

- eligible_pairs:

  Tibble or `NULL`. Pathway 2 only.

- outcome_col:

  Character or `NULL`. Pathway 2 only. Patient outcome column for R_ALL
  / R_NF split.

- nonfatal_values:

  Character vector. Outcome values for R_NF cohort.

- panel_eligibility:

  Named list. Eligibility thresholds; passed to
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- residual_structure:

  Character. Pathway 2. `"identity"` (default) or `"correlated"`. See
  [`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).

- estimand:

  Character. Estimand label. Default `"observed_stewardship_event_mix"`.

- prior_config:

  Named list. Pathway 2 only. Priors.

- sampler_config:

  Named list. Pathway 2 only. Structured sampler settings: `chains` (4),
  `iter_warmup` (1000), `iter_sampling` (1000), `seed` (123),
  `adapt_delta`, `max_treedepth`, `parallel_chains`.

- n_posterior_draws_for_profiles:

  Integer. Pathway 2. Number of posterior draws used for MVN simulation
  of profiles. Default `2000L`.

## Value

Named list: `profiles`, `eligibility`, `diagnostics`, `fitted_models`,
`config_used`.

## Details

For Pathway 2, pass `pathogen` to fit one pathogen at a time – the
recommended workflow. Run this function once per pathogen and collect
the results in the analysis repository.
