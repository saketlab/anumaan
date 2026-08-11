# Compute AMR-specific posterior predictive discrepancy statistics

Drives `ppc_draws$generate_state(s)` once per posterior state (see
[`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md)),
streaming each state's complete and observation-masked replicate through
every requested discrepancy statistic and accumulating only the compact
per-state scalar contributions – never an `S x N_events x D` array.
Statistic families (Stan User's Guide-style posterior predictive checks,
adapted to AMR resistance data):

- marginal:

  Per-class resistance rate, pooled and per eligible hospital x pathogen
  x class cell (reusing `fitted_model$eligibility_report$marginal`).

- resistant_count:

  Number of resistant classes per event, among whatever was tested for
  that event, summarised (mean/variance/median/ p90/max); plus, for
  fully-observed eligible panels, the full count distribution.

- pairwise:

  RR/RS/SR/SS co-resistance proportions for every sufficiently co-tested
  class pair (reusing `fitted_model$eligibility_report$pairwise`).

- complete_profile:

  All-susceptible/all-resistant/most-common-profile frequency,
  distinct-profile count, Shannon entropy, Simpson concentration, and
  the top-3 most frequent OBSERVED profiles
  (`profile_top_observed_frequency` – a purely data-dominant descriptive
  summary, deliberately NOT a claim of clinical importance; see the file
  header of `probit_posterior_predictive.R` for the distinction and the
  deferred `clinical_priority_profiles` extension point), restricted to
  adequately-supported fully-observed hospital x pathogen panels
  ([`enumerate_binary_profiles`](https://saketlab.github.io/anumaan/reference/enumerate_binary_profiles.md)-
  style class panel via
  [`.resolve_profile_class_panel()`](https://saketlab.github.io/anumaan/reference/dot-resolve_profile_class_panel.md)).

- hospital_heterogeneity:

  SD/IQR/MAD/range of per-hospital resistance proportions among
  hospitals meeting `min_hospital_support`, per class.

- admission_clustering, patient_clustering:

  Within-/between-group same-class agreement and Hamming distance, only
  when the relevant role is declared eligible (see `random_effect_roles`
  and `random_effect_eligibility_audit` below); skipped cleanly
  (`support_status != "supported"`), never as a failure, when
  unsupported.

## Usage

``` r
compute_probit_ppc_statistics(
  fitted_model_or_ppc_draws,
  n_states = 500L,
  seed = 123L,
  preserve_observation_mask = TRUE,
  statistics = c("marginal", "resistant_count", "pairwise", "complete_profile",
    "hospital_heterogeneity", "admission_clustering", "patient_clustering"),
  random_effect_roles = NULL,
  random_effect_eligibility_audit = NULL,
  min_hospital_support = 30L,
  min_complete_profile_events = 30L,
  ci_level = 0.95
)
```

## Arguments

- fitted_model_or_ppc_draws:

  Either a fitted model object (as for
  [`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md))
  or an `"amr_ppc_draws"` object already returned by that function – the
  latter avoids re-deriving posterior draw setup when statistics are
  recomputed for the same simulated replicates.

- n_states, seed, preserve_observation_mask:

  See
  [`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md).
  Ignored when `fitted_model_or_ppc_draws` is already an
  `"amr_ppc_draws"` object (its own settings are used).

- statistics:

  Character vector of statistic families to compute. Default: all of
  `"marginal"`, `"resistant_count"`, `"pairwise"`, `"complete_profile"`,
  `"hospital_heterogeneity"`, `"admission_clustering"`,
  `"patient_clustering"`.

- random_effect_roles:

  Optional named character vector/list mapping the name of a GROUPING
  COLUMN present in `fitted_model$event_metadata` (the NAME of each
  entry) to one of `"hospital"`, `"admission"`, `"patient"` (the VALUE –
  a closed set of ROLE labels, never an assumption about what any
  particular experiment's columns are actually named), e.g.
  `c(admission_id = "admission")`. Deliberately decoupled from whether
  that column is declared as a random-effect block in this particular
  fit – e.g. tagging a column that is NOT modelled as a random intercept
  with role `"admission"` is exactly how a fit lacking an admission
  random effect gets checked for whether it fails to reproduce observed
  within-admission clustering (see Stan User's Guide- style synthetic
  recovery tests in `test-probit-predictive-synthetic-recovery.R`).
  `NULL` (default) means no admission/patient-role clustering statistics
  can run (reported with `support_status = "role_not_mapped"`).

- random_effect_eligibility_audit:

  Optional data frame with columns `block` (role:
  `"hospital"`/`"admission"`/`"patient"`) and `status`
  (`"eligible_primary"`/`"eligible_sensitivity_only"`/ `"not_eligible"`)
  – e.g. the pathogen-specific random-effect eligibility audit produced
  by the analysis layer. When `NULL`, a generic structural fallback
  check is used instead (see
  [`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md)'s
  file header for the general never-hardcode-eligibility-in-the-package
  philosophy this mirrors).

- min_hospital_support:

  Integer; minimum tested events for a hospital to enter the
  hospital-heterogeneity comparison. Default 30.

- min_complete_profile_events:

  Integer; minimum fully-observed events for a hospital x pathogen panel
  to enter complete-profile / full resistant-count-distribution
  statistics. Default 30.

- ci_level:

  Numeric; determines the `replicated_q025`/ `replicated_q975` quantile
  levels (`(1-ci_level)/2` and its complement). Default 0.95.

## Value

A tibble with columns `statistic_name`, `stratum`, `observed_value`,
`replicated_mean`, `replicated_sd`, `replicated_q025`, `replicated_q50`,
`replicated_q975`, `ppc_tail_probability`, `ppc_two_sided`,
`n_replications`, `support_status`.

## References

Stan Development Team. "Posterior and Prior Predictive Checks." Stan
User's Guide.
<https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html>
