# anumaan 0.1.0.9044

* `prep_standardize_antibiotics()` now keeps combination antibiotics intact.
  Names such as "piperacillin-tazobactam" and "trimethoprim-sulfamethoxazole"
  were reduced to one component and given that component's class
  ("Penicillins", "Sulfonamides"). They now match their WHO combination rows
  ("Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal",
  "Sulfonamide-trimethoprim-combinations"), and reserve combinations such as
  ceftazidime-avibactam keep their AWaRe category. Matching treats `-`, `+`,
  `/`, `&` and "and" as one separator, ignores word order and tolerates
  one-letter misspellings in long names.
* A combination with no reference row keeps its full name, is left
  unclassified and is reported in a warning instead of being reduced to one
  component.
* Added `inst/extdata/antibiotic_combination_supplement.csv` for combinations
  absent from `WHO_aware_class.csv` (cefoperazone-sulbactam,
  ticarcillin-clavulanic acid). It is appended to the default reference.

# anumaan 0.1.0.9042

* Added the joint Bernoulli resistance-profile pipeline (`profile_*`), a
  log-linear model over all 2^D resistance states that marginalises untested
  antibiotic classes exactly. It covers cohort selection and support audits
  (`profile_select_cohort()`, `profile_audit_support()`,
  `profile_plan_precision()`, `profile_define_cohort()`), model design and
  fitting (`profile_build_design()`, `profile_fit_joint_bernoulli()`),
  diagnostics (`profile_check_sampler()`, `profile_ppc()`,
  `profile_prior_predictive()`), validation (`profile_validate_masked()`,
  `profile_validate_heldout()`, `profile_validate_newpatient()`,
  `profile_validate_temporal()`), estimation (`profile_estimate()`),
  aggregation (`profile_aggregate_for_daly()`) and reporting
  (`profile_release()`, `profile_render_report()`).
* Added `profile_audit_profile_distribution()` for per-profile prevalence
  auditing and `profile_aggregate_counts_with_allocation()` for discrete
  count aggregation.
* Patient holdout validation now supports stratified sampling and scores
  complete joint profiles. Sampling is fixed for single numeric patient IDs,
  singleton strata stay in training, and stratum holdout is capped.
* Masked validation no longer leaves an event with no observed classes, and
  average precision no longer depends on the order of tied scores.
* Estimation and masked-cell prediction share one compatible-state softmax,
  which prevents underflow in conditional profile probabilities.
* Validation refits now receive the original model specification through
  `design_spec`; unknown `design_spec` fields are an error.
* The burden eligibility flags from `profile_aggregate_for_daly()` are renamed
  `computational_checks_passed_YLL` and `computational_checks_passed_YLD`, so
  they are not read as a release decision.
* `prep_deduplicate_ast()` is generalised, and the `profile_*` functions
  delegate their AST conflict policy to it.
* Length of stay now has a single definition. Distribution-based LOS uses
  medians, and `daly_fit_los_distn_by_profile()` fits LOS at profile level.
  `prepare_los_data()` is removed.
* Removed unused multivariate-probit utilities:
  `plot_probit_diagnostics()`, `plot_probit_posterior_predictive_checks()`,
  `compute_probit_ppc_statistics()`, `compute_posterior_predictive_status()`,
  `compute_prior_predictive_status()`,
  `compute_prior_predictive_fingerprint()` and the
  `simulate_probit_*_predictive()` functions.
* `compute_hospital_daly()` moved from the plotting file to the YLL file; its
  behaviour is unchanged.

# anumaan 0.1.0.9036

* Prepared for CRAN submission: every exported example now runs, internal
  helpers no longer produce Rd pages, and console tables printed during
  preprocessing are emitted with `message()` so `suppressMessages()` silences
  them.
* Pathway 1 now works with osqp >= 1.0. osqp 1.0 removed the `$solve()`
  method, so every osqp solve errored and `estimate_profiles_convex()`,
  `compute_resistance_profiles()` and `bootstrap_profiles_convex()` silently
  returned uniform profile probabilities flagged as non-converged. They now call
  `osqp::solve_osqp()`, which works across osqp versions.
* `compute_marginals_from_data()` and `compute_pairwise_from_data()` no longer
  drop antibiotic classes whose names contain hyphens, slashes or spaces (e.g.
  "Third-generation-cephalosporins"), which `preprocess_for_profiles()` had
  renamed when pivoting to wide format.
* Organism, drug and syndrome names are now compared case-insensitively through
  one lower-case key. As a result `org_group` and `rr_pathogen` are populated for
  `run_preprocess()` output (previously "Other" and `NA` for every row), and
  `panel_map`, `syndrome_name` and `pathogen_name` arguments match regardless of
  case; outputs keep the data's own labels. The derived keys `rr_pathogen` and
  the `pathogen` column of `daly_load_rr_reference()` are now lower case.
* Fixed `inst/extdata/adjustment_for_CFR.csv`: an unquoted comma in "Infections
  of bones, joints, and related organs" corrupted that row.
* Functions that take a `seed` argument now restore the caller's random number
  generator state on exit.
* `preprocess_for_profiles()` no longer errors at the deduplication step; it now
  deduplicates isolate x drug rows before collapsing to antibiotic class level.
* `daly_calc_yll_associated()` and `daly_load_rr_reference()` now read the
  bundled CSV reference tables (the `.xlsx` files they expected are no longer
  shipped).

# anumaan 0.1.0.9021

* Expanded the specimen classification reference with stewardship sample labels
  seen in 2021-era extracts, including brain abscess, instrument, lung aspirate,
  lymph node, hair, and superficial biopsy.
* Collapsed `sterile_classification` outputs from detailed reference labels to
  `Sterile`, `Non-Sterile`, and `Others/Ambiguous` for downstream analysis.
* Made EDA and burden plotting functions usable without `ggpubr`; when `ggpubr`
  is not installed, `eda_theme()` falls back to `ggplot2::theme_bw()`.
* Refreshed vignettes and installation guidance for the 0.1.0.9021 workflow.
