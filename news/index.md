# Changelog

## anumaan 0.1.0.9045

- [`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md)
  has been rewritten. It now returns a list (`rr`, `hospital`, `draws`,
  `skipped`) instead of a data frame.
  - Hospitals enter as a fixed effect in one model per pathogen x class
    (`hospital_effect = "fixed"`, the new default), so the relative LOS
    is shared and every hospital keeps its own LOS level. `"none"` and
    `"stratified"` (one model per hospital, the previous behaviour)
    remain available. `RxHospital_p` tests whether the resistance effect
    differs between hospitals.
  - The unit is the admission (`admission_id_col`; `patient_id_col` is
    deprecated). Resistance is the worst phenotype of the pathogen’s own
    results in the admission, and covariates come from its earliest
    culture of that pathogen. Polymicrobial admissions now enter every
    pathogen’s models.
  - Returns adjusted `LOSR` and `LOSS` (model predictions averaged over
    each hospital’s resistant admissions with R = 1 and R = 0) and their
    difference, per hospital and overall, instead of raw means only.
  - Confidence intervals for `RR_LOS`, `LOSR`, `LOSS` and the excess
    from a bootstrap that resamples clusters (`cluster_id_col`,
    e.g. patients) within hospitals (`ci_method = "bootstrap"`), or from
    the delta method with a cluster-robust covariance (`"delta"`).
    `return_draws = TRUE` returns the bootstrap replicates.
  - New arguments `classes` (a panel per pathogen),
    `reference_hospital`, `intermediate_as_resistant`, `hai_value`,
    `cai_value`, `icu_values`, `flag_min`. `hai_threshold_hours` is
    deprecated.
- [`daly_assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/daly_assign_rr_to_profiles.md)
  accepts profiles as R/S strings (e.g. the `profile_delta` / `R_mean`
  output of
  [`profile_aggregate_for_daly()`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md))
  with `profile_class_order` giving the class of each letter, and
  accepts the
  [`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md)
  list directly. Each profile then takes the adjusted LOS of its
  dominant class at its own hospital. Ties in the max rule go to the
  first class in class order and are reported.
- [`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md)
  and
  [`daly_calc_yld_attributable()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_attributable.md)
  match each hospital to its own profiles and LOS when profiles are
  hospital-specific, and return percentile intervals when given
  `los_draws` (bootstrap draws from
  [`daly_fit_los_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_rr.md))
  and optionally `profile_draws` (profile posterior draws).
- [`daly_fit_mortality_rr()`](https://saketlab.github.io/anumaan/reference/daly_fit_mortality_rr.md)
  now builds class resistance from each pathogen’s own results.
  Previously one pathogen’s results could set another pathogen’s
  resistance in the same patient, and patients with several pathogens
  entered only one pathogen’s models.

## anumaan 0.1.0.9044

- [`prep_standardize_antibiotics()`](https://saketlab.github.io/anumaan/reference/prep_standardize_antibiotics.md)
  now keeps combination antibiotics intact. Names such as
  “piperacillin-tazobactam” and “trimethoprim-sulfamethoxazole” were
  reduced to one component and given that component’s class
  (“Penicillins”, “Sulfonamides”). They now match their WHO combination
  rows (“Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal”,
  “Sulfonamide-trimethoprim-combinations”), and reserve combinations
  such as ceftazidime-avibactam keep their AWaRe category. Matching
  treats `-`, `+`, `/`, `&` and “and” as one separator, ignores word
  order and tolerates one-letter misspellings in long names.
- A combination with no reference row keeps its full name, is left
  unclassified and is reported in a warning instead of being reduced to
  one component.
- Added `inst/extdata/antibiotic_combination_supplement.csv` for
  combinations absent from `WHO_aware_class.csv`
  (cefoperazone-sulbactam, ticarcillin-clavulanic acid). It is appended
  to the default reference.

## anumaan 0.1.0.9042

- Added the joint Bernoulli resistance-profile pipeline (`profile_*`), a
  log-linear model over all 2^D resistance states that marginalises
  untested antibiotic classes exactly. It covers cohort selection and
  support audits
  ([`profile_select_cohort()`](https://saketlab.github.io/anumaan/reference/profile_select_cohort.md),
  [`profile_audit_support()`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
  [`profile_plan_precision()`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md),
  [`profile_define_cohort()`](https://saketlab.github.io/anumaan/reference/profile_define_cohort.md)),
  model design and fitting
  ([`profile_build_design()`](https://saketlab.github.io/anumaan/reference/profile_build_design.md),
  [`profile_fit_joint_bernoulli()`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)),
  diagnostics
  ([`profile_check_sampler()`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md),
  [`profile_ppc()`](https://saketlab.github.io/anumaan/reference/profile_ppc.md),
  [`profile_prior_predictive()`](https://saketlab.github.io/anumaan/reference/profile_prior_predictive.md)),
  validation
  ([`profile_validate_masked()`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md),
  [`profile_validate_heldout()`](https://saketlab.github.io/anumaan/reference/profile_validate_heldout.md),
  [`profile_validate_newpatient()`](https://saketlab.github.io/anumaan/reference/profile_validate_newpatient.md),
  [`profile_validate_temporal()`](https://saketlab.github.io/anumaan/reference/profile_validate_temporal.md)),
  estimation
  ([`profile_estimate()`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)),
  aggregation
  ([`profile_aggregate_for_daly()`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md))
  and reporting
  ([`profile_release()`](https://saketlab.github.io/anumaan/reference/profile_release.md),
  [`profile_render_report()`](https://saketlab.github.io/anumaan/reference/profile_render_report.md)).
- Added
  [`profile_audit_profile_distribution()`](https://saketlab.github.io/anumaan/reference/profile_audit_profile_distribution.md)
  for per-profile prevalence auditing and
  [`profile_aggregate_counts_with_allocation()`](https://saketlab.github.io/anumaan/reference/profile_aggregate_counts_with_allocation.md)
  for discrete count aggregation.
- Patient holdout validation now supports stratified sampling and scores
  complete joint profiles. Sampling is fixed for single numeric patient
  IDs, singleton strata stay in training, and stratum holdout is capped.
- Masked validation no longer leaves an event with no observed classes,
  and average precision no longer depends on the order of tied scores.
- Estimation and masked-cell prediction share one compatible-state
  softmax, which prevents underflow in conditional profile
  probabilities.
- Validation refits now receive the original model specification through
  `design_spec`; unknown `design_spec` fields are an error.
- The burden eligibility flags from
  [`profile_aggregate_for_daly()`](https://saketlab.github.io/anumaan/reference/profile_aggregate_for_daly.md)
  are renamed `computational_checks_passed_YLL` and
  `computational_checks_passed_YLD`, so they are not read as a release
  decision.
- [`prep_deduplicate_ast()`](https://saketlab.github.io/anumaan/reference/prep_deduplicate_ast.md)
  is generalised, and the `profile_*` functions delegate their AST
  conflict policy to it.
- Length of stay now has a single definition. Distribution-based LOS
  uses medians, and
  [`daly_fit_los_distn_by_profile()`](https://saketlab.github.io/anumaan/reference/daly_fit_los_distn_by_profile.md)
  fits LOS at profile level. `prepare_los_data()` is removed.
- Removed unused multivariate-probit utilities:
  `plot_probit_diagnostics()`,
  `plot_probit_posterior_predictive_checks()`,
  `compute_probit_ppc_statistics()`,
  `compute_posterior_predictive_status()`,
  `compute_prior_predictive_status()`,
  `compute_prior_predictive_fingerprint()` and the
  `simulate_probit_*_predictive()` functions.
- [`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md)
  moved from the plotting file to the YLL file; its behaviour is
  unchanged.

## anumaan 0.1.0.9036

- Prepared for CRAN submission: every exported example now runs,
  internal helpers no longer produce Rd pages, and console tables
  printed during preprocessing are emitted with
  [`message()`](https://rdrr.io/r/base/message.html) so
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) silences
  them.
- Pathway 1 now works with osqp \>= 1.0. osqp 1.0 removed the `$solve()`
  method, so every osqp solve errored and
  [`estimate_profiles_convex()`](https://saketlab.github.io/anumaan/reference/estimate_profiles_convex.md),
  [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
  and
  [`bootstrap_profiles_convex()`](https://saketlab.github.io/anumaan/reference/bootstrap_profiles_convex.md)
  silently returned uniform profile probabilities flagged as
  non-converged. They now call
  [`osqp::solve_osqp()`](https://rdrr.io/pkg/osqp/man/solve_osqp.html),
  which works across osqp versions.
- [`compute_marginals_from_data()`](https://saketlab.github.io/anumaan/reference/compute_marginals_from_data.md)
  and
  [`compute_pairwise_from_data()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_from_data.md)
  no longer drop antibiotic classes whose names contain hyphens, slashes
  or spaces (e.g. “Third-generation-cephalosporins”), which
  [`preprocess_for_profiles()`](https://saketlab.github.io/anumaan/reference/preprocess_for_profiles.md)
  had renamed when pivoting to wide format.
- Organism, drug and syndrome names are now compared case-insensitively
  through one lower-case key. As a result `org_group` and `rr_pathogen`
  are populated for
  [`run_preprocess()`](https://saketlab.github.io/anumaan/reference/run_preprocess.md)
  output (previously “Other” and `NA` for every row), and `panel_map`,
  `syndrome_name` and `pathogen_name` arguments match regardless of
  case; outputs keep the data’s own labels. The derived keys
  `rr_pathogen` and the `pathogen` column of
  [`daly_load_rr_reference()`](https://saketlab.github.io/anumaan/reference/daly_load_rr_reference.md)
  are now lower case.
- Fixed `inst/extdata/adjustment_for_CFR.csv`: an unquoted comma in
  “Infections of bones, joints, and related organs” corrupted that row.
- Functions that take a `seed` argument now restore the caller’s random
  number generator state on exit.
- [`preprocess_for_profiles()`](https://saketlab.github.io/anumaan/reference/preprocess_for_profiles.md)
  no longer errors at the deduplication step; it now deduplicates
  isolate x drug rows before collapsing to antibiotic class level.
- [`daly_calc_yll_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yll_associated.md)
  and
  [`daly_load_rr_reference()`](https://saketlab.github.io/anumaan/reference/daly_load_rr_reference.md)
  now read the bundled CSV reference tables (the `.xlsx` files they
  expected are no longer shipped).

## anumaan 0.1.0.9021

- Expanded the specimen classification reference with stewardship sample
  labels seen in 2021-era extracts, including brain abscess, instrument,
  lung aspirate, lymph node, hair, and superficial biopsy.
- Collapsed `sterile_classification` outputs from detailed reference
  labels to `Sterile`, `Non-Sterile`, and `Others/Ambiguous` for
  downstream analysis.
- Made EDA and burden plotting functions usable without `ggpubr`; when
  `ggpubr` is not installed,
  [`eda_theme()`](https://saketlab.github.io/anumaan/reference/eda_theme.md)
  falls back to
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).
- Refreshed vignettes and installation guidance for the 0.1.0.9021
  workflow.
