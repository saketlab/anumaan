# Changelog

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
