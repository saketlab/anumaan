# Resolve the antibiotic-class panel used for profile enumeration

Defines the set of classes entering the \\2^D\\ profile enumeration for
one hospital x pathogen pair from the fit-time eligibility report
(`fitted_model$eligibility_report`), rather than "tested at least once".
A class is included only if its marginal
n_tested/n_resistant/n_susceptible thresholds were met at that hospital
for that pathogen – this applies to **every** residual structure,
including identity, so an identity-residual panel is never narrowed for
reasons that only matter for estimating a correlated residual matrix.
For correlated-residual fits specifically, classes that lack sufficient
pairwise co-testing with every other candidate class at that hospital
are additionally dropped iteratively (the class involved in the most
insufficient pairs is dropped first, repeated until all remaining pairs
clear the threshold), since \\\Omega\\ for an under-co-tested pair would
otherwise be prior-dominated.

## Usage

``` r
.resolve_profile_class_panel(
  class_cols,
  hospital,
  pathogen,
  eligibility_report,
  upper_re_col,
  pathogen_col,
  residual_structure = "identity"
)
```

## Value

List with:

- classes:

  Character vector of `class_cols` entering the panel, in `class_cols`
  order.

- excluded:

  Character vector of `class_cols` NOT entering the panel (empty if
  none).

- method:

  Character; `"marginal_only"` (identity residual, or no pairwise report
  available), `"marginal_plus_pairwise"` (correlated residual), or
  `"no_eligibility_report_available"` (defensive fallback for fit
  objects built some other way).

- reason:

  Character or `NA`; `NA` when nothing was excluded, otherwise a
  human-readable breakdown of which classes were excluded and why
  (insufficient marginal support vs. insufficient pairwise co-testing).

## Details

The eligibility report passed in is computed by
[`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md)
on `event_data` *after* all experiment-specific filtering (pathogen
filter, `eligible_pairs` semi-join, all-NA-event drop) – i.e. on the
exact row population that was fitted, not on an earlier unfiltered wide
table. See the eligibility-report construction in
[`fit_bayesian_multivariate_probit()`](https://saketlab.github.io/anumaan/reference/fit_bayesian_multivariate_probit.md).
