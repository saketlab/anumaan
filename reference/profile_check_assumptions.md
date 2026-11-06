# Structured Pre-Fit Assumption Register

Assembles the fixed 12-question assumption register this family's design
review requires (Section C of the implementation spec), auto-filling
whatever rows are genuinely computable from diagnostics ALREADY produced
elsewhere in this package – never re-deriving them – and leaving every
other row `status = "not_evaluated"`. **A statistically untestable
assumption (e.g. whether testing is informatively selective on the
unobserved resistance value, or whether shared interactions are
scientifically plausible) is never marked `"verified_by_data_contract"`
by this function**, no matter what inputs are supplied – only
`correct_unit`, `design_estimable`, and
`direct_support_present`/`priors_appropriate` have any automatic path to
that status, and only when the relevant hard check actually passed.

## Usage

``` r
profile_check_assumptions(
  support = NULL,
  fixed_effect_design = NULL,
  prior_predictive = NULL,
  reviewer_notes = list()
)
```

## Arguments

- support:

  Output of
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
  or `NULL` to leave `direct_support_present` unfilled.

- fixed_effect_design:

  Output of
  [`profile_build_fixed_effect_design`](https://saketlab.github.io/anumaan/reference/profile_build_fixed_effect_design.md),
  or `NULL` to leave `design_estimable` unfilled.

- prior_predictive:

  Output of
  [`profile_prior_predictive`](https://saketlab.github.io/anumaan/reference/profile_prior_predictive.md),
  or `NULL` to leave `priors_appropriate` unfilled.

- reviewer_notes:

  Named list, keyed by `assumption_id` (an error names the offending key
  if any is unrecognized – this catches a typo rather than silently
  dropping the note). Each element is a list with any of `note` (free
  text), `status` (one of `"verified_by_data_contract"`, `"diagnosed"`,
  `"assumed"`, `"not_evaluated"`, `"violated"`), `sensitivity_plan`.
  Applied AFTER the automatic fills, so a reviewer can override an
  automatic verdict with their own judgment (e.g. downgrading
  `"verified_by_data_contract"` to `"assumed"` if they have reservations
  the automatic check cannot see).

## Value

A data frame, one row per `assumption_id`: `claim`, `evidence`,
`diagnostic`, `limitation`, `sensitivity_plan`, `severity`
(`"not_assessed"`, `"low"`, `"moderate"`, `"blocking"`), `status`,
`reviewer_note`.

## Details

Rows this function CAN auto-fill, and from what:

- `correct_unit`:

  from `support` –
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
  `$class$n_duplicates_collapsed` (identical duplicates collapsed, never
  inflating a count). Reported as `"diagnosed"`, not
  `"verified_by_data_contract"` –
  [`profile_audit_support()`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
  default `conflict_action = "error"` means a CONFLICTING duplicate
  would already have stopped whatever call produced `support` (so its
  mere existence is good evidence), but this function cannot see which
  `conflict_action` that call actually used, so it stops short of
  "verified".

- `design_estimable`:

  from `fixed_effect_design` –
  [`profile_build_fixed_effect_design`](https://saketlab.github.io/anumaan/reference/profile_build_fixed_effect_design.md)'s
  return value; checks `rank == ncol(X)`.

- `direct_support_present`:

  from `support` –
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
  return value; flags any zero-count group x class cell. Always
  `status = "diagnosed"`, never `"verified"` – a nonzero count is not
  identifiability or adequate precision (see
  [`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)).

- `priors_appropriate`:

  from `prior_predictive` –
  [`profile_prior_predictive`](https://saketlab.github.io/anumaan/reference/profile_prior_predictive.md)'s
  return value; flags extreme (near-0 or near-max) entropy ratios.
  Always `status = "diagnosed"` – flags only extreme concentration, a
  reviewer still has to judge clinical plausibility of the actual
  implied margins/pair probabilities.

Every other row requires an explicit human judgment, supplied via
`reviewer_notes`.
