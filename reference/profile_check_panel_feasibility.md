# Pathogen x Candidate-Panel x Hospital Feasibility Table

Builds the feasibility table the way
[`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md)
says it has to be built: the `n_min` benchmark is computed ONCE per
candidate panel (shared across strata), and checked against each
stratum's OWN complete-event count for that exact panel. There is no
aggregation step, no coverage-fraction cutoff, and no auto-selected
"best" panel size – this reports `meets_benchmark` per
pathogen-panel-stratum row and leaves the choice of how many classes to
use, per pathogen, to be made by inspecting that table (together with
[`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)'s
Wilson/Goodman widths and clinical judgment about which classes matter),
not by this function.

## Usage

``` r
profile_check_panel_feasibility(
  data,
  group_col,
  candidate_panels,
  class_col = NULL,
  value_col = "antibiotic_value",
  resistant_values = c("R", "Resistant", "resistant"),
  susceptible_values = c("S", "Susceptible", "susceptible"),
  intermediate_values = c("I", "Intermediate", "intermediate"),
  event_id_col = NULL,
  epsilon,
  alpha = 0.05,
  simultaneous_across_groups = FALSE,
  unknown_action = c("error", "treat_as_untested"),
  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins",
    "first")
)
```

## Arguments

- data, group_col, class_col, value_col, resistant_values,
  susceptible_values, intermediate_values, event_id_col:

  Passed through to
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)
  once per candidate panel (its `class_cols` argument is the candidate
  panel itself).

- candidate_panels:

  A named list of character vectors, each a candidate set of classes to
  evaluate (e.g.
  `list(core4 = c("Carbapenems","Fluoroquinolones","Penicillins", "Aminoglycosides"), alt4 = c(...))`).
  Unnamed entries are labelled `panel_1`, `panel_2`, ...

- epsilon, alpha:

  Passed to
  [`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md).

- simultaneous_across_groups:

  Logical. When `TRUE`, applies the Bonferroni correction across the
  number of groups found for each panel (a joint guarantee across all of
  them). Default `FALSE` (each group checked against its own `1 - alpha`
  guarantee).

- unknown_action, conflict_action:

  Same meaning as in
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
  passed through once per candidate panel – the same AST policy is
  applied to every panel being screened.

## Value

A data frame, one row per candidate panel x group: `panel`, `classes`,
`D`, `K`, `group`, `n_complete`, `n_min`, `meets_benchmark`.

## Details

`candidate_panels` is supplied explicitly, never auto-generated: a panel
of a given size is not interchangeable with any other same-size panel
(A-B-C-D and A-B-D-E can have very different joint support even at the
same `D`), so which combinations to test is a decision for the caller –
a clinically pre-screened shortlist, or an exhaustive enumeration the
caller builds and passes in.
