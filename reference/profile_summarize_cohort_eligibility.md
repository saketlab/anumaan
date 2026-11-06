# Per-Group Cohort Eligibility Summary (One Report, Every Eligibility Fact)

Assembles, in one call, the single per-group eligibility report a
pipeline review actually wants to read: event/patient counts, the
complete/partial/zero-observed breakdown, BOTH independent eligibility
verdicts this family computes, and the full class/pair/higher-order
support tables – instead of a caller re-deriving this by hand from five
separate calls. This is purely a REPORTING layer: it calls
[`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
[`profile_audit_interaction_support`](https://saketlab.github.io/anumaan/reference/profile_audit_interaction_support.md),
[`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md),
[`profile_define_cohort`](https://saketlab.github.io/anumaan/reference/profile_define_cohort.md),
and
[`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md)
once each, on the SAME already-cohort-filtered `data` – it computes
nothing new beyond assembling their outputs into one table.

## Usage

``` r
profile_summarize_cohort_eligibility(
  data,
  group_col,
  panel,
  class_col = NULL,
  value_col = "antibiotic_value",
  resistant_values = c("R", "Resistant", "resistant"),
  susceptible_values = c("S", "Susceptible", "susceptible"),
  intermediate_values = c("I", "Intermediate", "intermediate"),
  intermediate_as_resistant = TRUE,
  event_id_col,
  patient_id_col = NULL,
  epsilon,
  alpha = 0.05,
  tiers = c(tier1 = 5, tier2 = 7.5, tier3 = 10),
  min_tier = "tier2",
  estimands = c("class", "pair", "profile"),
  conf_level = 0.95,
  max_interaction_order = NULL,
  min_co_tested = 5L,
  unknown_action = c("error", "treat_as_untested"),
  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins",
    "first")
)
```

## Arguments

- data:

  Data frame, already restricted to the modeling cohort (e.g.
  `profile_select_cohort(...)$data`) – every count here describes THIS
  cohort, never the full unfiltered dataset.

- group_col, class_col, value_col, resistant_values, susceptible_values,
  intermediate_values, intermediate_as_resistant, unknown_action,
  conflict_action:

  Same meaning as in
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md).

- panel:

  Character vector. The candidate panel's classes (what
  [`profile_check_panel_feasibility`](https://saketlab.github.io/anumaan/reference/profile_check_panel_feasibility.md)
  calls one entry of `candidate_panels`, and
  [`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md)
  calls `panel`).

- event_id_col:

  Character. Required (no default, unlike elsewhere in this family) –
  event counting is this function's whole job, and a fallback row-index
  id would silently miscount events for long-format data (one row per
  class per event, not one row per event).

- patient_id_col:

  Character or `NULL`. Column identifying the patient. When supplied,
  `n_patients` is the count of DISTINCT patients per group; when `NULL`
  (default), `n_patients` is `NA_integer_` – no other function in this
  family tracks patient identity, so it is never inferred from something
  else.

- epsilon:

  Numeric in `(0,1)`. Passed to
  [`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md)
  for `direct_support_pass`.

- alpha:

  Numeric in `(0,1)`. Passed to
  [`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md).
  Default `0.05`.

- tiers, min_tier, estimands, conf_level:

  Passed to
  [`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)/[`profile_define_cohort`](https://saketlab.github.io/anumaan/reference/profile_define_cohort.md)
  for `included_in_fit`. Defaults match those functions' own defaults.

- max_interaction_order:

  Integer or `NULL`. Passed to
  [`profile_audit_interaction_support`](https://saketlab.github.io/anumaan/reference/profile_audit_interaction_support.md)
  as `max_order`. Ignored (that table is returned empty) when
  `length(panel) < 2`.

- min_co_tested:

  Integer. Passed to
  [`profile_audit_interaction_support`](https://saketlab.github.io/anumaan/reference/profile_audit_interaction_support.md).
  Default `5L`.

## Value

A list:

- `summary`:

  one row per group: `group`, `n_pathogen_events`, `n_patients`,
  `n_complete`, `n_partial`, `n_zero_observed`, `direct_support_pass`,
  `included_in_fit`

- `class`:

  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
  `$class` table, on this cohort/panel

- `pair`:

  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
  `$pair` table, on this cohort/panel

- `interaction`:

  [`profile_audit_interaction_support`](https://saketlab.github.io/anumaan/reference/profile_audit_interaction_support.md)'s
  table (higher-order support); empty data frame if `length(panel) < 2`

## Details

The two eligibility verdicts are genuinely independent, not duplicates
of each other – see
[`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md)'s
docs for why neither subsumes the other:

- `direct_support_pass`:

  this group's `n_complete` meets
  [`profile_precision_bound`](https://saketlab.github.io/anumaan/reference/profile_precision_bound.md)'s
  distribution-free, Hoeffding/union-bound `n_min` for this exact panel.

- `included_in_fit`:

  this group clears `min_tier` on every requested estimand, per the
  Wilson/Goodman precision-tier roster decision in
  [`profile_define_cohort`](https://saketlab.github.io/anumaan/reference/profile_define_cohort.md).

A group can pass one and fail the other.

`n_zero_observed` needs care: the shared pivot this family uses
internally only produces a row for events with AT LEAST ONE tested panel
class, so a naive count over that pivot would silently miss
zero-observed events entirely rather than counting them. This function
reconstructs the full event universe from `data` itself (via
`group_col` + `event_id_col`) and left-joins the tested-class counts
onto it, so an event with nothing tested on this panel is correctly
tallied as `n_zero_observed`, not silently dropped.
