# Build the Wide Model Dataset for a Chosen Panel

Stage 3 of the resistance-profiling pipeline: turns long- or wide-format
AST data, restricted to one chosen panel of antibiotic classes, into the
`Y` matrix
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
expects (events x classes, entries in `c(0, 1, NA)`), plus everything
needed to trace and model it: group/event identifiers, a covariate
table, and an exclusion manifest for events with nothing in the panel
tested at all.

## Usage

``` r
profile_prepare_model_data(
  data,
  group_col,
  panel,
  class_col = NULL,
  value_col = "antibiotic_value",
  resistant_values = c("R", "Resistant", "resistant"),
  susceptible_values = c("S", "Susceptible", "susceptible"),
  intermediate_values = c("I", "Intermediate", "intermediate"),
  intermediate_as_resistant = TRUE,
  event_id_col = NULL,
  covariates = NULL,
  unknown_action = c("error", "treat_as_untested"),
  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins",
    "first"),
  covariate_conflict_action = c("error", "first")
)
```

## Arguments

- data, group_col, class_col, value_col, resistant_values,
  susceptible_values, intermediate_values, event_id_col:

  Same as
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md).
  `event_id_col` is required (not optional) for long-format `data` – a
  real per-event key is needed to merge covariates back onto the wide
  outcome matrix.

- panel:

  Character vector. The chosen panel's classes – becomes `colnames(Y)`,
  in this order.

- intermediate_as_resistant:

  Logical. Same meaning as in
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md);
  controls whether `Y` folds `"I"` into `1` (resistant) or leaves it
  untested (`NA`).

- covariates:

  Character vector of column names in `data` to carry into the returned
  covariate table, one value per event. Rows sharing an `event_id_col`
  value are checked for agreement on every requested covariate first;
  see `covariate_conflict_action`. Default `NULL` (no covariates carried
  – just `Y` and identifiers).

- unknown_action, conflict_action:

  Same meaning as in
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)
  – the shared AST policy for susceptibility RESULTS, applied
  identically here.

- covariate_conflict_action:

  Character, `"error"` (default) or `"first"`. A DIFFERENT conflict from
  `conflict_action` above: this one is about COVARIATE values (e.g. age,
  ward) disagreeing across the several raw rows that share one
  `event_id_col` value, not about AST results. `"error"` stops, naming
  every event id and covariate where rows disagree – silently keeping
  whichever row happened to come first is exactly the failure mode this
  guards against. `"first"` restores the old behavior (take the first
  row's value) for callers who have already verified, outside this
  function, that their covariates are event-invariant.

## Value

A list:

- `Y`:

  matrix, events x `length(panel)`, entries in `c(0, 1, NA)`, column
  names = `panel`

- `ids`:

  data frame, one row per event matching `Y`'s rows: `event_id`, `group`

- `covariates`:

  data frame, one row per event matching `Y`'s rows, the requested
  `covariates` columns (zero columns if none requested)

- `status`:

  character vector, one per event: `"complete"` (every panel class
  tested) or `"partial"`

- `manifest`:

  one-row data frame: `n_events_input`, `n_events_excluded_untested` (no
  panel class tested at all), `n_events_output` (`= nrow(Y)`),
  `n_complete`, `n_partial`

## Details

Covariates are returned as a plain table, not a model-matrix design –
identifiers (patient id, hospital, admission id) are carried for
grouping and provenance only and do not automatically become predictors
or random effects. Build the actual `X` passed to
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
from this table with
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) and
whatever formula the model specification calls for; that choice belongs
to the caller (the analysis repository's model config), not to this
function.

Every event counted as model-eligible here is the set that
[`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)/[`profile_check_panel_feasibility`](https://saketlab.github.io/anumaan/reference/profile_check_panel_feasibility.md)
should be RE-RUN against once this function has dropped events with
nothing tested in the panel at all (or, if a covariate's completeness is
later enforced downstream, re-run again on whatever survives that too) –
support computed only on the pre-cohort-selection or pre-covariate-
reconciliation data describes a bigger cohort than the one that actually
reaches the model.
