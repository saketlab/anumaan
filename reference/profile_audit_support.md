# Tabulate Class-Level Testing Support, by Group

Pure counting – no statistics. For every group (hospital, district,
whatever `group_col` holds) and every class in `class_cols`, counts how
many events were tested, how many were resistant, and how many were
susceptible. Also produces the same counts for every class pair
(co-tested) and for the full class set at once (complete panel), since a
hospital can have plenty of events yet too few co-tested pairs or
complete panels for a joint model.

## Usage

``` r
profile_audit_support(
  data,
  group_col,
  class_cols,
  class_col = NULL,
  value_col = "antibiotic_value",
  resistant_values = c("R", "Resistant", "resistant"),
  susceptible_values = c("S", "Susceptible", "susceptible"),
  intermediate_values = c("I", "Intermediate", "intermediate"),
  intermediate_as_resistant = TRUE,
  event_id_col = NULL,
  unknown_action = c("error", "treat_as_untested"),
  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins",
    "first")
)
```

## Arguments

- data:

  Data frame, one row per event x class (long) or one row per event with
  one column per class (wide) – auto-detected from whether `class_cols`
  resolve to one shared value column or to several distinct columns.

- group_col:

  Character. Column identifying the stratum (hospital, district, lab,
  ...). Not assumed to mean "hospital".

- class_cols:

  Character vector. Antibiotic classes in the panel. For long data,
  these are the values `class_col` takes. For wide data, these are
  column names, each holding an S/I/R-type result.

- class_col:

  Character. For long-format data, the column naming which class a row's
  result belongs to. Ignored for wide data (when `class_cols` resolve to
  real column names in `data`).

- value_col:

  Character. Column holding the susceptibility result. Resistant is
  matched case-insensitively against `resistant_values`; susceptible
  against `susceptible_values`; anything else (including `NA`) counts as
  untested for that class.

- resistant_values, susceptible_values, intermediate_values:

  Character vectors of values in `value_col` counted as resistant /
  susceptible / intermediate. Defaults cover the common
  \`"R"\`/\`"S"\`/\`"I"\` and
  \`"Resistant"\`/\`"Susceptible"\`/\`"Intermediate"\` spellings.
  Intermediate results still count toward `tested` regardless of
  `intermediate_as_resistant` – an ambiguous result is still a result,
  not a missing test.

- intermediate_as_resistant:

  Logical. When `TRUE` (default, matching the clinical convention
  already used by
  [`prep_recode_intermediate_ast()`](https://saketlab.github.io/anumaan/reference/prep_recode_intermediate_ast.md)
  elsewhere in `anumaan`), intermediate results are folded into
  `resistant` for the class-level counts and into the "R" side of the
  pairwise RR/RS/SR/SS breakdown. When `FALSE`, intermediate results are
  reported only in their own `intermediate` column/count and excluded
  from the binary resistant/susceptible tallies (class-level) and from
  the pairwise cell breakdown (they still count toward `co_tested`,
  since co-testing is about coverage, not the result).

- event_id_col:

  Character or `NULL`. If supplied, used to dedupe rows before counting
  complete-panel support (one row per event x class is assumed
  otherwise).

- unknown_action, conflict_action:

  Passed through to the shared AST policy in
  [`.profile_tested_long_wide`](https://saketlab.github.io/anumaan/reference/dot-profile_tested_long_wide.md):
  `unknown_action` (default `"error"`) controls what happens to a
  non-blank `value_col` value that matches none of `resistant_values`/
  `susceptible_values`/`intermediate_values`; `conflict_action` (default
  `"error"`) controls what happens when the same group x event x class
  cell has duplicate rows with DIFFERING results.

## Value

A list with three data frames:

- `class`:

  one row per group x class: `tested`, `resistant`, `susceptible`,
  `intermediate` (the raw intermediate count, always reported regardless
  of how it was folded into `resistant`) – kept for backward
  compatibility – plus the separated AST-policy counts `n_raw_tested`
  (every non-blank raw value seen, including unrecognized ones, before
  dedup), `n_binary_observed` (= `tested`, the resolved R/S/I count
  after dedup/conflict handling), `n_R`/`n_S`/ `n_I` (raw per-category
  counts, before `intermediate_as_resistant` folding), `n_untested`
  (genuinely blank/`NA`/declared untested-token rows), `n_invalid`
  (non-blank, unrecognized values – only nonzero when
  `unknown_action = "treat_as_untested"`; otherwise such a value stops
  the function), and `n_duplicates_collapsed` (identical-result
  duplicate rows collapsed to one)

- `pair`:

  one row per group x unordered class pair: `co_tested` (both classes
  have any result, R/I/S), and the `n_RR`/`n_RS`/`n_SR`/`n_SS` 2x2
  breakdown plus `near_separation` (`TRUE` when any one of those four
  cells is below 5, the same historical near-separation rule used for
  the J2-L hospital tiers)

- `profile`:

  one row per group: `complete_panel`, the count of events with every
  class in `class_cols` tested
