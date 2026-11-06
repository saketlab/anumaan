# (Internal) Classify and Pivot Long/Wide AST Data for the profile\_\* Family

Shared by every `profile_*` function that needs a group x event x class
S/R/I matrix:
[`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
[`profile_audit_interaction_support`](https://saketlab.github.io/anumaan/reference/profile_audit_interaction_support.md),
and
[`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md)
all call this once instead of each re-implementing the same long-to-wide
classify-and-pivot. Not exported – this is specifically the
classification step (mapping raw values to canonical R/S/I via
`resistant_values`/etc.), which is resistance-profile-domain logic, not
a generic preprocessing pivot; it is deliberately NOT the same job as
`anumaan`'s
[`prep_create_wide_ast_matrix()`](https://saketlab.github.io/anumaan/reference/prep_create_wide_ast_matrix.md),
which pivots raw values verbatim with no classification and is for
preprocessing-stage use.

## Usage

``` r
.profile_tested_long_wide(
  data,
  group_col,
  class_cols,
  class_col,
  value_col,
  resistant_values,
  susceptible_values,
  intermediate_values,
  event_id_col,
  untested_values = c("NT", "Not tested", "not tested", "NA", "N/A"),
  unknown_action = c("error", "treat_as_untested"),
  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins",
    "first"),
  intermediate_as_resistant = NULL
)
```

## Details

One AST policy, applied identically everywhere this helper is used. Raw
value parsing tries the caller's declared `resistant_values`/
`susceptible_values`/`intermediate_values` vocabulary first (exact
match, case-insensitive), then falls back to
[`prep_clean_ast_values`](https://saketlab.github.io/anumaan/reference/prep_clean_ast_values.md)'s
messy-string extraction (handles `"S (HIGH LEVEL)"`, keyword scans,
etc.) for anything that didn't match – the same parsing logic used
everywhere else in `anumaan`, not a separate, cruder exact-match-only
implementation kept just here. blank/`NA`/declared `untested_values`
tokens are UNTESTED, not an outcome; a non-blank value that matches none
of the above (declared vocabulary OR `prep_clean_ast_values`'s
extraction) is UNKNOWN, handled per `unknown_action` (`"error"`, the
default – a typo or an unanticipated lab code should be loud, not
silently swept into "untested"; `"treat_as_untested"` to opt out when
that is a deliberate, declared choice). Duplicate rows for the SAME
group x event x class are resolved by delegating to
[`prep_deduplicate_ast`](https://saketlab.github.io/anumaan/reference/prep_deduplicate_ast.md)
(via `conflict_action`, passed through as its `strategy`) – not a
separate, parallel implementation: identical duplicates always collapse
to one first (an `n_duplicates_collapsed` count is returned, never
silently inflating support), and duplicates with DIFFERING R/S/I results
are a genuine data conflict, never resolved by taking the first row
unless that policy is explicitly chosen – `"error"` (the default) stops
with every conflicting cell named; `"drop_conflicts"` removes every row
for a conflicting cell (never guesses which one is right);
`"resistant_wins"`/ `"susceptible_wins"`/`"first"` resolve by the same
clinical priority or row-order policy `prep_deduplicate_ast` already
offers elsewhere in the package. Grouping is always by
`(.group, .event)` jointly, so an `event_id_col` value that repeats
ACROSS groups (e.g. a hospital-local patient id reused at a different
hospital) is never merged across hospitals – group-scoped identifiers
are the default, not an opt-in.

When `intermediate_as_resistant` is supplied (not `NULL`), an additional
`.result_folded` column (and its own wide pivot,
`wide_folded`/`result_col_folded`) is built by delegating to
[`prep_recode_intermediate_ast`](https://saketlab.github.io/anumaan/reference/prep_recode_intermediate_ast.md)
– the SAME Colistin-aware clinical convention already used elsewhere in
`anumaan` (Colistin `I -> S`, every other class `I -> R`), rather than a
separate, cruder symmetric `I -> R` rule. `FALSE` leaves
`.result_folded` identical to `.result` (`I` stays its own category, as
before). `NULL` (the default) skips this step entirely – only the
callers that fold at all
([`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
[`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md))
pass a non-`NULL` value.
