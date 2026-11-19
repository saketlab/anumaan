# Audit the Observed Frequency of Every Resistance Profile, Per Group

Complements
[`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
aggregate `complete_panel` count with the full breakdown behind it: for
a `D`-class panel, how many of a group's complete-panel events fall into
each of the `2^D` specific resistance profiles (e.g. `"RRRR"`, `"RRSR"`,
..., `"SSSS"`), plus that profile's own empirical prevalence and Wilson
confidence interval.

## Usage

``` r
profile_audit_profile_distribution(
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
  conf_level = 0.95,
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

- conf_level:

  Numeric in (0, 1). Confidence level for the per-profile Wilson
  interval. Default 0.95.

- unknown_action, conflict_action:

  Passed through to the shared AST policy in
  [`.profile_tested_long_wide`](https://saketlab.github.io/anumaan/reference/dot-profile_tested_long_wide.md):
  `unknown_action` (default `"error"`) controls what happens to a
  non-blank `value_col` value that matches none of `resistant_values`/
  `susceptible_values`/`intermediate_values`; `conflict_action` (default
  `"error"`) controls what happens when the same group x event x class
  cell has duplicate rows with DIFFERING results.

## Value

A data frame with one row per group x profile (every group crossed with
every one of the `2^D` profiles enumerated by
[`enumerate_binary_profiles`](https://saketlab.github.io/anumaan/reference/enumerate_binary_profiles.md),
in the same class order as `class_cols`): `group`, `profile` (e.g.
`"RRSR"`), `n_hr` (observed complete-event count for that exact
profile), `n_h` (the group's total complete-panel count – identical to
`profile_audit_support(...)$profile$complete_panel`, repeated per
profile row for convenience), `prevalence` (`n_hr / n_h`, `NA` when
`n_h = 0`), `ci_lower`, `ci_upper` (the Wilson score interval around
`prevalence` at `conf_level`).

## Details

This is a pre-fit, data-only tabulation – it counts observed complete
events, nothing else. It is deliberately a separate function from
[`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)
(whose profile-level check is one conservative worst-case-\\p=0.5\\
Goodman bound over all `2^D` categories at once, from the total count
alone, never the per-profile split) and from
[`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md)
(which compares a FITTED model's posterior-replicated profiles against
these same observed counts, and therefore requires a completed fit).
Three distinct questions, three distinct outputs – this function answers
"how is the data actually distributed", not "is the total sample large
enough in the worst case" or "does the fitted model reproduce it".

Every one of the `2^D` profiles is reported for every group, zero-filled
when unobserved – a zero count is not treated as evidence of a zero
probability; the Wilson interval at `n_hr = 0` still returns a valid
(nondegenerate) upper bound. No arbitrary per-profile minimum count is
imposed; precision is read directly off the reported interval width,
which can differ enormously between the most and least common profile at
the same hospital even when the aggregate `complete_panel` count is
large.
