# Audit Direct Co-Tested Support for Each Interaction Order of a Candidate Panel

A candidate panel's model includes one interaction term per class subset
of order 2 up to `max_order` (pairs, triples, ...). This counts, for
each such subset and each group, how many events have EVERY class in
that subset tested together – i.e. whether that interaction term has any
direct co-tested support at all, as opposed to being estimated almost
entirely from the model's shared cross-class prior (`eta_sd` in
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)).

## Usage

``` r
profile_audit_interaction_support(
  data,
  group_col,
  panel,
  class_col = NULL,
  value_col = "antibiotic_value",
  resistant_values = c("R", "Resistant", "resistant"),
  susceptible_values = c("S", "Susceptible", "susceptible"),
  intermediate_values = c("I", "Intermediate", "intermediate"),
  event_id_col = NULL,
  max_order = NULL,
  min_co_tested = 5L,
  unknown_action = c("error", "treat_as_untested"),
  conflict_action = c("error", "drop_conflicts", "resistant_wins", "susceptible_wins",
    "first")
)
```

## Arguments

- data, group_col, class_col, value_col, resistant_values,
  susceptible_values, intermediate_values, event_id_col:

  Same as
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md).

- panel:

  Character vector. The candidate panel's classes.

- max_order:

  Integer or `NULL`. Highest interaction order to check. Default:
  `length(panel)` (every order up to the complete panel).

- min_co_tested:

  Integer. An operational warning threshold, not an identifiability
  criterion: an interaction term is flagged `has_direct_support = FALSE`
  for a group when fewer than this many events in that group have every
  class in the subset tested together. Default `5`, the same historical
  near-separation threshold used elsewhere in this family
  ([`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)'s
  `near_separation` flag).

- unknown_action, conflict_action:

  Same meaning as in
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)
  – the same shared AST policy in
  [`.profile_tested_long_wide`](https://saketlab.github.io/anumaan/reference/dot-profile_tested_long_wide.md),
  applied identically here.

## Value

A data frame, one row per (order, class subset, group): `order`,
`classes` (the subset, "; "-joined), `group`, `n_co_tested`,
`has_direct_support`.

## Details

**This is a structural-support audit, not an identifiability check, and
does not claim to be one.** `n_co_tested >= min_co_tested` establishes
neither statistical identification nor adequate information for a
subset's interaction term – it is an operational, configurable warning
threshold on raw co-testing counts, nothing more. A formal
identifiability assessment would examine a FITTED model's observed-data
information matrix (or the rank of the observed-probability map) for
unidentified or very weak directions, and whether posterior conclusions
are prior-dominated rather than data-driven – that requires an actual
fit, not just the raw co-testing pattern, and is a separate, later
concern (closer to `profile_sensitivity` once a fit exists, not yet
built into this family). This function is the cheap, pre-fit screen
only: subsets with zero or sparse direct co-tested support are flagged
before any time is spent fitting at all, nothing stronger is claimed.
