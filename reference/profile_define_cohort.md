# Build an Eligible-Group Roster From Precision-Evaluated Support

The function that was missing from the J2-L history: takes the pass/fail
table from
[`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)
and actually decides which groups (hospitals, ...) are in, producing
both the roster vector and a traceable manifest – replacing a hand-typed
list with a real, reproducible, data-derived decision.

## Usage

``` r
profile_define_cohort(
  thresholds,
  min_tier = "tier2",
  estimands = c("class", "pair", "profile"),
  mode = c("all", "complete_case")
)
```

## Arguments

- thresholds:

  Output of
  [`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)
  (the list form, with `$class`/`$pair`/`$profile`).

- min_tier:

  Character. Which `passes_<tier>` column must be `TRUE` – one of the
  names used in `tiers` when
  [`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md)
  was run (e.g. `"tier1"`, `"tier2"`, ...). Default `"tier2"` (the
  middle of the default three-tier set, `tier1=5/tier2=7.5/tier3=10`
  points).

- estimands:

  Character vector, any of `"class"`, `"pair"`, `"profile"`. A group
  must pass `min_tier` on every requested estimand, for every
  class/pair/the profile, to be included. Default all three (the
  strictest, most defensible choice – narrow this explicitly if e.g.
  only class margins matter for a given use).

- mode:

  Character, `"all"` or `"complete_case"`. Recorded in the manifest only
  – this function does not itself subset event-level data; it documents
  which cohort definition the roster is meant for.

## Value

A list:

- `roster`:

  character vector of included group names

- `manifest`:

  data frame, one row per group, with a logical `included` column and a
  `reason` string naming the first failing estimand/class when excluded
