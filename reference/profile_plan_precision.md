# Evaluate Precision Against a Set of Named Tiers

Takes the raw counts from
[`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)
and computes, for each row, how wide a confidence interval a reported
proportion from that count would actually have – then marks whether that
width clears each requested precision tier. Class-level rows use the
Wilson score interval (single proportion); pairwise and profile rows use
the conservative simultaneous Goodman interval (several proportions
reported from one group at once), matching the same two interval
families used in `methodology/j2_support_precision/` in
anumaan-analysis.

## Usage

``` r
profile_plan_precision(
  support,
  tiers = c(tier1 = 5, tier2 = 7.5, tier3 = 10),
  conf_level = 0.95,
  n_categories = NULL
)
```

## Arguments

- support:

  One element of the list returned by
  [`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md)
  (i.e. `$class`, `$pair`, or `$profile`), or the whole list (all three
  are evaluated).

- tiers:

  Named numeric vector, half-width in percentage points for each named
  tier, ordered from most to least demanding. Default
  `c(tier1 = 5, tier2 = 7.5, tier3 = 10)` – the same three
  percentage-point bars used historically for J2-L, just named plainly
  by rank (`tier1` = narrowest interval required = hardest to clear)
  instead of labels like "strict"/"moderate"/"loose" that imply an
  interpretation this function doesn't make. Any number of tiers is fine
  – add a `tier4`, rename them, whatever a given use needs.

- conf_level:

  Numeric. Confidence level for the interval. Default 0.95.

- n_categories:

  Integer. Number of simultaneous categories for the Goodman interval on
  pair/profile rows (e.g. 4 for the RR/RS/SR/SS cells of one pair, or
  `2^D` for a full D-class profile). Required when evaluating
  `$pair`/`$profile` directly; ignored for `$class`.

## Value

The same data frame(s) as the input, with an added `interval_half_width`
column (percentage points, at the worst-case/conservative proportion for
Goodman rows) and one logical `passes_<tier>` column per entry in
`tiers` (e.g. `passes_tier1`, `passes_tier2`, ...).
