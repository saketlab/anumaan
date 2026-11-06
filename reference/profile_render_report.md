# Fixed-Order Evidence Report, With a Release Gate

Assembles whatever this family's other functions have ALREADY produced
into the fixed 10-section report order (see this file's header), writes
every tabular piece supplied to `output_dir` as CSV/RDS plus an artifact
manifest (relative paths and md5 content hashes), and renders one
navigable HTML file. Computes nothing new except the release gate (see
`release_status` below) and the manifest – every number in the report
comes from an argument the caller already computed elsewhere in this
family.

## Usage

``` r
profile_render_report(
  cohort_summary = NULL,
  assumptions = NULL,
  prior_predictive = NULL,
  design = NULL,
  sampler_check = NULL,
  ppc = NULL,
  heldout = NULL,
  estimates = NULL,
  sensitivity = NULL,
  release_table = NULL,
  named_use_review = NULL,
  output_dir
)
```

## Arguments

- cohort_summary:

  Output of
  [`profile_summarize_cohort_eligibility`](https://saketlab.github.io/anumaan/reference/profile_summarize_cohort_eligibility.md),
  or `NULL`. Feeds sections 1-2.

- assumptions:

  Output of
  [`profile_check_assumptions`](https://saketlab.github.io/anumaan/reference/profile_check_assumptions.md),
  or `NULL`. Feeds section 3, and the release gate (any
  `status == "violated"` row blocks release).

- prior_predictive:

  Output of
  [`profile_prior_predictive`](https://saketlab.github.io/anumaan/reference/profile_prior_predictive.md),
  or `NULL`. Feeds section 3.

- design:

  Output of
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md),
  or `NULL`. Feeds section 4 (`$manifest`, the term dictionary).

- sampler_check:

  Output of
  [`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md),
  or `NULL`. Feeds section 5, and the release gate (a
  `status == "fail"`, or this argument being entirely absent, blocks
  release – computational status unknown is never treated as
  computational status fine).

- ppc:

  Output of
  [`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md),
  or `NULL`. Feeds section 6. Any family not supplied (e.g.
  `count_groups` never computed) is reported as "Not available", never
  silently absent from the report's structure.

- heldout:

  Output of
  [`profile_validate_masked`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md),
  [`profile_validate_heldout`](https://saketlab.github.io/anumaan/reference/profile_validate_heldout.md),
  [`profile_validate_newpatient`](https://saketlab.github.io/anumaan/reference/profile_validate_newpatient.md),
  or
  [`profile_validate_temporal`](https://saketlab.github.io/anumaan/reference/profile_validate_temporal.md),
  or `NULL`. Feeds section 7.

- estimates:

  Output of
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md),
  or `NULL`. Feeds section 8.

- sensitivity:

  Output of
  [`profile_sensitivity`](https://saketlab.github.io/anumaan/reference/profile_sensitivity.md),
  or `NULL`. Feeds section 9.

- release_table:

  Output of
  [`profile_release`](https://saketlab.github.io/anumaan/reference/profile_release.md),
  or `NULL`. Feeds sections 9-10 (the per-group decision table).

- named_use_review:

  List or `NULL`. The human reviewer's recorded decision:
  `list(approved = TRUE/FALSE, named_use = "<text>", reviewer = "<text>", date = "<text>", notes = "<text>")`.
  `NULL` (default) means no review has been recorded yet –
  release_status is then at best `"exploratory"`, never `"approved"`, no
  matter how clean the rest of the evidence is.

- output_dir:

  Character. Directory for the HTML report, every supplied table's
  CSV/RDS, and `manifest.csv`.

## Value

A list: `release_status` (`"blocked"`, `"exploratory"`, or
`"approved"`), `release_reasons` (character vector, empty only when
`"approved"`), `report_path` (the written HTML file), `manifest` (data
frame: `file`, `relative_path`, `md5`).
