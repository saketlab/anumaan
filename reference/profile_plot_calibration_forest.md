# Observed-vs-Model Calibration Forest Plot

One CI forest plot per
[`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md)
level: the model's posterior mean + credible interval (horizontal error
bar, blue) beside the observed value (a red x at the same row) – a
reader sees at a glance whether the empirical point falls inside the
model's interval, exactly what `interval_contains_observed` encodes
numerically. `status == "not_computable"` rows are dropped, not drawn as
empty/misleading rows. Faceted by the entity (class, class pair, or
state) so many hospitals compare at a glance within one panel.

## Usage

``` r
profile_plot_calibration_forest(
  ppc_table,
  level = c("marginal", "pairwise", "complete_profile"),
  title_base = "",
  max_facets = 12L
)
```

## Arguments

- ppc_table:

  One of
  [`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md)'s
  three tables (`$marginal`, `$pairwise`, or `$complete_profile`).

- level:

  Character, one of `"marginal"`, `"pairwise"`, `"complete_profile"` –
  which table `ppc_table` is, since the observed/model/entity column
  names differ by level.

- title_base:

  Character. Prefixed to the plot title.

- max_facets:

  Integer or `NULL`. For `"complete_profile"` (up to `2^D` facets,
  unreadable fast for `D > 4`), keeps only the `max_facets` states with
  the largest mean model probability across groups. `NULL`: keep all.
  Ignored for `"marginal"`/`"pairwise"`. Default `12L`.

## Value

A ggplot object, or `NULL` if nothing in `ppc_table` was evaluated.
