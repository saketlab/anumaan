# Build a Covariate Design Matrix With a Config-Chosen Fixed-Effect Stratifier

Builds the `X` matrix
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
expects, for the ONE joint model this family fits across every
configured stratum at once – never a silent separate fit per stratum.
The interaction coefficients
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
builds (`eta`) have no stratum index at all (they are defined purely
over the state space), so they are automatically SHARED across every
level of whatever stratifier this function encodes here; nothing
downstream needs to know what that stratifier was.

## Usage

``` r
profile_build_fixed_effect_design(
  ids,
  group_col = "group",
  covariates = NULL,
  contrast = c("reference", "sum_to_zero"),
  reference_level = NULL
)
```

## Arguments

- ids:

  Data frame, one row per event, matching
  [`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md)'s
  `$ids` row order (its own `$ids` is the intended input, but any data
  frame with a stratum column works).

- group_col:

  Character or `NULL`. Column in `ids` naming the fixed-effect
  stratifier – hospital, state, district, country, whatever the config
  specifies. Default `"group"`, matching
  `profile_prepare_model_data()$ids`'s own column name. `NULL` means NO
  fixed effects (fully pooled) – an explicit, valid choice, not
  something to avoid by default.

- covariates:

  Data frame or `NULL`, same row order as `ids`
  (`profile_prepare_model_data()$covariates` is the intended input).
  Additional covariates (age, sex, ...) included alongside the
  stratifier. Default `NULL` (stratifier + intercept only).

- contrast:

  Character, `"reference"` (default) or `"sum_to_zero"`. `"reference"`
  drops one level (`reference_level`) and keeps an intercept – standard
  treatment coding. `"sum_to_zero"` uses deviation-from-grand-mean
  coding instead, with no dropped reference level.

- reference_level:

  Character or `NULL`. Which level of `group_col` is dropped, when
  `contrast = "reference"`. `NULL` (default) uses the first level
  alphabetically – R's own default, made explicit and recorded here
  rather than left implicit.

## Value

A list: `X` (the design matrix), `term_names` (`colnames(X)`, in order),
`group_col` (echoed, `NULL` if no stratifier), `reference_level` (the
dropped level, `NA` for `sum_to_zero` or no stratifier), `contrast`
(echoed), `levels` (every level of `group_col` present, in the order
used), `rank` (`qr(X)$rank` – always `== ncol(X)`, or the function has
already stopped).

## Details

The stratifier itself is never assumed to be "hospital": it is whatever
column the caller's analysis config names via `group_col` – hospital,
state, district, country, or `NULL` for no fixed effects at all (a fully
pooled single-intercept model). That choice belongs to whoever prepares
the config in `anumaan-analysis`, exactly like `group_col` is already
generic everywhere else in this family
([`profile_audit_support`](https://saketlab.github.io/anumaan/reference/profile_audit_support.md),
[`profile_prepare_model_data`](https://saketlab.github.io/anumaan/reference/profile_prepare_model_data.md),
...) – this function does not special-case hospitals, and does not
decide the granularity for the caller.

Guards against the failure mode an unconstrained design invites: a full
set of stratum dummy columns PLUS an intercept is rank-deficient
(aliased) – `contrast = "reference"` drops one declared/default level
instead, and `contrast = "sum_to_zero"` uses deviation-from- grand-mean
coding with no dropped reference level at all. Either way, the returned
design's rank is checked and a rank deficiency (from this or from
collinear `covariates`) stops with a specific diagnosis rather than
silently fitting an unidentified model.
