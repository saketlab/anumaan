# Build Stan-Ready Design for a D-Class Joint Bernoulli Resistance Model

The general skeleton: takes a wide S/R outcome matrix and a
caller-supplied covariate design, and returns everything a Stan program
needs to fit the saturated (or order-capped) log-linear joint Bernoulli
model with exact marginalization of untested classes. `D` – and
therefore every interaction count – comes entirely from `ncol(Y)`; this
function does not hardcode a panel size or a covariate set. The
interaction coefficients (`eta`) have NO stratum index, so they are
shared across every level of whatever fixed-effect stratifier is encoded
into `X` (hospital, state, district, country, or none) – this is the
SAME fit regardless of how many strata `X` represents, never a separate
fit per stratum (see
[`profile_build_fixed_effect_design`](https://saketlab.github.io/anumaan/reference/profile_build_fixed_effect_design.md),
the recommended way to build `X` with that stratifier correctly
specified).

## Usage

``` r
profile_build_design(
  Y,
  X,
  class_names = NULL,
  max_order = NULL,
  beta_sd = NULL,
  interaction_sd = function(k) 1.5/(k - 1),
  grainsize = 64L
)
```

## Arguments

- Y:

  Numeric/integer matrix or data frame, `N` events x `D` classes,
  entries in `c(0, 1, NA)` (susceptible / resistant / untested). Column
  names become the class labels unless `class_names` is supplied.

- X:

  Numeric matrix, `N` x `K`: the covariate design for the main effects.
  Entirely up to the caller – build it with
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) from
  whatever formula is appropriate (age, sex, hospital, ICU, whatever).
  Must have `nrow(X) == nrow(Y)`.

- class_names:

  Character vector length `D`, optional (default: `colnames(Y)`).

- max_order:

  Integer, passed to
  [`profile_bernoulli_terms`](https://saketlab.github.io/anumaan/reference/profile_bernoulli_terms.md).
  Default: fully saturated (`= D`).

- beta_sd:

  Numeric scalar or length-`K` vector: prior SD for the main-effect
  coefficients. Default `2.5` for an intercept-like first column, `1.5`
  for the rest, matching J2-L's convention, but fully overridable.

- interaction_sd:

  Numeric scalar, or a function `order -> sd`, giving the prior SD for
  each interaction order's global scalars. Default halves the SD per
  additional order (`function(k) 1.5 / (k - 1)`), matching J2-L's 1.5 /
  0.75(\*) / 0.5 / 0.35 decreasing schedule in spirit ((\*) J2-L's
  order-2 term is actually a separate ICU/location effect, not part of
  this skeleton's covariate block – fold it into `X` instead).

- grainsize:

  Integer. Stan `reduce_sum` grainsize. Default 64L.

## Value

A named list, directly usable as `cmdstanr`/`rstan` Stan data, plus a
`manifest` data frame (one row per parameter, from
[`profile_bernoulli_terms`](https://saketlab.github.io/anumaan/reference/profile_bernoulli_terms.md))
and `design` (the state enumeration + per-order design matrices, for
inspection/plotting). Stan list elements: `N_events, D, K, S` (`= 2^D`),
`P` (total interaction-term count), `X`, `B` (`S x D` main-effect
design), `V` (`S x P` concatenated higher-order design), `n_compatible`,
`compatible_state` (`N_events x S`, only the first `n_compatible[e]`
entries of row `e` are used), `beta_sd` (length `K`), `eta_sd` (length
`P`).

## See also

The reference Stan model this pairs with:
`system.file("stan", "joint_bernoulli_general.stan", package = "anumaan")`
