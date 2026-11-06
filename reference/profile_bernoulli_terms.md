# Full Parameter Manifest for a D-Class Joint Bernoulli Model

Enumerates every term the model will fit: the `D` covariate-dependent
main effects, plus one global scalar per class-combination for every
interaction order from 2 up to `max_order` (default: fully saturated,
i.e. up to `D`). This is the thing that makes the model "automatic" in
the number of antibiotic classes chosen – adding or removing a class
changes `D`, which changes every count here without the caller
hand-editing interaction terms.

## Usage

``` r
profile_bernoulli_terms(D, class_names = NULL, max_order = D)
```

## Arguments

- D:

  Integer. Number of antibiotic classes.

- class_names:

  Character vector of length `D`, optional. Used to label interaction
  terms (e.g. `"Aminoglycosides:Carbapenems"`).

- max_order:

  Integer, `2 <= max_order <= D`. Cap the highest interaction order fit
  (e.g. `max_order = 2` for pairwise-only – useful when `D` is large
  enough that the full `2^D` state space or `choose(D, D)` terms become
  impractical). Default: fully saturated (`= D`).

## Value

A data frame, one row per interaction term, columns `order` (1 = main
effect), `term` (label), and `classes` (list-column of class-name
vectors). Row count for orders 2..max_order is
`sum_{k=2}^{max_order} choose(D, k)`.
