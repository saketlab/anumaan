# Interaction-Order Design Matrix Over the State Space

For interaction order `k`, returns a `2^D x choose(D, k)` matrix whose
column for class-combination `c` is the elementwise product of those `k`
classes' sign columns in `states` – the standard saturated log-linear /
Ising-model design. `k = 1` just returns `states` itself (the
main-effect design, `B` in the Stan model).

## Usage

``` r
profile_interaction_design(states, k)
```

## Arguments

- states:

  Output of
  [`profile_enumerate_states`](https://saketlab.github.io/anumaan/reference/profile_enumerate_states.md).

- k:

  Integer. Interaction order, `1 <= k <= ncol(states)`.

## Value

Numeric matrix, `nrow(states) x choose(ncol(states), k)`, with
`attr(, "combn")` recording which classes (by column index of `states`)
each output column corresponds to.
