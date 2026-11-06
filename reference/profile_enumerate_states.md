# Enumerate the 2^D Joint Resistance States

Builds every combination of resistant (+1) / susceptible (-1) across `D`
antibiotic classes. Row order is the canonical state ordering used by
[`profile_interaction_design`](https://saketlab.github.io/anumaan/reference/profile_interaction_design.md)
and
[`profile_compatible_states`](https://saketlab.github.io/anumaan/reference/profile_compatible_states.md)
– do not reorder after generating it for a given model, the Stan data
has to agree with this ordering.

## Usage

``` r
profile_enumerate_states(D)
```

## Arguments

- D:

  Integer. Number of antibiotic classes in the panel.

## Value

Integer matrix, `2^D` rows x `D` columns, entries in `c(-1, 1)`
(susceptible / resistant).
