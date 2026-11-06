# Which Joint States Are Compatible With a Partially-Observed Event

For one event with possibly-missing classes, returns the indices (into
`states`, i.e. into `1:2^D`) of every joint state consistent with what
was actually tested. This is the exact-marginalization machinery:
untested classes are summed out by including every state that agrees
with the tested classes' signs, regardless of what the untested classes
would have been.

## Usage

``` r
profile_compatible_states(states, observed_sign)
```

## Arguments

- states:

  Output of
  [`profile_enumerate_states`](https://saketlab.github.io/anumaan/reference/profile_enumerate_states.md).

- observed_sign:

  Integer vector length `D`, entries in `c(-1, 0, 1)`: `-1`/`1` for
  susceptible/resistant, `0` (or `NA`) for untested.

## Value

Integer vector of row indices into `states`.
