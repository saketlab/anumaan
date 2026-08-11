# Run `expr` under a temporary, locally-seeded RNG state, restoring the caller's global RNG state (`.Random.seed`) afterwards – so a predictive-check generator never disturbs any OTHER stochastic code running in the same session before or after it is used.

Run `expr` under a temporary, locally-seeded RNG state, restoring the
caller's global RNG state (`.Random.seed`) afterwards – so a
predictive-check generator never disturbs any OTHER stochastic code
running in the same session before or after it is used.

## Usage

``` r
.ppc_with_local_seed(seed, expr)
```
