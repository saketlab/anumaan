# Deterministic per-state seed, a function of (seed, s) only

Ensures a given state's replicate depends ONLY on (seed, s, theta_s) –
never on call order, on interleaving with other generator objects built
from the same or a different seed, or on any other RNG consumption
elsewhere in the calling session between generator construction and the
first `generate_state(s)` call. Small differences in the seed passed to
[`set.seed()`](https://rdrr.io/r/base/Random.html) produce uncorrelated
Mersenne Twister streams, so a simple affine combination is sufficient
here (this is not a cryptographic hash, just a reproducible
index-to-seed map).

## Usage

``` r
.ppc_state_seed(seed, s)
```
