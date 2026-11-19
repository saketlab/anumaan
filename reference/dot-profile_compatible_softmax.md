# (Internal) Numerically Stable Softmax Over a Compatible-State Subset

Normalizes raw state scores within the compatible subset `idx`,
stabilized against that subset's own maximum. Normalizing over all
states first and then restricting can underflow every compatible entry
to 0 when those states score far below the global maximum. Shared by
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
and the masked-cell predictor so both paths use the same computation.

## Usage

``` r
.profile_compatible_softmax(score, idx)
```

## Arguments

- score:

  Numeric vector of length S (unnormalized log-scores).

- idx:

  Integer indices of compatible states (non-empty).

## Value

Numeric vector, same length as `idx`, summing to 1.
