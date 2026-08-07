# Conditional profile probabilities under a correlated residual, via Gibbs

For one hospital-pathogen panel at one posterior draw, samples the
missing latent class dimensions conditional on the observed resistance
SIGNS of the tested classes (we only ever observe
\\\mathrm{sign}(Z\_{ed})\\, never the exact latent value), via a Gibbs
sampler on the truncated multivariate normal \\Z_e \sim N(\mu_e,
\Omega)\\, vectorized across all events in the panel simultaneously (the
per-dimension conditional-normal regression coefficients depend only on
\\\Omega\\ – fixed across events at a given draw – so one matrix
operation updates every event's Gibbs state at once).

## Usage

``` r
.gibbs_conditional_profile_probs(
  mu_hp,
  Omega_hp,
  obs_panel,
  profile_bin,
  n_burnin = 10L,
  n_kept = 20L
)
```

## Arguments

- mu_hp:

  Numeric matrix, n_hp x Dp. Linear predictor for this draw.

- Omega_hp:

  Numeric matrix, Dp x Dp. Residual correlation matrix for this draw
  (already symmetrized/regularised by the caller).

- obs_panel:

  Numeric matrix, n_hp x Dp, values in {0, 1, NA}.

- profile_bin:

  Integer matrix, n_profiles x Dp, 0/1 (from
  [`enumerate_binary_profiles()`](https://saketlab.github.io/anumaan/reference/enumerate_binary_profiles.md));
  row \\j\\ must correspond to integer code \\j - 1\\ in binary (column
  \\d\\ = bit \\d - 1\\), which is how
  [`enumerate_binary_profiles()`](https://saketlab.github.io/anumaan/reference/enumerate_binary_profiles.md)
  orders rows.

- n_burnin, n_kept:

  Integer. Gibbs iterations discarded as burn-in and iterations kept for
  the empirical profile-probability count, respectively, at this single
  posterior draw. Total Gibbs cost scales with
  `(n_burnin + n_kept) * n_posterior_draws_for_profiles`, so
  correlated-residual profile generation should typically use fewer
  posterior draws than an identity-residual fit of the same size.

## Value

Numeric matrix, n_hp x n_profiles: empirical conditional profile
probability per event, rows summing to 1.

## Details

At every iteration, dimensions with an observed value are drawn from a
truncated normal restricted to the half-line matching the observed sign
(inverse-CDF sampling); dimensions without an observed value are drawn
from the unrestricted conditional normal. Because the observed
dimensions are truncation-constrained at every iteration, every kept
iteration's pattern is automatically consistent with the observed cells
– profiles inconsistent with observed data are never visited, rather
than assigned zero probability after the fact. Counting each kept
iteration's full binary pattern against the enumerated profile list
gives an empirical conditional profile-probability vector per event that
sums to exactly 1 by construction (every iteration falls into exactly
one profile).

`Dp == 1` (a single-class panel) has no joint structure to speak of and
is handled directly via \\\Phi(\mu\_{ed})\\, bypassing Gibbs entirely.
