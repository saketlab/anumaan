# Identity-residual posterior predictive draw: Y_rep ~ Bernoulli(Phi(mu))

Distributionally identical to drawing Z ~ N(mu, 1) and thresholding at
0, but avoids materialising the unneeded N_events x D latent Z matrix
and mirrors what
[`validate_marginal_calibration()`](https://saketlab.github.io/anumaan/reference/validate_marginal_calibration.md)
already does with `stats::pnorm(mu)` for the identity case.

## Usage

``` r
.ppc_generate_identity(mu_s)
```
