# Shared posterior-draw setup for posterior predictive simulation

Mirrors the draw-subsampling / mu-reconstruction pattern of
`.probit_validation_draws_setup()` in `daly_resistance_validation.R`
(same
[`posterior::as_draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)
subsampling, same `.arr()` unpacking, same
[`re_contribution()`](https://saketlab.github.io/anumaan/reference/re_contribution.md)-based
mu reconstruction), but is intentionally a SEPARATE helper rather than a
shared/reused one: it needs the raw `L_Omega` Cholesky factor (not the
reconstructed/regularised `Omega` matrix the validation module needs for
`pbivnorm`), and `daly_resistance_validation.R` currently has zero test
coverage, so keeping predictive-check code isolated from it avoids
introducing any risk of silently changing existing validation behaviour.

## Usage

``` r
.probit_predictive_draws_setup(fitted_model, n_states, seed)
```
