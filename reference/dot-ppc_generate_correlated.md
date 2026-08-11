# Correlated-residual posterior predictive draw: Z_rep ~ MVN(mu, Omega), Y_rep = I(Z_rep \> 0)

Reuses the state's already-stored Cholesky factor `L_omega_s` directly
(one Cholesky per state, per Part 20's performance requirement) via a
single D x N_events matrix multiply, vectorised across every event in
that state – no per-event loop.

## Usage

``` r
.ppc_generate_correlated(mu_s, L_omega_s)
```
