# Draw one prior predictive parameter state and reconstruct mu, mirroring the Stan fitting model's transformed-parameters block EXACTLY: `re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi]`, `mu = X_event beta + re_contribution(re_effect, flat_re_idx)`.

Draw one prior predictive parameter state and reconstruct mu, mirroring
the Stan fitting model's transformed-parameters block EXACTLY:
`re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi]`,
`mu = X_event beta + re_contribution(re_effect, flat_re_idx)`.

## Usage

``` r
.prior_draw_state(setup, lkj_arr, s)
```
