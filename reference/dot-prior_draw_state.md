# Draw one prior predictive parameter state and reconstruct mu, mirroring the Stan fitting model's transformed-parameters block EXACTLY: `re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi]`, `mu = X_event beta + re_contribution(re_effect, flat_re_idx)`.

Draw one prior predictive parameter state and reconstruct mu, mirroring
the Stan fitting model's transformed-parameters block EXACTLY:
`re_effect[, lo:hi] = diag_pre_multiply(tau_re[r], L_corr_re[r]) * z_re[, lo:hi]`,
`mu = X_event beta + re_contribution(re_effect, flat_re_idx)`.

## Usage

``` r
.prior_draw_state(setup, lkj_re_arr, lkj_residual_arr, s)
```

## Arguments

- setup:

  Resolved prior-predictive setup.

- lkj_re_arr:

  LKJ Cholesky draws for random-effect correlation blocks.

- lkj_residual_arr:

  LKJ Cholesky draws for the residual correlation.

- s:

  Prior-state draw index.
