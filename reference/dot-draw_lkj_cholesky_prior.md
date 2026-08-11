# Draw `n_states * n_matrices` iid LKJ-distributed correlation-matrix Cholesky factors via Stan's own `lkj_corr_cholesky_rng()`, in a single `fixed_param` sample() call (no NUTS, no warmup).

Draw `n_states * n_matrices` iid LKJ-distributed correlation-matrix
Cholesky factors via Stan's own `lkj_corr_cholesky_rng()`, in a single
`fixed_param` sample() call (no NUTS, no warmup).

## Usage

``` r
.draw_lkj_cholesky_prior(D, eta, n_matrices, n_states, seed)
```

## Value

`array(dim = c(n_states, n_matrices, D, D))`.
