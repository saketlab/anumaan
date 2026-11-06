# Minimum Complete-Event Count for a Simultaneous Profile-Precision Bound

A distribution-free sample-size benchmark for estimating ALL `K = 2^D`
joint-profile probabilities of a `D`-class panel at once, each within
absolute error `epsilon`, from independent complete events. Derived from
Hoeffding's inequality applied to each of the `K` profile indicators
plus a union bound across profiles (optionally also across `n_strata`
strata, e.g. hospitals, via a further Bonferroni correction):
\$\$n\_{\min} = \lceil \log(2 K S / \alpha) / (2 \epsilon^2) \rceil\$\$
where \\S\\ = `n_strata` (or 1 if `NULL`).

## Usage

``` r
profile_precision_bound(D, epsilon, alpha = 0.05, n_strata = NULL)
```

## Arguments

- D:

  Integer. Number of antibiotic classes in the candidate panel.

- epsilon:

  Numeric in `(0,1)`. Maximum acceptable absolute error per profile
  probability – a scientific/clinical choice, not something this
  function derives.

- alpha:

  Numeric in `(0,1)`. Significance level (so `1 - alpha` is the
  per-stratum confidence). Default `0.05`.

- n_strata:

  Integer or `NULL`. When supplied, applies a Bonferroni correction
  across this many strata (e.g. hospitals) so the `1 - alpha` guarantee
  holds jointly across all of them, not just for each one individually.
  `NULL` (default) gives a per-stratum guarantee only – the usual case
  when each hospital is checked on its own against a shared benchmark.

## Value

A list: `D`, `K` (`= 2^D`), `epsilon`, `alpha`, `n_strata`, `n_min`.

## Details

This is a different quantity from the Wilson/Goodman interval widths in
[`profile_plan_precision`](https://saketlab.github.io/anumaan/reference/profile_plan_precision.md):
Wilson/Goodman answer "how wide is the interval around one estimated
proportion", while this answers "how many complete events guarantee
every one of the `K` profile estimates is simultaneously within
`epsilon` of truth, with probability `1 - alpha`". Neither subsumes the
other; both are reported.

The bound is deliberately conservative (distribution-free, no credit for
partial events, no correction for clustering or non-representative
complete-case selection) and assumes independent, representative
complete events – stated limitations, not footnotes to ignore. It
assesses direct data support only; a stratum below the bound may still
be modelable via partial events and shared cross-class structure, but
that requires a separate model-based assessment (`profile_sensitivity`),
not this function.
