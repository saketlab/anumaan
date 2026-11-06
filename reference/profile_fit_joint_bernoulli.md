# Fit the Joint Bernoulli Resistance-Profile Model

Compiles `inst/stan/joint_bernoulli_general.stan` (or a caller- supplied
alternative) and samples it against a Stan-ready design from
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md),
saving every artifact a later review needs: the fit object (`fit.rds`),
a parameter summary (`parameter_summary.csv`), per-chain sampler
diagnostics (`sampler_diagnostics.rds`), a snapshot of the exact sampler
config used, hashed for reproducibility (`config_snapshot.rds`), and a
short human-readable run log (`run_log.txt`). cmdstanr's own per-chain
CSV files are also written into `output_dir` (via `output_dir` forwarded
to `$sample()`).

## Usage

``` r
profile_fit_joint_bernoulli(
  design,
  output_dir,
  sampler_config = list(),
  stan_file = NULL,
  save_warmup = FALSE
)
```

## Arguments

- design:

  Output of
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (reads `design$stan_data`; `design$manifest` is recorded in the config
  snapshot for provenance).

- output_dir:

  Character. Directory for all artifacts (created if missing).

- sampler_config:

  Named list. Any of `chains` (default `4L`), `iter_warmup` (`1000L`),
  `iter_sampling` (`1000L`), `seed` (`123L`), `parallel_chains`
  (default: `= chains`), `adapt_delta` (`0.9`), `max_treedepth` (`10L`),
  `threads_per_chain` (`1L` – `> 1` enables `reduce_sum` threading in
  the Stan model's `grainsize` machinery). Unset entries use these
  defaults.

- stan_file:

  Character or `NULL`. Default: the package's
  `inst/stan/joint_bernoulli_general.stan`.

- save_warmup:

  Logical, forwarded to `$sample()`. Default `FALSE`.

## Value

A list: `fit` (the `cmdstanr` `CmdStanMCMC` object), `summary` (data
frame, `fit$summary()`), `diagnostics` (list from
`fit$diagnostic_summary()`), `config` (the resolved sampler config, the
design manifest, and `config_hash`), `output_dir`.
