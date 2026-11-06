# SBC + Recovery/Robustness Simulation Harness

SBC + Recovery/Robustness Simulation Harness

## Usage

``` r
profile_run_simulation_study(
  design,
  n_sim,
  mode = c("sbc", "recovery"),
  true_params = NULL,
  missingness_fn = NULL,
  sampler_config = list(),
  output_dir,
  ci_level = 0.95,
  seed = 123L
)
```

## Arguments

- design:

  A template design from
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (its `X`, `B`, `V`, `beta_sd`, `eta_sd`, `D`, `K`, `P`, `S` are the
  structural skeleton every simulated dataset reuses – only `Y` is
  re-simulated).

- n_sim:

  Integer. Number of simulated datasets to fit. Each is a SEPARATE, real
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  run – genuinely expensive; size it to your compute budget and to the
  Monte Carlo precision you actually need (not a universal default).

- mode:

  Character, `"sbc"` or `"recovery"`.

- true_params:

  (mode `"recovery"` only) list with `beta` (`K x D` matrix) and `eta`
  (length-`P` vector) – FIXED true parameters every simulated dataset
  uses. Required for `"recovery"`; must be `NULL` for `"sbc"` (which
  draws its own true parameters from the declared prior, fresh per
  dataset).

- missingness_fn:

  Function or `NULL`. `function(Y)` -\> `Y` with some entries set to
  `NA`, applied to each simulated COMPLETE dataset before fitting – for
  declaring a sparse- mask or selective-testing scenario explicitly.
  Default `NULL` (every simulated dataset is fully observed).

- sampler_config, output_dir:

  Passed to every simulated fit.

- ci_level:

  Numeric. Default `0.95`. Used for `mode = "recovery"`'s coverage
  check.

- seed:

  Integer.

## Value

A list:

- `mode`:

  echoed

- `results`:

  one row per simulated dataset, INCLUDING failed ones (never silently
  dropped from the denominator): `sim_id`, `ok`, `failure_reason`,
  `sampler_status` (from
  [`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md),
  `NA` if the fit itself failed). For `mode = "sbc"`: one `rank_<param>`
  column per `beta`/`eta` parameter (an integer in `0..n_draws_used`,
  `NA` if not computable) and `n_draws_used` (after the ESS-based
  thinning mitigation for correlated draws). For `mode = "recovery"`:
  one `bias_<param>`, `covered_<param>`, `width_<param>` triple per
  parameter.

- `summary`:

  `mode = "sbc"`: a uniformity check per parameter (chi-squared
  goodness-of-fit of the rank histogram, a descriptive flag, NOT a
  formal hypothesis test with a universal threshold).
  `mode = "recovery"`: mean bias, empirical coverage, mean width, and
  failure RATE (over `n_sim`, not over converged fits only) per
  parameter.
