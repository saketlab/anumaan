# Check Sampler Validity and Render Diagnostic Plots

The numeric half of "did this Stan run actually work": checks the
standard HMC/NUTS validity criteria (divergences, max-treedepth hits,
E-BFMI, split-Rhat, bulk/tail ESS) against widely used thresholds from
the Stan/`posterior` documentation, and renders the three standard
diagnostic plots (trace, rank, energy) via `bayesplot` – the same
`mcmc_trace()`/`mcmc_rank_overlay()` functions already used elsewhere in
`anumaan` (`plot_worst_offender_diagnostics.R`), plus
`mcmc_nuts_energy()` for the energy check that file doesn't currently
render.

## Usage

``` r
profile_check_sampler(
  fit_result,
  rhat_max = 1.01,
  ess_min = 400L,
  ebfmi_min = 0.2,
  variables = NULL,
  output_dir = NULL,
  ci_level = 0.95,
  derived_estimates = NULL,
  mcse_budget = NULL
)
```

## Arguments

- fit_result:

  The list returned by
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  (reads `$fit`, `$summary`, `$diagnostics`), or a raw `cmdstanr`
  `CmdStanMCMC` object (summary/diagnostics are then computed fresh from
  it).

- rhat_max:

  Numeric. Split-Rhat above this value fails a parameter. Default `1.01`
  (Vehtari et al. 2021's revised threshold).

- ess_min:

  Integer. Bulk or tail ESS below this fails a parameter. Default `400L`
  (the Stan-documented minimum for stable 95 interval estimates).

- ebfmi_min:

  Numeric. E-BFMI below this, per chain, fails that chain. Default `0.2`
  (Betancourt 2017's threshold).

- variables:

  Character vector or `NULL`. Which parameters' Rhat/ESS to check and
  plot. Default `NULL`: every parameter in the summary except `lp__`.

- output_dir:

  Character or `NULL`. If supplied, writes `sampler_check.rds`
  (status/reasons/checks, not the plot objects) and each plot as its own
  PDF into this directory.

- ci_level:

  Numeric in (0, 1). Credible-interval level the quantile-MCSE checks
  (both `mcse` and `derived_mcse`) are computed at. Default `0.95` (MCSE
  of the 0.025 and 0.975 quantiles –
  [`posterior::mcse_quantile()`](https://mc-stan.org/posterior/reference/mcse_quantile.html)'s
  `ess_tail()` alone summarizes a DIFFERENT combination of tails and
  does not replace checking the actual interval endpoints you report).

- derived_estimates:

  Named list or `NULL`. Each element a numeric matrix of posterior draws
  (`n_draws x n_quantities`) for a DERIVED estimand not among
  `variables` – e.g. per-hospital resistance probabilities from
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)'s
  `class_probability_draws`/`profile_probability_draws`, reshaped to one
  column per hospital x class/state. When supplied, mean and quantile
  MCSE are computed per column and returned as `derived_mcse` – raw
  sampler diagnostics on `beta`/ `eta` do not by themselves say whether
  a DERIVED probability is precise enough to report. Default `NULL`
  (skipped). NOTE: unlike `mcse` (computed on the fit's own chain-aware
  draws), `derived_estimates` columns are treated as one pooled sequence
  –
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  does not track which chain each subsampled draw came from, so no
  between-chain correction is applied here; this MCSE is a reasonable,
  slightly optimistic approximation, not the exact chain-aware quantity
  `mcse` reports.

- mcse_budget:

  Numeric or `NULL`. The MCSE budget each `derived_estimates` column's
  mean and quantile MCSE is checked against – e.g. `0.1 * epsilon`,
  where `epsilon` is the approved probability half-width target (a
  project convention, not something this function derives). `NULL`
  (default): MCSE is still reported, just not pass/failed against
  anything. Extending sampling addresses a numerically UNSTABLE estimate
  (MCSE too wide for the budget); it does NOT make a long but genuinely
  imprecise scientific interval (wide because of real data/model
  uncertainty) narrower – these are different problems, and this check
  only diagnoses the former.

## Value

A list:

- `status`:

  `"pass"`, `"warning"`, or `"fail"`. `"warning"` means only ESS
  (precision, not sampler validity) failed to clear `ess_min`; any
  divergence, treedepth, E-BFMI, or Rhat failure is `"fail"`.
  MCSE-budget failures (on `derived_estimates`) do NOT affect `status` –
  they are a separate, reporting-level judgment the caller makes
  explicitly, not a sampler-validity verdict.

- `reasons`:

  character vector of failed checks' names (empty when
  `status == "pass"`)

- `checks`:

  data frame, one row per check: `check`, `value`, `threshold`, `status`

- `mcse`:

  data frame, one row per `variables` entry: `variable`, `mcse_mean`,
  `mcse_lower_quantile`, `mcse_upper_quantile` (at `ci_level`'s
  endpoints)

- `derived_mcse`:

  empty (zero rows) unless `derived_estimates` is supplied; one row per
  (estimate name, column): `estimate_name`, `column`, `mcse_mean`,
  `mcse_lower_quantile`, `mcse_upper_quantile`, `budget`,
  `passes_budget` (`NA` when `mcse_budget` is `NULL`)

- `plots`:

  list of 3 `ggplot` objects: `trace`, `rank`, `energy`

## Details

This is deliberately upstream of any scientific judgment about the MODEL
(calibration, held-out performance, which classes/panels to trust) –
those are `profile_ppc`, `profile_validate_masked`,
`profile_validate_heldout`, and `profile_sensitivity` (elsewhere in this
family). A sampler can pass every check here and still fit a badly
miscalibrated model; this function only answers whether the SAMPLER did
its job, which is a precondition for trusting anything downstream of it.
