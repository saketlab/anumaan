# Sensitivity Contrast Across Alternative Model Fits

Compares
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
output across two or more already-completed fits of the SAME design
(same panel, same covariates, same `Y`/`ids`) made under different
modelling assumptions – different interaction-prior scales (`eta_sd` via
[`profile_build_design()`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)'s
`interaction_sd`), different `max_order`, different missingness
handling, etc.

## Usage

``` r
profile_sensitivity(
  fit_results,
  design,
  Y,
  ids,
  reference = NULL,
  estimand = "observed_preserving_conditional",
  level = c("marginal", "complete_profile"),
  n_draws = 500L,
  seed = 123L
)
```

## Arguments

- fit_results:

  Named list, length `>= 2`, each element a
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
  return value (or raw `CmdStanMCMC`) fitted against the SAME
  `design`/`Y`/ `ids` under a different assumption. Names become this
  comparison's columns (e.g.
  `list(default_priors = fit1, independence_priors = fit2)`).

- design, Y, ids:

  Same as
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  – the SAME values for every fit (only the posterior draws differ
  between fits; `B`/`V`/`X` are structural and do not depend on the
  prior that was used when fitting).

- reference:

  Character or `NULL`. Name of the `fit_results` entry everything else
  is contrasted against. Default `NULL`: the first name in
  `fit_results`.

- estimand:

  Passed to
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  for every fit. Default `"observed_preserving_conditional"`.

- level:

  Character, `"marginal"` or `"complete_profile"`. `"marginal"`
  collapses
  [`profile_estimate()`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)'s
  per-state probabilities to per-class resistance probability first
  (easier to read across many classes); `"complete_profile"` contrasts
  the full per-state probabilities directly.

- n_draws, seed:

  Passed to
  [`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
  for every fit.

## Value

A data frame, one row per group x class (`level = "marginal"`) or group
x state (`level = "complete_profile"`): each fit's value (one column per
name in `fit_results`), one `<name>_minus_<reference>` SIGNED difference
column per non-reference fit (direction matters – a model that moves
every group's estimate up under the alternative assumption is a
different finding from one that moves some up and some down by the same
average magnitude, which `sensitivity_contrast` alone cannot
distinguish), and `sensitivity_contrast` – the maximum ABSOLUTE
difference from the reference across the other fits. This single number
answers "how much does this group's estimate move" but not "in which
direction" – use the signed columns for that.

## Details

This function does NOT run any fits itself, by design, matching this
package's standing rule that it never silently chooses a model for the
caller: the caller fits each alternative explicitly (one
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
call per candidate assumption, each its own tracked run with its own
sampler diagnostics to check via
[`profile_check_sampler`](https://saketlab.github.io/anumaan/reference/profile_check_sampler.md)
before trusting it in this comparison) and passes the resulting
`fit_results` in.
