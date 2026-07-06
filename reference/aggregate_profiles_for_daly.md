# Aggregate Posterior Profile Draws into R_ALL and R_NF Summaries

Summarises draw-level R_ALL and R_NF values from
[`compute_event_profile_probabilities()`](https://saketlab.github.io/anumaan/reference/compute_event_profile_probabilities.md)
into posterior mean and credible interval per hospital x pathogen x
profile combination. All \\2^D\\ profiles are present in the output
(inherited from the complete enumeration in the simulation step).

## Usage

``` r
aggregate_profiles_for_daly(
  profile_output,
  hospital_col = "hospital",
  pathogen_col = "pathogen",
  estimand = "observed_stewardship_event_mix",
  ci_level = 0.95
)
```

## Arguments

- profile_output:

  List returned by
  [`compute_event_profile_probabilities()`](https://saketlab.github.io/anumaan/reference/compute_event_profile_probabilities.md).

- hospital_col:

  Character. Hospital column in `aggregate_draws`. Must match the upper
  RE column used during fitting. Default `"hospital"`.

- pathogen_col:

  Character. Pathogen column. Default `"pathogen"`.

- estimand:

  Character. Estimand label to attach to output. Default
  `"observed_stewardship_event_mix"`.

- ci_level:

  Numeric. Credible interval coverage. Default `0.95`.

## Value

Tibble with one row per hospital x pathogen x profile: R_ALL and R_NF
posterior mean and credible interval.
