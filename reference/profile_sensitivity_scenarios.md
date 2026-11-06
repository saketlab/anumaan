# Named Sensitivity-Scenario Catalog

The 10 named scenarios a scientific interpretation of this family's
output should check, as a plain reference table – NOT something this
function runs. Running any of them means fitting the compared models
explicitly (one
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
call per side, under whatever changed assumption the scenario names) and
passing the results to
[`profile_sensitivity`](https://saketlab.github.io/anumaan/reference/profile_sensitivity.md);
this function never searches for a narrow interval or picks a model for
the caller, matching
[`profile_sensitivity`](https://saketlab.github.io/anumaan/reference/profile_sensitivity.md)'s
own standing rule.

## Usage

``` r
profile_sensitivity_scenarios()
```

## Value

A data frame, one row per scenario: `scenario_id`, `name`, `comparable`
(what must stay the same between the two sides for the contrast to mean
anything), `tests` (what the comparison actually tells you).
