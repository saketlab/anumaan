# Static (non-stochastic) per-event setup for prior predictive simulation – mirrors the STATIC portion of `.probit_predictive_draws_setup()` (posterior draw arrays are irrelevant here; a prior predictive draw needs only the real design structure – X_event, grouping indices, class panel – as conditioning predictors, per Part 10).

Static (non-stochastic) per-event setup for prior predictive simulation
– mirrors the STATIC portion of
[`.probit_predictive_draws_setup()`](https://saketlab.github.io/anumaan/reference/dot-probit_predictive_draws_setup.md)
(posterior draw arrays are irrelevant here; a prior predictive draw
needs only the real design structure – X_event, grouping indices, class
panel – as conditioning predictors, per Part 10).

## Usage

``` r
.prior_predictive_draws_setup(fitted_model, prior_config_override)
```
