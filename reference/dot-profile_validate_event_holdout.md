# (Internal) Refit on a Training Subset and Score Whole Held-Out Events

Shared machinery for
[`profile_validate_newpatient`](https://saketlab.github.io/anumaan/reference/profile_validate_newpatient.md)
and
[`profile_validate_temporal`](https://saketlab.github.io/anumaan/reference/profile_validate_temporal.md):
refits on `Y[train_idx,]`/ `X[train_idx,]` via
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md),
then predicts EVERY held-out event's per-class probabilities using ONLY
its covariates (the `"unconditional"` estimand via
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
– no AST information from that event is used, matching what is actually
available at deployment), and scores against whichever classes were
ACTUALLY tested in the real (unmasked) `Y` for those events via the SAME
scoring machinery
[`profile_validate_masked`](https://saketlab.github.io/anumaan/reference/profile_validate_masked.md)
uses. A naturally untested cell contributes no score.

## Usage

``` r
.profile_validate_event_holdout(
  Y,
  X,
  ids,
  output_dir,
  holdout_idx,
  sampler_config,
  n_draws,
  seed,
  validation_mode
)
```
