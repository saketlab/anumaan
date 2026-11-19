# (Internal) Joint-Profile Log Loss and Multiclass Brier, for Held-Out Events With a Fully Known True Profile

Cell-level scores (AUROC, Brier, average precision, calibration – see
[`.profile_score_mask_predictions`](https://saketlab.github.io/anumaan/reference/dot-profile_score_mask_predictions.md))
can agree perfectly between two models that disagree completely about
the joint profile: two models both giving P(classD=R)=P(classC=R)=0.5
can place all of their JOINT mass on entirely different complete states.
This scores the WHOLE profile as one unit instead, for exactly the
held-out events where the true complete profile is actually known (every
class genuinely tested) – partially-tested held-out events contribute no
row here, since there is no single true state to score against.

## Usage

``` r
.profile_joint_scores(est, Y, holdout_idx, ids)
```

## Details

`est$profile_probability` (from
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md),
`estimand = "unconditional"`) is already the posterior-MEAN probability
per event x state – i.e. already \\\bar q\_{e,s}\\, averaged over draws;
no further draw-level averaging happens here.
