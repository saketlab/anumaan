# (Internal) Calibration Intercept/Slope

From a logistic regression of truth on logit(predicted probability) –
slope 1/intercept 0 is perfect calibration. NA when unestimable
(one-class truth, separation, glm non-convergence) rather than a
fabricated number.

## Usage

``` r
.calibration_slope_intercept(y, p)
```
