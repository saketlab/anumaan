# (Internal) Binned Calibration Curve

Mean predicted probability versus observed frequency, per bin of
predicted probability. Reports every non-empty bin, including its own
`n` – a thin bin's point is still shown, not silently smoothed over.

## Usage

``` r
.calibration_bins(y, p, n_bins = 10L)
```
