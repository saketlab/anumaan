# (Internal) Log Loss

LogLoss = -1/n \* sum\[y\*log(q) + (1-y)\*log(1-q)\] – lower is better,
the opposite sign convention from this file's own per-row `log_score`
(kept for backward compatibility: higher is better there). Both are
reported; do not confuse the two.

## Usage

``` r
.log_loss(y, p)
```
