# (Internal) Average Precision

Area under the recall-weighted precision step function, AP = sum_k
(Recall_k - Recall\_{k-1}) \* Precision_k. Depends on resistant
prevalence (unlike AUROC); undefined, not zero, with zero positive
cases. Matches the standard step-function definition (e.g.
scikit-learn's `average_precision_score`); ties at a threshold are not
specially collapsed – the telescoping sum is invariant to tie order, so
this is not an approximation.

## Usage

``` r
.average_precision(y, p)
```
