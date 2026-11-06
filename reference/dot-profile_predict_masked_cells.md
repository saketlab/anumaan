# (Internal) Predict Masked Cells' Resistance Probability

For each masked cell, the event's other classes (genuinely untested, OR
masked for this same check) are set to sign `0`; the posterior
predictive probability of the masked class being resistant is
\\P(\text{masked}=R \mid \text{other observed}) = \sum\_{s \in
\text{compatible} \cap R} p_s \big/ \sum\_{s \in \text{compatible}}
p_s\\, using
[`profile_compatible_states`](https://saketlab.github.io/anumaan/reference/profile_compatible_states.md)
on the event's own sign vector with the masked class(es) opened up –
exact, not simulated. Averaged over posterior draws.

## Usage

``` r
.profile_predict_masked_cells(fit_result, design, Y, masked_df, n_draws, seed)
```
