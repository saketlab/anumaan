# (Internal) Build Validated profile_build_design() Arguments for a Refit

Shared by every validation refit (new-patient, temporal, masked-cell).
Only `class_names`, `max_order`, `beta_sd` and `interaction_sd` are
accepted; any other name in `design_spec` is an error rather than being
silently ignored (a misspelled field would otherwise fall back to the
default unnoticed). A key is forwarded only when supplied: passing
`interaction_sd = NULL` would override
[`profile_build_design()`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)'s
own non-NULL default with NULL. An empty `design_spec` is allowed but
announced, since it means the refit uses package defaults, which may
differ from the model being validated.

## Usage

``` r
.profile_design_args(Y, X, design_spec, caller)
```
