# Compute each event's total random-effect contribution from a flattened re_effect\[D, total_re_levels\] matrix (as emitted by the generic Stan models' generated re_effect, or reconstructed R-side from posterior draws of z_re/tau_re/L_corr_re). This is the SINGLE generic helper every downstream mu-reconstruction site should call instead of hand-summing hospital_effect/patient_effect/admission_effect.

Compute each event's total random-effect contribution from a flattened
re_effect\[D, total_re_levels\] matrix (as emitted by the generic Stan
models' generated re_effect, or reconstructed R-side from posterior
draws of z_re/tau_re/L_corr_re). This is the SINGLE generic helper every
downstream mu-reconstruction site should call instead of hand-summing
hospital_effect/patient_effect/admission_effect.

## Usage

``` r
re_contribution(re_effect, flat_group_index)
```

## Arguments

- re_effect:

  D x total_re_levels numeric matrix (one posterior draw).

- flat_group_index:

  n_events x R integer matrix (see
  prepare_random_effects()\$flat_group_index), or a single event's
  R-length integer vector.

## Value

If flat_group_index is a matrix: n_events x D matrix of summed RE
contributions. If a vector: length-D vector for one event.
