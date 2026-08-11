# Build the generic random-effect representation for an arbitrary number of blocks

Flattens an arbitrary number of declared random-intercept blocks into
the representation used by the generic Stan models and by every
downstream mu-reconstruction site (profile generation, DALY aggregation,
validation).

## Usage

``` r
prepare_random_effects(
  data,
  random_effects,
  min_repeated_levels = NULL,
  singleton_threshold = 0.9,
  on_mostly_singleton = c("warn", "stop")
)
```

## Arguments

- data:

  Data frame containing every declared block's group_col.

- random_effects:

  Character vector (legacy) or list-of-blocks (see
  .normalize_random_effects_spec()).

- min_repeated_levels:

  Optional integer. If supplied, a block whose number of levels with
  \>=2 observations falls below this is checked INDEPENDENTLY of
  singleton_threshold below (previously this was bugged: nested inside
  the singleton_fraction check, so it was silently skipped whenever
  singleton_fraction fell at or below the threshold even if
  n_repeated_levels was itself too low). Triggers a stop (via
  on_mostly_singleton = "stop") or warning. NULL (default) never checks
  this – eligibility gating belongs in the analysis layer (see
  anumaan-analysis's random_effect_eligibility.R), not hardcoded into
  the package.

- singleton_threshold:

  Numeric in (0, 1\]. A block whose fraction of singleton levels
  (exactly 1 observation) exceeds this triggers a separate warning/stop,
  evaluated independently of min_repeated_levels. Default 0.90.

- on_mostly_singleton:

  "warn" (default) or "stop" – applies to BOTH the min_repeated_levels
  check and the singleton_threshold check.

## Value

An object of class "amr_random_effects": block metadata, level maps,
per-event flattened group indices, nesting diagnostics.
