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
  \>=2 observations falls below this is a hard stop (via
  on_mostly_singleton = "stop") rather than a warning. NULL (default)
  never stops here – eligibility gating belongs in the analysis layer
  (see anumaan-analysis's random_effect_eligibility.R), not hardcoded
  into the package.

- on_mostly_singleton:

  "warn" (default) or "stop" when a block is more than 90 percent
  singleton levels (1 observation).

## Value

An object of class "amr_random_effects": block metadata, level maps,
per-event flattened group indices, nesting diagnostics.
