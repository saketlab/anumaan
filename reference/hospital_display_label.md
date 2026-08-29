# Human-readable hospital/site display label

Replaces underscores with spaces only – does not alter the underlying
\`center_name\`/hospital identifier used anywhere else (joins, grouping,
file paths all keep the raw value; this is a plotting-display transform
only).

## Usage

``` r
hospital_display_label(x)
```

## Arguments

- x:

  Character vector of raw hospital/site names (e.g.
  `"AIIMS_trauma_center"`).

## Value

Character vector, e.g. `"AIIMS trauma center"`.
