# Prioritize Resistance

Helper function that applies hierarchy + RR ranking to select the most
important resistance class per event.

## Usage

``` r
prioritize_resistance(data, event_col, class_col, rr_col = NULL, hierarchy)
```

## Arguments

- data:

  Data frame

- event_col:

  Character. Event ID column.

- class_col:

  Character. Class column.

- rr_col:

  Character or NULL. RR column. If NULL, uses hierarchy only.

- hierarchy:

  Named numeric vector. Hierarchy mapping.

## Value

Data frame with one row per event (highest priority class)
