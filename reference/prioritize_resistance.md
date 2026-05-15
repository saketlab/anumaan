# Prioritize Resistance Class

Internal helper for
[`select_resistance_class()`](https://saketlab.github.io/anumaan/reference/select_resistance_class.md).
Applies beta-lactam hierarchy + RR ranking to select one row per event.

## Usage

``` r
prioritize_resistance(data, event_col, class_col, rr_col = NULL, hierarchy)
```

## Arguments

- data:

  Data frame.

- event_col:

  Character. Event ID column.

- class_col:

  Character. Class column.

- rr_col:

  Character or NULL. RR column; if NULL uses hierarchy only.

- hierarchy:

  Named numeric vector mapping class name to rank.

## Value

Data frame with one row per event (highest priority class selected).
