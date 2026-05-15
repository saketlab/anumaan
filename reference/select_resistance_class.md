# Select Resistance Class for Burden Attribution

Selects a single resistance class per event using beta-lactam hierarchy
and relative risk (RR) values. Prevents double-counting in DALY burden
estimation by choosing the most clinically relevant resistant class.

## Usage

``` r
select_resistance_class(
  data,
  event_col = "event_id",
  class_col = "antibiotic_class",
  susceptibility_col = "antibiotic_value",
  rr_col = "rr_value",
  hierarchy = NULL,
  filter_resistant = TRUE
)
```

## Arguments

- data:

  Data frame with class-level resistance and RR columns.

- event_col:

  Character. Event ID column. Default "event_id".

- class_col:

  Character. Antibiotic class column. Default "antibiotic_class".

- susceptibility_col:

  Character. Susceptibility column. Default "antibiotic_value".

- rr_col:

  Character. RR value column. Default "rr_value". If missing, only
  hierarchy is used.

- hierarchy:

  Named numeric vector. Custom hierarchy (class name -\> rank). If NULL,
  uses `get_beta_lactam_hierarchy()`.

- filter_resistant:

  Logical. If TRUE, only consider resistant (R) classes. Default TRUE.

## Value

Data frame filtered to one resistance class per event.

## Details

Selection order: 1. Beta-lactam hierarchy rank (Carbapenems \> 4GC \>
3GC \> ...) 2. RR value (higher relative risk = higher priority) 3.
Alphabetical tie-breaker
