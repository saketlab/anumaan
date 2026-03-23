# Select Resistance Class

Selects a single resistance class per event using beta-lactam hierarchy
and relative risk (RR) values. Prevents double-counting in burden
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

  Data frame with resistance and RR information

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
  uses default from get_beta_lactam_hierarchy().

- filter_resistant:

  Logical. If TRUE, only consider resistant (R) classes. Default TRUE.

## Value

Data frame filtered to one resistance class per event

## Details

Selection logic: 1. Filter to resistant classes only (R) 2. Apply
beta-lactam hierarchy (Carbapenems \> 4GC \> 3GC \> ...) 3. Within same
hierarchy rank, prioritize by RR value (higher RR first) 4. If tied,
select alphabetically for reproducibility

## Examples

``` r
if (FALSE) { # \dontrun{
# Select single resistance class per event
selected <- select_resistance_class(data)

# Include susceptible classes too
selected <- select_resistance_class(data, filter_resistant = FALSE)
} # }
```
