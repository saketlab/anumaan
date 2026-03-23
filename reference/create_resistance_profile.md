# Create Resistance Profile

Generates a resistance profile string for each event summarizing
resistance patterns. Useful for identifying common resistance
phenotypes.

## Usage

``` r
create_resistance_profile(
  data,
  event_col = "event_id",
  antibiotic_col = "antibiotic_normalized",
  susceptibility_col = "antibiotic_value",
  format = "resistant_list",
  class_col = "antibiotic_class"
)
```

## Arguments

- data:

  Data frame with resistance data

- event_col:

  Character. Event ID column. Default "event_id".

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_normalized".

- susceptibility_col:

  Character. Susceptibility column. Default "antibiotic_value".

- format:

  Character. Output format: - "resistant_list": List resistant drugs
  only (default) - "full_pattern": Full S/R pattern string -
  "class_summary": Resistant classes only

- class_col:

  Character. Class column (required for format = "class_summary").
  Default "antibiotic_class".

## Value

Data frame with resistance_profile column added

## Examples

``` r
if (FALSE) { # \dontrun{
# List resistant drugs
data_with_profile <- create_resistance_profile(data)

# Full S/R pattern
data_with_profile <- create_resistance_profile(
  data,
  format = "full_pattern"
)

# Resistant classes only
data_with_profile <- create_resistance_profile(
  data,
  format = "class_summary"
)
} # }
```
