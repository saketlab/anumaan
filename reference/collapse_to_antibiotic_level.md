# Collapse to Antibiotic Level (OPTIONAL - Run When YOU Decide)

\*\*IMPORTANT\*\*: This function removes duplicate tests. Only run this
when you have reviewed your data and decided to collapse duplicates.

## Usage

``` r
collapse_to_antibiotic_level(
  data,
  event_col = "event_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized",
  susceptibility_col = "antibiotic_value",
  aggregation_rule = "any_R"
)
```

## Arguments

- data:

  Data frame with susceptibility results

- event_col:

  Character. Event ID column. Default "event_id".

- organism_col:

  Character. Organism column. Default "organism_normalized".

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_normalized".

- susceptibility_col:

  Character. Susceptibility column (S/I/R). Default "antibiotic_value".

- aggregation_rule:

  Character. Rule for aggregation: "any_R" (default, any R -\> R),
  "most_resistant" (R \> I \> S), "most_common" (mode).

## Value

Aggregated data frame (one row per event-organism-antibiotic)

## Details

Aggregates multiple test results for the same organism-antibiotic
combination within an event. Uses "any R -\> R" logic where resistance
in any test results in resistant classification.

\*\*When to use\*\*: After you've cleaned and normalized data, if you
have multiple tests for the same patient-organism-antibiotic and want
one result per combination.

\*\*What it removes\*\*: Duplicate rows based on event_id + organism +
antibiotic

## Examples

``` r
if (FALSE) { # \dontrun{
# Check for duplicates first
data %>%
  group_by(event_id, organism_normalized, antibiotic_normalized) %>%
  filter(n() > 1)

# Then decide to collapse using any R -> R
collapsed <- collapse_to_antibiotic_level(data)

# Or use most common result
collapsed <- collapse_to_antibiotic_level(
  data,
  aggregation_rule = "most_common"
)
} # }
```
