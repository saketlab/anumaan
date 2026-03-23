# Deduplicate Events

Removes duplicate isolate tests within the same event. Keeps the first
occurrence of each organism-antibiotic combination per event.

## Usage

``` r
deduplicate_events(
  data,
  event_col = "event_id",
  organism_col = "organism_normalized",
  antibiotic_col = "antibiotic_normalized",
  keep = "first"
)
```

## Arguments

- data:

  Data frame with event_id column

- event_col:

  Character. Event ID column. Default "event_id".

- organism_col:

  Character. Organism column. Default "organism_normalized".

- antibiotic_col:

  Character. Antibiotic column. Default "antibiotic_normalized".

- keep:

  Character. Which duplicate to keep: "first", "last", or "all". Default
  "first".

## Value

Deduplicated data frame
