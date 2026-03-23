# Create Event IDs from Patient-Level Data

Converts patient-level isolate data to event-level data by grouping
repeat cultures from the same infection episode.

## Usage

``` r
create_event_ids(
  data,
  patient_col = "patient_id",
  date_col = "date_of_culture",
  organism_col = "organism_normalized",
  specimen_col = "specimen_type",
  antibiotic_col = "antibiotic_name",
  value_col = "antibiotic_value",
  gap_days = 14
)
```

## Arguments

- data:

  Data frame (long format – one row per antibiotic test).

- patient_col:

  Patient ID column. Default "patient_id".

- date_col:

  Culture date column. Default "date_of_culture".

- organism_col:

  Organism column. Default "organism_normalized".

- specimen_col:

  Specimen type column. Default "specimen_type".

- antibiotic_col:

  Antibiotic name column. Default "antibiotic_name".

- value_col:

  Susceptibility result column (S/I/R). Default "antibiotic_value".

- gap_days:

  Days threshold: gap \> gap_days triggers a new event. Default 14.

## Value

Original data frame with event_id column added.

## Details

Event classification rules (in order of precedence):

\| Scenario \| Result \| \|—————————————————\|—————-\| \| Different body
sites (any date) \| Separate events\| \| Same site, same day, same
organism, same ABG \| One event \| \| Same site, same day, same
organism, diff ABG \| Separate events\| \| Same site, same organism,
within gap_days, same ABG \| Same event \| \| Same site, same organism,
within gap_days, ABG changed \| New event \| \| Same site, same
organism, after \> gap_days \| New event \| \| Same site, same day,
different organisms \| Separate events\|

Event IDs are numbered GLOBALLY per patient in chronological order, so a
patient's events across all organisms/sites form one consistent
sequence.
