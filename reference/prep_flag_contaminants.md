# Flag Contaminant Organisms

Identifies likely contaminant organisms using multi-path logic (provided
flags, device-based context, or heuristic organism+specimen matching).

## Usage

``` r
prep_flag_contaminants(
  data,
  method = "auto",
  context = c("general", "icu_device"),
  organism_col = "organism_normalized",
  specimen_col = "specimen_type"
)
```

## Arguments

- data:

  Data frame.

- method:

  Character. "auto" (try all methods), "device_based", "heuristic", or
  "provided". Default "auto".

- context:

  Character. "general" or "icu_device". Default "general".

- organism_col:

  Character. Normalized organism column. Default "organism_normalized".

- specimen_col:

  Character. Specimen type column. Default "specimen_type".

## Value

Data frame with `is_contaminant`, `contaminant_confidence`, and
`contaminant_method` columns.

## Details

When `context = "icu_device"` (AIIMS ICU BSI data), CoNS in blood with a
central line present AND line in place \> 2 days are classified as
CLABSI pathogens rather than contaminants.
