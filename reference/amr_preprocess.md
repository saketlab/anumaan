# AMR Data Preprocessing Pipeline

Master function that orchestrates the complete AMR data preprocessing
pipeline. Applies standardization, enrichment, and derivation in
sequence.

## Usage

``` r
amr_preprocess(
  data,
  config = NULL,
  phases = "all",
  verbose = TRUE,
  validate = TRUE,
  generate_report = TRUE
)
```

## Arguments

- data:

  Data frame. Raw AMR dataset.

- config:

  Configuration object from amr_config(). If NULL, uses defaults.

- phases:

  Character vector. Which phases to run: "standardize", "enrich",
  "derive", or "all" (default). Allows partial pipeline execution.

- verbose:

  Logical. Print detailed progress messages. Default TRUE.

- validate:

  Logical. Run validation checks before and after. Default TRUE.

- generate_report:

  Logical. Generate preprocessing report. Default TRUE.

## Value

List with class "amr_result": - data: Preprocessed data frame - config:
Configuration used - log: List of processing logs and summaries -
report: Preprocessing report (if generate_report = TRUE) - metadata:
Pipeline execution metadata

## Details

Pipeline phases: 1. \*\*Standardization\*\*: Column mapping, value
normalization, date parsing 2. \*\*Enrichment\*\*: Derive missing
optional variables (Age, LOS, infection type) 3. \*\*Derivation\*\*:
Create analytical variables (event IDs, MDR/XDR, weights)

## Examples

``` r
if (FALSE) { # \dontrun{
# Full pipeline with defaults
result <- amr_preprocess(raw_data)
clean_data <- result$data

# Custom configuration
config <- amr_config(
  hai_cutoff = 3,
  mdr_definition = "CDC",
  fuzzy_match = TRUE
)
result <- amr_preprocess(raw_data, config = config)

# Run only standardization phase
result <- amr_preprocess(raw_data, phases = "standardize")

# Summary
summary(result)
} # }
```
