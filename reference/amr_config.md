# Create AMR Preprocessing Configuration

Creates a configuration object that controls the behavior of the
preprocessing pipeline. Users can override default settings for their
specific use case.

## Usage

``` r
amr_config(
  column_mappings = NULL,
  fuzzy_match = TRUE,
  strict_validation = FALSE,
  date_columns = NULL,
  hai_cutoff = 2,
  infer_department = TRUE,
  event_gap_days = 14,
  mortality_window = 14,
  age_bins = "GBD_standard",
  contaminant_method = "auto",
  mdr_definition = "CDC",
  intermediate_as_resistant = TRUE,
  verbose = TRUE
)
```

## Arguments

- column_mappings:

  Named list of column name mappings. Default uses
  `default_column_mappings`.

- fuzzy_match:

  Logical. Enable fuzzy matching for column names. Default TRUE.

- strict_validation:

  Logical. If TRUE, stops execution if required fields are missing. If
  FALSE, issues warnings. Default FALSE.

- date_columns:

  Character vector of date column names to parse. Default NULL
  auto-detects date-like columns by name pattern.

- hai_cutoff:

  Numeric. Number of days after admission to classify as
  Hospital-Acquired Infection (HAI). Default 2.

- infer_department:

  Logical. Attempt to infer hospital department from other variables if
  missing. Default TRUE.

- event_gap_days:

  Numeric. Minimum days between events to create new event_id. Default
  14.

- mortality_window:

  Numeric. Days after culture to classify death as infection-related.
  Default 14.

- age_bins:

  Character. Age binning strategy: "GBD_standard", "pediatric", or
  "geriatric". Can also be a custom vector. Default "GBD_standard".

- contaminant_method:

  Character. Method for contaminant classification: "auto" (cascade
  through available methods), "device_based", "heuristic", "provided".
  Default "auto".

- mdr_definition:

  Character or numeric. MDR/XDR definition (used for both): "CDC",
  "WHO", or numeric threshold for number of resistant classes. Default
  "CDC".

- intermediate_as_resistant:

  Logical. Treat Intermediate (I) as Resistant (R) except for special
  cases like Colistin. Default TRUE.

- verbose:

  Logical. Print progress messages during processing. Default TRUE.

## Value

An `amr_config` object (list with class)

## Examples

``` r
if (FALSE) { # \dontrun{
# Use defaults
config <- amr_config()

# Customize for specific hospital
config <- amr_config(
  hai_cutoff = 3,
  mdr_definition = 5,
  strict_validation = FALSE
)
} # }
```
