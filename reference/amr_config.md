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
  date_columns = c("date_of_admission", "date_of_culture", "date_of_final_outcome",
    "DOB"),
  hai_cutoff = 2,
  infer_department = TRUE,
  event_gap_days = 14,
  mortality_window = 14,
  age_bins = "GBD_standard",
  contaminant_method = "auto",
  mdr_definition = "CDC",
  xdr_definition = "CDC",
  map_icd10 = TRUE,
  rr_table = "GBD_2021",
  organism_map = "default",
  antibiotic_map = "WHO_2023",
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

  Character vector of date column names to parse. Default includes
  standard date fields.

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

  Character or numeric. MDR definition: "CDC", "WHO", or numeric
  threshold for number of resistant classes. Default "CDC".

- xdr_definition:

  Character or numeric. XDR definition. Default "CDC".

- map_icd10:

  Logical. Attempt ICD-10 code mapping. Default TRUE.

- rr_table:

  Character. Name of built-in RR table to use, or path to custom RR
  table. Default "GBD_2021".

- organism_map:

  Character or named vector. "default" for built-in mapping, or custom
  named vector. Default "default".

- antibiotic_map:

  Character. Built-in antibiotic classification: "WHO_2023" or path to
  custom. Default "WHO_2023".

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
