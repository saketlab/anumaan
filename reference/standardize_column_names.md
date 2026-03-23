# Standardize Column Names to Package Convention

Maps incoming dataset column names to standardized names used throughout
the package. Supports exact matching and optional fuzzy matching for
unmatched columns.

## Usage

``` r
standardize_column_names(
  data,
  mapping = default_column_mappings,
  fuzzy_match = TRUE,
  fuzzy_threshold = 0.3,
  interactive = FALSE
)
```

## Arguments

- data:

  A data frame with raw column names

- mapping:

  Named list where names are standard column names and values are
  character vectors of acceptable aliases. Default uses
  `default_column_mappings`.

- fuzzy_match:

  Logical. If TRUE, attempts fuzzy matching for unmapped columns using
  string distance. Default TRUE.

- fuzzy_threshold:

  Numeric. Maximum string distance (0-1) for fuzzy matching. Lower
  values require closer matches. Default 0.3.

- interactive:

  Logical. If TRUE and fuzzy matches found, prompts user for
  confirmation. Default FALSE (auto-accept).

## Value

A list with components:

- data: Data frame with standardized column names

- mapping_log: List documenting which columns were mapped and how

- unmapped: Character vector of columns that couldn't be mapped

## Examples

``` r
if (FALSE) { # \dontrun{
raw_data <- data.frame(
  PatientID = 1:10,
  Organism = rep("E. coli", 10),
  Drug = rep("Ampicillin", 10)
)
result <- standardize_column_names(raw_data)
clean_data <- result$data
} # }
```
