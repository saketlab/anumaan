# Remove Duplicate Rows

Identifies and removes exact duplicate rows from the dataset. Optionally
keeps first or last occurrence.

## Usage

``` r
remove_duplicate_rows(data, keep = "first", subset = NULL, report = TRUE)
```

## Arguments

- data:

  Data frame

- keep:

  Character. Which duplicate to keep: "first" (default), "last", or
  "none".

- subset:

  Character vector. Column names to check for duplicates. If NULL,
  checks all columns. Default NULL.

- report:

  Logical. If TRUE, prints detailed duplicate report. Default TRUE.

## Value

Data frame with duplicates removed

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove exact duplicates (all columns)
clean_data <- remove_duplicate_rows(data)

# Remove duplicates based on specific columns
clean_data <- remove_duplicate_rows(
  data,
  subset = c("patient_id", "date_of_culture", "organism_normalized")
)
} # }
```
