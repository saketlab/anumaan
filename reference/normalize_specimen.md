# Normalize Specimen/Sample Type

Normalizes specimen/sample type names and adds sample_category and
sterile_classification from the reference CSV file.

## Usage

``` r
normalize_specimen(data, specimen_col = "specimen_type", add_categories = TRUE)
```

## Arguments

- data:

  Data frame with specimen column

- specimen_col:

  Character. Specimen column name. Default "specimen_type".

- add_categories:

  Logical. Add sample_category and sterile_classification. Default TRUE.

## Value

Data frame with specimen_normalized, sample_category, and
sterile_classification columns

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(specimen = c("Blood", "Urine", "CSF"))
result <- normalize_specimen(data, specimen_col = "specimen")
} # }
```
