# Filter Minimally Usable Records

Keeps records that have the bare minimum fields needed for any analysis:
a patient identifier, a culture date, an organism name, and at least one
non-NA AST result. This is a looser filter than
[`prep_filter_analysis_ready()`](https://saketlab.github.io/anumaan/reference/prep_filter_analysis_ready.md).

## Usage

``` r
prep_filter_minimally_usable(
  data,
  patient_col = "patient_id",
  culture_date_col = "culture_date",
  organism_col = "organism_name",
  ast_col = "ast_value_harmonized"
)
```

## Arguments

- data:

  Data frame after preprocessing.

- patient_col:

  Character. Patient ID column. Default "patient_id".

- culture_date_col:

  Character. Culture date column. Default "culture_date".

- organism_col:

  Character. Organism column. Default "organism_name".

- ast_col:

  Character. Harmonized AST column. Default "ast_value_harmonized".

## Value

Filtered data frame with only minimally usable rows.
