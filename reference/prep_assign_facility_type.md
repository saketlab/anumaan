# Assign a facility type / sector column from a user-supplied mapping

Adds a column that labels each row's facility with a category defined by
the caller (for example "Government" vs "Private", "Urban" vs "Rural",
or a region). The function is deliberately generic: it contains no
hardcoded facility names, so any dataset can be classified by passing an
appropriate `mapping`.

## Usage

``` r
prep_assign_facility_type(
  data,
  facility_col = "center_name",
  mapping,
  new_col = "facility_type",
  default = NA_character_,
  case_insensitive = TRUE
)
```

## Arguments

- data:

  A data frame.

- facility_col:

  Name of the column holding facility names. Default `"center_name"`.

- mapping:

  A named character vector where names are facility values and values
  are the type/sector to assign, e.g.
  `c(AFMC = "Government", SGRH = "Private")`.

- new_col:

  Name of the column to create. Default `"facility_type"`.

- default:

  Value assigned to facilities not found in `mapping`. Default
  `NA_character_`.

- case_insensitive:

  Logical. If `TRUE` (default), matching ignores case and surrounding
  whitespace, so `"AFMC"`, `"afmc"` and `" AFMC "` all match.

## Value

`data` with `new_col` added.

## Examples

``` r
df <- data.frame(center_name = c("AFMC", "SGRH", "Unknown"))
m <- c(AFMC = "Government", SGRH = "Private")
prep_assign_facility_type(df, mapping = m)
#>   center_name facility_type
#> 1        AFMC    Government
#> 2        SGRH       Private
#> 3     Unknown          <NA>
```
