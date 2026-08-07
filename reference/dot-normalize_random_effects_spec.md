# Normalize a random_effects specification into the canonical block-list form used throughout the package.

Accepts EITHER: - the legacy character vector, e.g. c("center_name",
"readmission_id") (temporary backward compatibility – each column
becomes its own block, using the COLUMN NAME itself as the block name,
never a numbered label); or - the new list-of-blocks form, e.g.
list(list(name = "hospital", group_col = "center_name", terms =
"intercept"), list(name = "admission", group_col = "admission_key",
terms = "intercept"))

## Usage

``` r
.normalize_random_effects_spec(random_effects)
```
