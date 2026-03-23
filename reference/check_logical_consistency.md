# Check Logical Consistency

Validates logical relationships in the data (e.g., date sequences, age
consistency, valid value ranges).

## Usage

``` r
check_logical_consistency(data, checks = "all", stop_on_failure = FALSE)
```

## Arguments

- data:

  Data frame

- checks:

  Character vector. Which checks to perform: - "date_sequence":
  admission \< culture \< outcome - "age_range": Age between 0-120 -
  "age_dob_match": Age matches DOB - "outcome_consistency": Died
  patients have outcome date - "all": All checks (default)

- stop_on_failure:

  Logical. Stop on inconsistency. Default FALSE.

## Value

List with consistency check results: - consistent: Logical -
issues_found: Data frame of inconsistent rows - summary: Character
vector of issue summaries

## Examples

``` r
if (FALSE) { # \dontrun{
consistency <- check_logical_consistency(data, checks = "all")

if (!consistency$consistent) {
  print(consistency$summary)
  View(consistency$issues_found)
}
} # }
```
