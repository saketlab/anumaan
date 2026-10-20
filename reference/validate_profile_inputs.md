# Validate Inputs for Resistance Profile Estimation (Pathway 1)

Checks that all mandatory columns exist, detects long vs wide format,
validates AST content, and verifies that columns required for any
requested stratification are present. Returns a tibble of check results
and stops on any failed mandatory check.

## Usage

``` r
validate_profile_inputs(
  data,
  col_map = list(isolate_col = "isolate_id", pathogen_col = "pathogen", ast_col =
    "ast_value", patient_col = "patient_id", date_col = "date_of_culture", geography_col
    = "state", specimen_col = "specimen_type", age_col = "age", dob_col = "dob",
    antibiotic_col = "antibiotic_name", class_col = "antibiotic_class", location_col =
    NULL, outcome_col = NULL),
  stratify_by = NULL,
  outcome_col = NULL
)
```

## Arguments

- data:

  Data frame to validate.

- col_map:

  Named list mapping logical roles to actual column names. Defaults use
  standard package column names. Set any entry to `NULL` to skip that
  check (only valid for optional roles).

- stratify_by:

  Character vector or `NULL`. Requested stratification dimensions.
  Accepted values: `"geography"`, `"year"`, `"outcome"`. Triggers
  additional column checks for each dimension. Default `NULL`.

- outcome_col:

  Character or `NULL`. Outcome column name. Required when `"outcome"` is
  in `stratify_by`. Can also be supplied via `col_map$outcome_col`.
  Default `NULL`.

## Value

A tibble with columns `check_name`, `status` (`"pass"` / `"warn"` /
`"fail"`), `message`, and `n_affected`. The attribute `detected_format`
(`"long"` or `"wide"`) is attached. Stops with an informative message if
any mandatory check fails.

## Details

**Mandatory columns (always):**

- `isolate_col` – unique isolate identifier; the counting unit for all
  resistance statistics. One patient may have multiple isolates.

- `pathogen_col` – organism / pathogen name.

- `ast_col` – susceptibility result; must contain S / I / R.

- `patient_col` – patient identifier (retained for context, not
  counted).

- `date_col` – culture date; used to derive `year` when `"year"` is in
  `stratify_by`.

- `geography_col` – geographic unit (state / district).

- `specimen_col` – specimen type (blood, urine, etc.).

- `age_col` OR `dob_col` – at least one must be present.

**Mandatory in long format only:**

- `antibiotic_col` – antibiotic / drug name column.

**Optional (warn if absent, do not stop):**

- `class_col` – antibiotic class; derived from `antibiotic_col` via WHO
  AWaRe reference if absent.

- `location_col` – ward / ICU / OPD; retained for subgroup outputs.

- `outcome_col` – patient outcome; required for outcome stratification.

## Examples

``` r
ast_long <- data.frame(
  isolate_id = rep(sprintf("ISO%02d", 1:6), each = 2),
  patient_id = rep(sprintf("P%02d", 1:6), each = 2),
  pathogen = rep(c("Klebsiella pneumoniae", "Escherichia coli"), each = 6),
  antibiotic_name = rep(c("Meropenem", "Ceftriaxone"), 6),
  ast_value = rep(c("R", "S", "S", "R"), 3),
  date_of_culture = as.Date("2024-01-01") + rep(0:5, each = 2),
  specimen_type = "Blood",
  state = "Maharashtra",
  age = 40
)
report <- validate_profile_inputs(ast_long)
#> 
#> [validate_profile_inputs] 11 passed | 1 warnings | 0 failed
#>   [!] Antibiotic class column 'antibiotic_class' absent. Will derive classes from antibiotic names via WHO AWaRe reference.
names(report)
#> [1] "check_name" "status"     "message"    "n_affected"
```
