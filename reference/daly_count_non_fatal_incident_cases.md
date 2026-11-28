# Count Non-Fatal Incident Cases by Syndrome and Hospital

Counts the non-fatal incident cases of an infectious syndrome L at each
hospital h, \\I\_{h,L}\\: the number of distinct **discharged
admissions** with syndrome L. This is the direct count – no CFR, no
`CR_L` adjustment and no pathogen weighting. It is the incidence used in
the YLD equations (see
[`daly_calc_yld_associated()`](https://saketlab.github.io/anumaan/reference/daly_calc_yld_associated.md)).

## Usage

``` r
daly_count_non_fatal_incident_cases(
  data,
  syndrome_col,
  syndrome_name = NULL,
  admission_col,
  outcome_col,
  discharged_value = "Discharged",
  facility_col = NULL,
  facility_name = NULL
)
```

## Arguments

- data:

  Data frame with one or more rows per admission.

- syndrome_col:

  Character. Column holding the syndrome of each admission (e.g.
  `"chosen_syndrome"`).

- syndrome_name:

  Character vector or `NULL`. Syndrome(s) to count (matched ignoring
  case and spacing). `NULL` counts every syndrome present.

- admission_col:

  Character. Admission identifier column (e.g. `"readmission_id"`).
  Cases are distinct admissions, not rows or patients.

- outcome_col:

  Character. Final outcome column.

- discharged_value:

  Character. Value of `outcome_col` marking a discharged (non-fatal)
  admission. Default `"Discharged"`.

- facility_col:

  Character or `NULL`. Hospital column. When given, counts are returned
  per hospital; `NULL` pools all records.

- facility_name:

  Character or `NULL`. Restrict to one hospital.

## Value

Data frame with `facility_col` (when given), `syndrome_col` and `N_NF_L`
(distinct discharged admissions). Hospital x syndrome combinations with
no discharged admission are not listed.

## Details

Each admission is counted once per syndrome, however many isolates,
pathogens or antibiotic rows it has. Admissions with a missing or blank
syndrome are not counted under any syndrome; a message reports how many
discharged admissions that excludes. A warning is raised if an admission
carries more than one syndrome, because it is then counted under each.

## Examples

``` r
adm <- data.frame(
  readmission_id = c("A1", "A1", "A2", "A3", "A4", "A5"),
  hospital = c("H1", "H1", "H1", "H1", "H2", "H2"),
  syndrome = c("BSI", "BSI", "BSI", "UTI", "BSI", ""),
  outcome = c("Discharged", "Discharged", "Death", "Discharged", "Discharged", "Discharged")
)
daly_count_non_fatal_incident_cases(adm, "syndrome", NULL, "readmission_id", "outcome",
  facility_col = "hospital"
)
#> 1 discharged admission(s) have no syndrome and are not counted under any syndrome.
#> Non-fatal incident cases: 3 discharged admission(s) across 2 syndrome(s) and 2 hospital(s).
#>   hospital syndrome N_NF_L
#> 1       H1      BSI      1
#> 2       H1      UTI      1
#> 3       H2      BSI      1
```
