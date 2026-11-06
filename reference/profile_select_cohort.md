# Select and Record the Modeling Cohort From Raw Event Data

Stage 0 of the resistance-profiling pipeline. Every eligibility,
precision, and panel-feasibility calculation downstream of this function
must be computed on the SAME cohort that will eventually be fit – never
on the full dataset while fitting only a restricted subset (e.g.
discharged-only events). This function makes that restriction an
explicit, traceable step instead of an assumption buried in a script.

## Usage

``` r
profile_select_cohort(data, filters, event_id_col = NULL)
```

## Arguments

- data:

  Data frame of events/isolates.

- filters:

  Named list. Each name is a column in `data`; each value is a
  character/numeric vector of values that column must match (via `%in%`)
  for a row to survive. Applied in list order (e.g.
  `list(final_outcome = "Discharged", center_name = my_hospitals)`).

- event_id_col:

  Character or `NULL`. Reserved for future event-level manifest detail;
  currently unused beyond validating it is present if supplied.

## Value

A list:

- `data`:

  the filtered data frame

- `manifest`:

  one row per filter applied, in order: `column`, `values`
  (comma-joined), `n_before`, `n_after`, `n_excluded`

## Details

Filters are applied in the order supplied, each an exact-match inclusion
rule against one column. The returned manifest shows how many rows
survived EACH filter in sequence, so a caller can see exactly how much a
given restriction removed on its own, rather than only one before/after
total that hides which criterion did the work.

This filters rows; it does not attach a clinical label to the result.
Restricting to `final_outcome == "Discharged"` describes exactly that
cohort, not automatically "all non-fatal infections" – transfers,
ongoing admissions, or other surviving-but-not-discharged outcomes are
excluded by that specific restriction, and the caller should label the
resulting cohort accordingly rather than assuming the broader label.
