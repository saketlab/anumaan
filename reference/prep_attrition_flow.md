# Track Attrition Through the Pipeline

Records a stage entry in an attrition flow table. Call once per pipeline
stage to build up a cumulative record of how many rows and patients are
retained after each step.

## Usage

``` r
prep_attrition_flow(
  flow,
  data,
  stage_name,
  reason = "",
  patient_col = "patient_id",
  event_col = "event_id"
)
```

## Arguments

- flow:

  Data frame of previous attrition stages (or NULL to start).

- data:

  Current data frame at this stage.

- stage_name:

  Character. Short label for this stage.

- reason:

  Character. Description of what changed at this stage.

- patient_col:

  Character. Patient ID column. Default "patient_id".

- event_col:

  Character. Event ID column (optional). Default "event_id".

## Value

Updated attrition flow data frame.

## Details

Usage pattern:

      flow <- NULL
      flow <- prep_attrition_flow(flow, data_raw,    "raw_intake",    "All raw records")
      flow <- prep_attrition_flow(flow, data_dated,  "dates_parsed",  "After date coercion")
      flow <- prep_attrition_flow(flow, data_ready,  "analysis_ready","After all filters")
      print(flow)
