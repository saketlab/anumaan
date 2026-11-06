# Run-State Monitor for a profile_fit_joint_bernoulli() Run

Reads whatever cmdstan has ALREADY WRITTEN to `output_dir` – never
re-runs or re-samples anything. Because
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)
calls cmdstanr's `$sample()` synchronously, this function cannot observe
a run live from inside the SAME call that started it; point it at the
same `output_dir` from a SEPARATE session/process while the fit is
running elsewhere for a genuine (coarse) mid-run progress read.

## Usage

``` r
profile_monitor_run(output_dir, expected_chains = NULL)
```

## Arguments

- output_dir:

  Character. The SAME `output_dir` passed to
  [`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md).

- expected_chains:

  Integer or `NULL`. Expected chain count – when supplied, a chain whose
  CSV does not exist at all is still reported as its own row
  (`phase = "not_started"`), rather than silently absent from the table.

## Value

A data frame, one row per found (or, with `expected_chains`, expected)
chain: `chain_id`, `csv_file` (`NA` if that chain's CSV does not exist),
`phase` (`"not_started"`, `"warmup"`, `"sampling"`, `"completed"`),
`last_iteration` (count of valid data rows written so far – NEVER counts
a truncated trailing line), `num_warmup`, `num_samples` (as configured,
from the CSV's own comment header), `last_modified` (file mtime – the
closest honest proxy to "last update time" this interface exposes),
`n_divergent`, `max_observed_treedepth` (both `NA` before any row is
written).

## Details

`phase` is inferred, not read directly – cmdstan's CSV records
`save_warmup`, so when that is `FALSE` (cmdstanr's default, and
[`profile_fit_joint_bernoulli`](https://saketlab.github.io/anumaan/reference/profile_fit_joint_bernoulli.md)'s
default), a chain with ZERO data rows written yet is reported `"warmup"`
(warmup draws are never written to the CSV in that mode, so no rows is
the ONLY signal available that warmup is still in progress – it is
indistinguishable from "about to start" with this interface, and this
function says so rather than claiming otherwise).
