# Plot posterior predictive checks for a fitted probit model

Compact multi-page PDF: observed-vs-replicated marginal resistance,
number-resistant-per-event distribution, observed-vs-replicated pairwise
RR/RS/SR/SS, complete-profile concentration/entropy, hospital
heterogeneity, admission/patient clustering, and (when `ppc_replicates`
with materialised replicate arrays is supplied) a small-multiples page
comparing the observed data against `n_small_multiples` randomly
selected posterior replications. A good posterior predictive plot makes
it visually apparent whether the observed data look typical among the
replicated datasets.

## Usage

``` r
plot_probit_posterior_predictive_checks(
  ppc_statistics,
  ppc_replicates = NULL,
  output_pdf_path,
  title_base,
  n_small_multiples = 10L,
  small_multiple_seed = 123L
)
```

## Arguments

- ppc_statistics:

  Tibble returned by
  [`compute_probit_ppc_statistics`](https://saketlab.github.io/anumaan/reference/compute_probit_ppc_statistics.md).

- ppc_replicates:

  Optional `"amr_ppc_draws"` object returned by
  [`simulate_probit_posterior_predictive`](https://saketlab.github.io/anumaan/reference/simulate_probit_posterior_predictive.md)
  with `return_replicates = TRUE` – enables the small-multiples page.
  `NULL` (default) skips that page.

- output_pdf_path:

  Path to the PDF to write.

- title_base:

  Character; prefixed to every page's title (e.g.
  `"<experiment_id>\n<pathogen>"`).

- n_small_multiples:

  Integer; number of replicate states to show on the small-multiples
  page. Default 10.

- small_multiple_seed:

  Integer seed for selecting which replicate states appear on the
  small-multiples page. Default 123.

## Value

Invisibly returns the path to the saved PDF, or `NULL` if skipped.
