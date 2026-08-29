# Resistant-class-count summary statistics (Phase B redesign): mean, median, 90th percentile and maximum are counts of resistant classes and share one axis; variance is a dispersion statistic on a different scale (roughly the square of a count) and must not share that axis, even though the old .ppc_dotplot_page() call put all five on one row range labelled "Value". Split into two plots: Panel A (location/tail summaries) and Panel B (variance). There is no hospital/class dimension here (each statistic is a single dataset-wide row, stratum == "all_events"), so neither panel needs faceting.

Resistant-class-count summary statistics (Phase B redesign): mean,
median, 90th percentile and maximum are counts of resistant classes and
share one axis; variance is a dispersion statistic on a different scale
(roughly the square of a count) and must not share that axis, even
though the old .ppc_dotplot_page() call put all five on one row range
labelled "Value". Split into two plots: Panel A (location/tail
summaries) and Panel B (variance). There is no hospital/class dimension
here (each statistic is a single dataset-wide row, stratum ==
"all_events"), so neither panel needs faceting.

## Usage

``` r
.ppc_plot_resistant_count_summary(st, title_base)
```
