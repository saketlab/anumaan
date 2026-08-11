# Within-/between-cluster same-class agreement and Hamming distance for an arbitrary event grouping (used for both admission and patient clustering).

`cluster_within_group_same_class_agreement`: pooled agreement rate
(fraction of within-group event-pair x class comparisons, restricted to
classes BOTH events tested, where the two values agree).
`cluster_within_group_hamming_distance`: mean, over within-group event
pairs, of that PAIR's own mismatch rate across its jointly-tested
classes (a per-pair-normalised Hamming distance – distinct from the
pooled agreement statistic above, which does not average per pair
first). `cluster_within_between_agreement_diff`: within-group agreement
minus agreement computed over a deterministic random sample of
between-group event pairs (sampling all possible between-group pairs is
infeasible for large N; the sample size scales with the number of
within-group pairs).

## Usage

``` r
.ppc_cluster_pair_stats(
  setup,
  role,
  cluster_id,
  restrict_diff_id = NULL,
  min_group_events = 2L,
  max_between_pairs = 5000L
)
```
