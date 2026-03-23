# Compute YLDs Attributable to Resistance

Multiplies `YLD_k` (from
[`calculate_YLD()`](https://saketlab.github.io/anumaan/reference/calculate_YLD.md))
by the LOS-based PAF (from
[`compute_paf_los()`](https://saketlab.github.io/anumaan/reference/compute_paf_los.md)).

## Usage

``` r
compute_yld_attributable(
  yld_k_tbl,
  paf_los_list,
  pathogen_col = "pathogen",
  yld_col = "YLD"
)
```

## Arguments

- yld_k_tbl:

  Data frame from
  [`calculate_YLD()`](https://saketlab.github.io/anumaan/reference/calculate_YLD.md)
  containing at least a pathogen column and a YLD column.

- paf_los_list:

  Named list from
  [`compute_paf_los()`](https://saketlab.github.io/anumaan/reference/compute_paf_los.md).

- pathogen_col:

  Character. Pathogen column in `yld_k_tbl`. Default `"pathogen"`.

- yld_col:

  Character. YLD column in `yld_k_tbl`. Default `"YLD"`.

## Value

`yld_k_tbl` augmented with columns `PAF_k`, `denominator`, and
`YLD_attributable`.

## Details

YLD_attributable_k = YLD_k \* PAF_k

Answers: "How much disability burden exists \*only because\* infections
were resistant instead of susceptible?" This is a counterfactual – it
measures the excess burden driven purely by resistance.

Note: YLD_attributable_k \< YLD_associated_k always, because PAF_k =
Fraction_k \* (1 - 1/E_RR_k) \< Fraction_k.
