# Compute YLDs Associated with Resistance

Multiplies `YLD_k` (from
[`calculate_YLD()`](https://saketlab.github.io/anumaan/reference/calculate_YLD.md))
by the associated burden score/fraction (from
[`compute_fraction_associated()`](https://saketlab.github.io/anumaan/reference/compute_fraction_associated.md)).

## Usage

``` r
compute_yld_associated(
  yld_k_tbl,
  fraction_assoc_list,
  pathogen_col = "pathogen",
  yld_col = "YLD",
  probability_col = "probability",
  rr_profile_col = "RR_LOS_profile"
)
```

## Arguments

- yld_k_tbl:

  Data frame from
  [`calculate_YLD()`](https://saketlab.github.io/anumaan/reference/calculate_YLD.md)
  containing at least a pathogen column and a YLD column.

- fraction_assoc_list:

  Named list from
  [`compute_fraction_associated()`](https://saketlab.github.io/anumaan/reference/compute_fraction_associated.md).

- pathogen_col:

  Character. Pathogen column in `yld_k_tbl`. Default `"pathogen"`.

- yld_col:

  Character. YLD column in `yld_k_tbl`. Default `"YLD"`.

- probability_col:

  Character. Profile probability column in `per_profile`. Default
  `"probability"`.

- rr_profile_col:

  Character. Profile LOS/RR column in `per_profile`. Default
  `"RR_LOS_profile"`.

## Value

`yld_k_tbl` augmented with columns `Fraction_k`, `R_k_delta`,
`LOS_k_delta`, and `YLD_associated`.

## Details

YLD_associated_k = YLD_k \* Fraction_k
