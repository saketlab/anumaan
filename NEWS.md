# anumaan 0.1.0.9021

* Expanded the specimen classification reference with stewardship sample labels
  seen in 2021-era extracts, including brain abscess, instrument, lung aspirate,
  lymph node, hair, and superficial biopsy.
* Collapsed `sterile_classification` outputs from detailed reference labels to
  `Sterile`, `Non-Sterile`, and `Others/Ambiguous` for downstream analysis.
* Made EDA and burden plotting functions usable without `ggpubr`; when `ggpubr`
  is not installed, `eda_theme()` falls back to `ggplot2::theme_bw()`.
* Refreshed vignettes and installation guidance for the 0.1.0.9021 workflow.
