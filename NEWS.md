# RCtest 1.2.

* Corrected `compute_kupiec()` so that VaR violations are evaluated against
  an explicitly supplied `realized` outcome series rather than against the
  benchmark forecast column.

* Updated `compute_kupiec()` documentation and examples.

* Updated `run_comprehensive_erc_analysis()`.

* Added direct numerical validation of the Kupiec LR-UC statistic against
  an independent implementation.

* Added regression tests for forecast evaluation.

* Added reproducibility scripts.


# RCtest 1.1.

* Added a required `realized` argument to `estimate_forecast_variance()`.

* Corrected variance scaling in SPA, CPA, CDF-based Reality Check, KLIC,
  and ZP procedures.

* Corrected bootstrap recentering in
  `white_reality_check_cdf_approx()`.

* Corrected one-sided alternative directions in
  `compute_per_model_statistics()`.

* Updated package documentation and description text.


# RCtest 1.0.

* Initial CRAN release.