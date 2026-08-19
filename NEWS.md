# RCtest 1.1.0

* `estimate_forecast_variance()` now requires a `realized` argument.
  Forecast errors are computed relative to the realized outcome instead
  of the benchmark forecast. 
* Fixed a redundant variance division in `superior_predictive_ability_test()`,
  `white_reality_check_conditional()`, `white_reality_check_cdf_approx()`,
  `kullback_leibler_test()`, and `reality_check_zp_test()`, which inflated
  test statistics by a factor of roughly sqrt(sample size).
* Fixed `white_reality_check_cdf_approx()`'s bootstrap null-cenetring,
  which previously caused the test to reject H0 in almost all samples
  regardless of the data. This is a separate fix from the one above.
* Fixed the one-sided tail direction in `compute_per_model_statistics()`
  (`H1 = "more"`/`"less"`). Two-sided results are unaffected.
* Minor editorial corrections to package description text.