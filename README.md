# RCtest
R package for forecast evaluation
This is a mirror of the CRAN R package https://cran.r-project.org/web/packages/RCtest/index.html
## Installation

```r
# From CRAN (stable):
install.packages("RCtest")

# From GitHub (development version):
# install.packages("remotes")
remotes::install_github("YOUR_USERNAME/RCtest")
```

## Tests implemented

| Test | Function | Reference |
|------|----------|-----------|
| White's Reality Check | `white_reality_check()` | White (2000) |
| Superior Predictive Ability | `superior_predictive_ability_test()` | Hansen (2005) |
| Conditional Predictive Ability | `white_reality_check_conditional()` | Giacomini & White (2006) |
| Expected Loss CDF Comparison | `white_reality_check_cdf_approx()` | Corradi & Swanson (2006) |
| KLIC density test | `kullback_leibler_test()` | Corradi & Swanson (2006) |
| ZP Quantile Loss | `reality_check_zp_test()` | Corradi & Swanson (2006) |
| CRPS computation | `compute_crps()` | Gneiting & Raftery (2007) |
| Kupiec UC coverage | `compute_kupiec()` | Kupiec (1995) |
| Diebold–Mariano (per model) | `compute_per_model_statistics()` | Diebold & Mariano (1995) |

All bootstrap inference uses the **Moving Block Bootstrap** (`mbb_resample_data()`) with **Newey–West HAC covariance** (`estimate_long_run_covariance()`).

## Dataset

The package ships with the `metals` dataset: a 165 × 15 matrix of monthly World Bank base metals price index forecasts (aluminium, copper, lead, nickel, tin, and zinc; base year 2010 = 100, USD) covering **March 2011 – November 2024**. Columns 1–14 contain point forecasts from 14 competing Bayesian and shrinkage-based models — variants of Bayesian Dynamic Mixture Models (BDMM), Dynamic Model Averaging (DMA), Time-Varying Parameter (TVP) models, Bayesian LASSO, Bayesian RIDGE, AR(1), and Linear Regression. Column 15 contains the **Historical Average (HA) benchmark**. The dataset is described in Drachal & Jędrzejewska (2025); underlying raw data are available on [figshare](https://doi.org/10.6084/m9.figshare.28382480.v1).

## Funding
This research was funded by the National Science Centre, Poland, grant no. 2022/45/B/HS4/00510.

## License
GPL-3
