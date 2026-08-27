# RCtest
R package for forecast evaluation
This repository contains the development version of the RCtest package.
The stable release is available on CRAN: <https://cran.r-project.org/package=RCtest>.

## Installation

```r
# From CRAN:
install.packages("RCtest")

# From GitHub:
# install.packages("remotes")
remotes::install_github("JoannaJedrzejewska/RCtest")
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

The package ships with the `metals` dataset: a 165 × 15 matrix of monthly World Bank base metals price index forecasts (aluminium, copper, lead, nickel, tin, and zinc; base year 2010 = 100, USD) covering **March 2011 – November 2024**. The **HA** column (Historical Average) contains the observed base-metals price-index series and is used as the realized outcome in the package examples. The remaining 14 columns contain point forecasts from Bayesian Dynamic Mixture Model (BDMM) variants, Dynamic Model Averaging (DMA), Time-Varying Parameter (TVP) models, Bayesian LASSO, Bayesian RIDGE, AR(1), and linear regression. In the illustrative forecast-comparison analysis, the **AR(1)** column is the **benchmark** forecast and the other 13 forecast columns are competing models. The dataset is described in Drachal & Jędrzejewska (2025); underlying raw data are available on [figshare](https://doi.org/10.6084/m9.figshare.28382480.v1).

## Funding
This research was funded by the National Science Centre, Poland, grant no. 2022/45/B/HS4/00510.

## License
GPL-3
