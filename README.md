<!-- README.md is generated from README.Rmd. Please edit that file -->

# bekkIRF

`bekkIRF` computes simulation-based impulse response functions for fitted
BEKK models. It is designed as a companion package for BEKK impulse response
analysis with symmetric and asymmetric volatility dynamics.

The package currently supports:

- simulation-based variance impulse response functions (VIRFs), correlation
  impulse response functions (CIRFs), skewness impulse response functions
  (SIRFs), kurtosis impulse response functions (KIRFs), and weights impulse
  response functions (WIRFs) for optimal portfolio weights,
- empirical and structural shocks,
- spectral and Cholesky matrix roots,
- BEKK, diagonal BEKK, and scalar BEKK outputs from `BEKKs`,
- parameter bootstrap based confidence intervals,
- S3 `print()`, `summary()`, and `plot()` methods for IRF objects.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("pak")
pak::pak("Ruven94/bekkIRF")
```

## Basic workflow

Estimate a BEKK model with `BEKKs`, compute IRFs with `compute_irf()`,
and inspect the result with `summary()` or `plot()`.

``` r
library(bekkIRF)
library(BEKKs)

data(gold_msci_returns)

spec <- BEKKs::bekk_spec(
  model = list(type = "bekk", asymmetric = TRUE)
)

fit <- BEKKs::bekk_fit(spec, gold_msci_returns)

irf <- compute_irf(
  fit,
  shock_type = "empirical",
  time = 444,
  root_type = "spectral",
  simsamp = 10000,
  n.ahead = 100,
  calc_virf = TRUE,
  calc_cirf = TRUE,
  calc_sirf = TRUE,
  calc_kirf = TRUE,
  calc_wirf = TRUE
)

print(irf)
summary(irf)
plot(irf, type = "VIRF")
```

## Bootstrap confidence intervals

Bootstrap confidence intervals are obtained in two steps. First, use
`bekk_bootstrap()` to generate bootstrap parameter draws. Then pass the
resulting `"bekkBootstrap"` object to `compute_irf()`.

``` r
boot <- bekk_bootstrap(
  fit,
  bekk_spec_model = spec,
  bootsamp = 999,
  cores = parallel::detectCores() - 1
)

irf_boot <- compute_irf(
  fit,
  shock_type = "empirical",
  time = 444,
  simsamp = 10000,
  n.ahead = 100,
  bekk_bootstrap = boot,
  ci_level = 0.95
)

plot(irf_boot, type = "CIRF", ci = TRUE)
```

## Package data

The package ships a small empirical data set:

``` r
library(bekkIRF)
data(gold_msci_returns)
head(gold_msci_returns)
```

`gold_msci_returns` contains centered daily log returns for MSCI World
Developed Markets and gold. The raw series were obtained from Refinitiv
Workspace and processed by the package author. The return sample runs from
2007-01-03 to 2025-11-28. Return dates are available as row names and via:

``` r
range(attr(gold_msci_returns, "dates"))
```

## References

Hafner, C. M. and Herwartz, H. (2006). Volatility impulse responses for
multivariate GARCH models: An exchange rate illustration. *Journal of
International Money and Finance*, 25(5), 719-740.

Hafner, C. M. and Herwartz, H. (2023). Correlation impulse response functions.

Hafner, C. M. and Herwartz, H. (2023). Asymmetric volatility impulse response
functions.
