#' bekkIRF: Impulse Response Functions for BEKK Models
#'
#' `bekkIRF` computes simulation-based impulse response functions for fitted
#' BEKK models, including symmetric and asymmetric specifications estimated
#' with the `BEKKs` package. The package focuses on paired shock and no-shock
#' path simulations and optional bootstrap confidence intervals.
#'
#' @details
#' The usual workflow is:
#'
#' \enumerate{
#'   \item Estimate a BEKK, diagonal BEKK, or scalar BEKK model with
#'     [BEKKs::bekk_fit()].
#'   \item Compute impulse response functions with [compute_irf()].
#'   \item Optionally obtain bootstrap parameter draws with [bekk_bootstrap()]
#'     and pass them back to [compute_irf()] for bootstrap confidence intervals.
#'   \item Inspect results with `print()`, `summary()`, and `plot()`.
#' }
#'
#' The package supports variance impulse response functions (VIRFs),
#' correlation impulse response functions (CIRFs), skewness impulse response
#' functions (SIRFs), kurtosis impulse response functions (KIRFs), and weights
#' impulse response functions (WIRFs) for optimal portfolio weights whenever
#' the corresponding quantities are requested in [compute_irf()].
#'
#' @references
#' Hafner, C. M. and Herwartz, H. (2006). Volatility impulse responses for
#' multivariate GARCH models: An exchange rate illustration. *Journal of
#' International Money and Finance*, 25(5), 719-740.
#'
#' Hafner, C. M. and Herwartz, H. (2023). Correlation impulse response
#' functions.
#'
#' Hafner, C. M. and Herwartz, H. (2023). Asymmetric volatility impulse
#' response functions.
#'
#' @seealso [compute_irf()], [bekk_bootstrap()], [plot.bekkIRF()],
#'   [gold_msci_returns]
#' @author
#' Ruven Zapf \email{Ruven.Zapf@uni-goettingen.de}
#'
#' Helmut Herwartz \email{hherwartz@uni-goettingen.de}
#' @useDynLib bekkIRF, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

utils::globalVariables(c("horizon", "lower", "upper"))
