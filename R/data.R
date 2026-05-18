#' MSCI World Developed Markets and gold returns
#'
#' A two-column matrix of centered daily log returns used for package examples,
#' validation scripts, and the package application study. The matrix keeps the
#' return dates as row names and in the `"dates"` attribute. The return sample
#' runs from 2007-01-03 to 2025-11-28.
#'
#' @format A numeric matrix with 4,920 dated rows and 2 columns:
#' \describe{
#'   \item{msci.Ret}{MSCI World Developed Markets returns.}
#'   \item{gold.Ret}{Gold returns.}
#' }
#' The matrix has daily trading dates as `rownames(gold_msci_returns)` and as
#' `attr(gold_msci_returns, "dates")`.
#'
#' @source Refinitiv Workspace data processed by the package author.
"gold_msci_returns"
