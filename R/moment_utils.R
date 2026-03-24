#' Compute fourth-order moment matrix
#'
#' Computes the empirical fourth-order moment matrix used in higher-moment
#' impulse response analysis. The matrix is defined as
#' \deqn{\Psi = \frac{1}{T} \sum_{t=1}^T (\xi_t \xi_t') \otimes (\xi_t \xi_t').}
#'
#' @param xi A numeric matrix of dimension `T x K` containing standardized
#'   innovations.
#'
#' @returns A numeric matrix of dimension `K^2 x K^2`.
#' @keywords internal
compute_psi_kurtosis <- function(xi) {

  # --- Checks --------------------------------------------------------------
  if (!is.matrix(xi)) {
    stop("`xi` must be a matrix.")
  }

  if (!is.numeric(xi)) {
    stop("`xi` must be numeric.")
  }

  T <- nrow(xi)
  K <- ncol(xi)

  if (T == 0) {
    stop("`xi` must have at least one row.")
  }

  # --- Computation ---------------------------------------------------------
  Psi_sum <- matrix(0, K * K, K * K)

  for (t in seq_len(T)) {
    xi_t <- xi[t, , drop = FALSE]   # 1 x K
    A  <- t(xi_t) %*% xi_t            # K x K outer product

    Psi_sum <- Psi_sum + kronecker(A, A)
  }

  Psi <- Psi_sum / T
  return(Psi)
}


#' Compute third-order moment matrix
#'
#' Computes the empirical third-order moment matrix used in higher-moment
#' impulse response analysis. The matrix is defined as
#' \deqn{\Psi = \frac{1}{T} \sum_{t=1}^T (\xi_t \xi_t') \otimes \xi_t.}
#'
#' @param xi A numeric matrix of dimension `T x K` containing standardized
#'   innovations.
#'
#' @returns A numeric matrix of dimension `K^2 x K`.
#' @keywords internal
compute_psi_skewness <- function(xi) {

  # --- Checks --------------------------------------------------------------
  if (!is.matrix(xi)) {
    stop("`xi` must be a matrix.")
  }

  if (!is.numeric(xi)) {
    stop("`xi` must be numeric.")
  }

  T <- nrow(xi)
  K <- ncol(xi)

  if (T == 0) {
    stop("`xi` must have at least one row.")
  }

  # --- Computation ---------------------------------------------------------
  Psi_sum <- matrix(0, K * K, K)

  for (t in seq_len(T)) {
    xi_t <- xi[t, , drop = FALSE]   # 1 x K
    A  <- t(xi_t) %*% xi_t            # K x K outer product
    B  <- as.numeric(xi_t)          # K vector

    Psi_sum <- Psi_sum + kronecker(A, B)
  }

  Psi <- Psi_sum / T
  return(Psi)
}
