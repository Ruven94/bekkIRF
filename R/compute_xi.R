#' Compute standardized innovations
#'
#' Computes the standardized innovation vectors \eqn{\xi_t} from a sequence
#' of observations and conditional covariance matrices.
#'
#' For each time point \eqn{t}, the function solves
#' \deqn{Q_t \xi_t = e_t,}
#' where \eqn{e_t} is the observed return vector and \eqn{Q_t} is a matrix
#' root of the conditional covariance matrix \eqn{H_t}. The matrix root is
#' obtained via [matroot()].
#'
#' @param H_t A numeric matrix of dimension `N x K^2` containing the stacked
#'   conditional covariance matrices.
#' @param data A numeric matrix of dimension `N x K` containing centered
#'   return observations.
#' @param root_type Character string specifying the matrix root used in
#'   [matroot()]. Either `"spec"` or `"chol"`.
#'
#' @returns A numeric matrix of dimension `N x K` containing the standardized
#'   innovations.
#' @keywords internal

compute_xi <- function(H_t, data, root_type = c("spec", "chol")) {
  root_type <- match.arg(root_type)

  # --- Checks --------------------------------------------------------------
  if (!is.matrix(H_t)) {
    stop("`H_t` must be a matrix.")
  }

  if (!is.matrix(data)) {
    stop("`data` must be a matrix.")
  }

  N <- nrow(data)
  K <- ncol(data)

  if (nrow(H_t) != N) {
    stop("`H_t` and `data` must have the same number of rows.")
  }

  if (ncol(H_t) != K^2) {
    stop("`H_t` must have `ncol(data)^2` columns.")
  }

  # --- Compute xi ----------------------------------------------------------
  xi <- matrix(0, nrow = N, ncol = K)

  for (i in seq_len(N)) {
    H_i <- matrix(H_t[i, ], nrow = K, ncol = K)
    Q_i <- matroot(H_i, type = root_type)

    xi[i, ] <- drop(solve(Q_i, data[i, ]))
  }

  return(xi)
}
