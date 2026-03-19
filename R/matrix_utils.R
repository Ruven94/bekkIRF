#' Compute a matrix square root
#'
#' Computes either the symmetric square root of a positive semi-definite
#' matrix via spectral decomposition or the lower-triangular Cholesky root
#' of a positive definite matrix.
#'
#' @param mat A symmetric numeric matrix. Must be positive semi-definite
#'   for `type = "spec"` and positive definite for `type = "chol"`.
#' @param type Character string specifying the matrix root. Either
#'   `"spec"` for the symmetric spectral root or `"chol"` for the
#'   Cholesky root.
#' @param tol Numeric tolerance for small negative eigenvalues caused by
#'   numerical imprecision.
#'
#' @returns A matrix square root of `mat` based on the chosen type.
#' @keywords internal

matroot <- function(mat, type = c("spec", "chol"), tol = 1e-12) {
  type <- match.arg(type)

  # --- Checks --------------------------------------------------------------
  if (!is.matrix(mat)) {
    stop("`mat` must be a matrix.")
  }

  if (nrow(mat) != ncol(mat)) {
    stop("`mat` must be square.")
  }

  if (!isTRUE(all.equal(mat, t(mat), tolerance = tol))) {
    stop("`mat` must be symmetric.")
  }

  # --- Root ----------------------------------------------------------------
  if (type == "spec") {
    ev <- eigen(mat, symmetric = TRUE)

    values <- ev$values
    values[values < 0 & values > -tol] <- 0

    if (any(values < 0)) {
      stop("`mat` is not positive semi-definite.")
    }

    root <- ev$vectors %*% diag(sqrt(values), nrow = length(values)) %*% t(ev$vectors)
  }

  if (type == "chol") {
    root <- tryCatch(
      t(chol(mat)),
      error = function(e) {
        stop("`mat` must be positive definite for `type = \"chol\"`.")
      }
    )
  }

  return(root)
}
