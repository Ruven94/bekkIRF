#' Compute theoretical symmetric BEKK VIRFs
#'
#' Computes the closed-form volatility impulse response function for a
#' symmetric BEKK model as in Hafner and Herwartz (2006). The function supports
#' structural standardized shocks supplied by the user and empirical shocks
#' obtained from the fitted model's standardized innovations.
#'
#' @param bekk_model A fitted symmetric BEKK model object with entries `H_t`,
#'   `data`, and either matrix parameters `A`, `G`, or scalar-BEKK parameters
#'   `a`, `g`.
#' @param root_type Matrix root used to map standardized shocks into returns.
#'   Use `"spectral"`/`"spec"` or `"cholesky"`/`"chol"`.
#' @param shock_type Either `"structural"` for a user supplied standardized
#'   shock or `"empirical"` to use `xi[time, ]`.
#' @param shock Numeric standardized shock vector. Required for
#'   `shock_type = "structural"`.
#' @param time Initial time index used for `H_t` and empirical shocks.
#' @param n.ahead Number of horizons.
#' @param format Output format. `"vech"` returns the lower triangular part of
#'   each VIRF matrix, while `"array"` returns the full `K x K x n.ahead` array.
#'
#' @returns A matrix of dimension `n.ahead x K * (K + 1) / 2` for
#'   `format = "vech"` or an array of dimension `K x K x n.ahead` for
#'   `format = "array"`.
#' @keywords internal
compute_theoretical_virf <- function(bekk_model,
                                     root_type = c("spectral", "cholesky", "spec", "chol"),
                                     shock_type = c("structural", "empirical"),
                                     shock = NULL,
                                     time = NULL,
                                     n.ahead = 100,
                                     format = c("vech", "array")) {
  root_type <- match.arg(root_type)
  shock_type <- match.arg(shock_type)
  format <- match.arg(format)

  root_type <- switch(
    root_type,
    spectral = "spec",
    cholesky = "chol",
    spec = "spec",
    chol = "chol"
  )

  if (isTRUE(bekk_model$asymmetric)) {
    stop("Closed-form theoretical VIRFs are implemented only for symmetric BEKK models.")
  }

  required_fields <- c("H_t", "data")
  missing_fields <- setdiff(required_fields, names(bekk_model))
  if (length(missing_fields) > 0) {
    stop(
      "`bekk_model` is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      "."
    )
  }

  data <- as.matrix(bekk_model$data)
  H_t <- as.matrix(bekk_model$H_t)
  N <- nrow(data)
  K <- ncol(data)
  A <- bekk_parameter_matrix(bekk_model, "A", "a", K)
  G <- bekk_parameter_matrix(bekk_model, "G", "g", K)

  if (nrow(H_t) != N || ncol(H_t) != K^2) {
    stop("`bekk_model$H_t` must have dimensions `nrow(data) x ncol(data)^2`.")
  }

  if (nrow(A) != K || ncol(A) != K || nrow(G) != K || ncol(G) != K) {
    stop("`bekk_model$A` and `bekk_model$G` must be square matrices compatible with `data`.")
  }

  if (is.null(time) || length(time) != 1L || is.na(time)) {
    stop("`time` must be a single valid row index.")
  }

  time <- as.integer(time)
  if (time < 1L || time > N) {
    stop("`time` must be between 1 and `nrow(bekk_model$data)`.")
  }

  if (length(n.ahead) != 1L || is.na(n.ahead) || n.ahead <= 0) {
    stop("`n.ahead` must be a positive integer.")
  }
  n.ahead <- as.integer(n.ahead)

  if (shock_type == "empirical") {
    xi <- compute_xi(H_t, data, root_type = root_type)
    shock <- xi[time, ]
  } else if (is.null(shock)) {
    stop("`shock` must be supplied for `shock_type = \"structural\"`.")
  }

  shock <- as.numeric(shock)
  if (length(shock) != K) {
    stop("`shock` must have length equal to `ncol(bekk_model$data)`.")
  }

  H_0 <- matrix(H_t[time, ], nrow = K, ncol = K)
  Q_0 <- matroot(H_0, type = root_type)

  shock_outer <- tcrossprod(Q_0 %*% shock)
  virf_current <- t(A) %*% (shock_outer - H_0) %*% A

  out <- array(0, dim = c(K, K, n.ahead))
  out[, , 1L] <- virf_current

  if (n.ahead > 1L) {
    for (h in 2:n.ahead) {
      virf_current <- t(A) %*% virf_current %*% A +
        t(G) %*% virf_current %*% G
      out[, , h] <- virf_current
    }
  }

  if (format == "array") {
    return(out)
  }

  vech_idx <- lower.tri(matrix(NA, K, K), diag = TRUE)
  virf_vech <- matrix(0, nrow = n.ahead, ncol = K * (K + 1) / 2)

  for (h in seq_len(n.ahead)) {
    virf_vech[h, ] <- out[, , h][vech_idx]
  }

  colnames(virf_vech) <- apply(which(vech_idx, arr.ind = TRUE), 1L, function(idx) {
    paste0("H", idx[1L], idx[2L])
  })

  virf_vech
}
