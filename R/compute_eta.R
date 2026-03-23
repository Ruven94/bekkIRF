#' Compute asymmetric return component
#'
#' Computes the asymmetric return component used in the asymmetric
#' BEKK recursion. The return vector is retained only if the sign pattern
#' specified by `signs` is jointly satisfied; otherwise a zero vector is
#' returned.
#'
#' @param e_t A numeric vector of returns at time t.
#' @param signs A numeric vector of length `length(e_t)` consisting only of
#'   `1` and `-1`. A value of `1` requires a positive return, while `-1`
#'   requires a negative return.
#'
#' @returns A numeric vector of the same length as `e_t` representing the
#'   asymmetric return component.
#' @keywords internal

compute_eta <- function(e_t, signs) {
  if (!is.numeric(e_t)) {
    stop("`e_t` must be numeric.")
  }

  if (!is.numeric(signs)) {
    stop("`signs` must be numeric.")
  }

  if (length(signs) != length(e_t)) {
    stop("`signs` must have the same length as `e_t`.")
  }

  if (!all(signs %in% c(-1, 1))) {
    stop("`signs` must contain only -1 and 1.")
  }

  if (all(signs * e_t > 0)) {
    return(e_t)
  }

  rep(0, length(e_t))
}
