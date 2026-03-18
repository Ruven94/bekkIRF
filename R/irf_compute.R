#' Compute an impulse response function
#'
#' @param x A numeric input
#'
#' @returns A numeric value.
#' @export
#'
#' @examples
#' x <- 3
#' irf_compute(x)

irf_compute <- function(x){
  return(x^2 + x)
}
