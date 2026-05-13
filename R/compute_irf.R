
#' Compute simulation-based BEKK impulse response functions
#'
#' Simulates paired non-shock and shock BEKK paths from one initial covariance
#' matrix and averages the resulting differences over Monte Carlo paths.
#'
#' @param bekk_model A fitted BEKK model object with entries `H_t`, `data`,
#'   `C0`, `A`, `G`, and optionally `B`, `asymmetric`, and `signs`.
#' @param root_type Matrix root used to map standardized innovations into
#'   shocks. Use `"spectral"`/`"spec"` or `"cholesky"`/`"chol"`.
#' @param shock_type Either `"structural"` for a user supplied standardized
#'   shock or `"empirical"` to use `xi[time, ]`.
#' @param shock Numeric standardized shock vector. Required for
#'   `shock_type = "structural"`.
#' @param time Initial time index used for `H_0` and empirical shocks.
#' @param simsamp Number of Monte Carlo paths.
#' @param n.ahead Number of horizons.
#' @param seed Integer random seed.
#' @param calc_virf,calc_cirf,calc_kirf,calc_sirf,calc_wirf Logical flags for
#'   the requested IRF types.
#' @param bekk_bootstrap Optional `"bekkBootstrap"` object. If supplied,
#'   `compute_irf()` recomputes the selected IRFs for each converged bootstrap
#'   parameter draw and returns pointwise bootstrap confidence intervals.
#' @param ci_level Confidence level for central bootstrap intervals, e.g.
#'   `0.95`, or a length-two numeric vector of quantile probabilities.
#' @param bootstrap_cores Number of parallel workers used for bootstrap IRF
#'   paths. Ignored if `bekk_bootstrap = NULL`.
#' @param bootstrap_progress Logical. If `TRUE`, print a text progress bar for
#'   bootstrap IRF paths.
#' @param bootstrap_chunk_size Number of bootstrap IRF paths submitted per
#'   parallel batch. The default uses `bootstrap_cores`.
#'
#' @returns A list with mean IRF matrices for the selected IRF types.
#' @export
compute_irf <- function(bekk_model,
                        root_type = c("spectral", "cholesky", "spec", "chol"),
                        shock_type = c("structural", "empirical"),
                        shock = NULL,
                        time = NULL,
                        simsamp = 10000,
                        n.ahead = 100,
                        seed = 123,
                        calc_virf = TRUE,
                        calc_cirf = TRUE,
                        calc_kirf = FALSE,
                        calc_sirf = FALSE,
                        calc_wirf = FALSE,
                        bekk_bootstrap = NULL,
                        ci_level = 0.95,
                        bootstrap_cores = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
                        bootstrap_progress = TRUE,
                        bootstrap_chunk_size = NULL) {
  root_type <- match.arg(root_type)
  shock_type <- match.arg(shock_type)

  root_type <- switch(
    root_type,
    spectral = "spec",
    cholesky = "chol",
    spec = "spec",
    chol = "chol"
  )
  root_type_id <- switch(root_type, spec = 0L, chol = 1L)

  required_fields <- c("H_t", "data", "C0", "A", "G")
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
  K <- ncol(data)
  N <- nrow(data)

  if (is.null(time) || length(time) != 1L || is.na(time)) {
    stop("`time` must be a single valid row index.")
  }

  time <- as.integer(time)
  if (time < 1L || time > N) {
    stop("`time` must be between 1 and `nrow(bekk_model$data)`.")
  }

  if (simsamp <= 0) {
    stop("`simsamp` must be positive.")
  }

  if (n.ahead <= 0) {
    stop("`n.ahead` must be positive.")
  }

  xi <- compute_xi(H_t, data, root_type = root_type)
  psi_skew <- compute_psi_skewness(xi)
  psi_kurt <- compute_psi_kurtosis(xi)

  C <- as.matrix(bekk_model$C0)
  A <- as.matrix(bekk_model$A)
  G <- as.matrix(bekk_model$G)

  asym <- isTRUE(bekk_model$asymmetric)
  if (asym) {
    if (is.null(bekk_model$B)) {
      stop("`bekk_model$B` is required for asymmetric BEKK models.")
    }
    B <- as.matrix(bekk_model$B)
  } else {
    B <- matrix(0, K, K)
  }

  signs <- bekk_model$signs
  if (is.null(signs)) {
    signs <- rep(-1, K)
  }
  signs <- as.numeric(signs)

  H_0 <- matrix(H_t[time, ], nrow = K, ncol = K)

  if (shock_type == "empirical") {
    shock <- xi[time, ]
  } else if (is.null(shock)) {
    stop("`shock` must be supplied for `shock_type = \"structural\"`.")
  }

  shock <- as.numeric(shock)

  out <- compute_irf_core_cpp(
    H_0 = H_0,
    shock = shock,
    xi = xi,
    C = C,
    A = A,
    G = G,
    B = B,
    signs = signs,
    psi_kurt = psi_kurt,
    psi_skew = psi_skew,
    root_type = root_type_id,
    asym = asym,
    simsamp = as.integer(simsamp),
    n_ahead = as.integer(n.ahead),
    seed = as.integer(seed),
    calc_virf = isTRUE(calc_virf),
    calc_cirf = isTRUE(calc_cirf),
    calc_kirf = isTRUE(calc_kirf),
    calc_sirf = isTRUE(calc_sirf),
    calc_wirf = isTRUE(calc_wirf)
  )

  class(out) <- c("bekkIRF", "list")

  if (is.null(bekk_bootstrap)) {
    return(out)
  }

  ci_probs <- compute_irf_ci_probs(ci_level)
  bootstrap_irf <- compute_irf_bootstrap_paths(
    bekk_bootstrap = bekk_bootstrap,
    main_out = out,
    H_0 = H_0,
    shock = shock,
    xi = xi,
    signs = signs,
    psi_kurt = psi_kurt,
    psi_skew = psi_skew,
    root_type_id = root_type_id,
    asym = asym,
    simsamp = as.integer(simsamp),
    n_ahead = as.integer(n.ahead),
    seed = as.integer(seed),
    calc_virf = isTRUE(calc_virf),
    calc_cirf = isTRUE(calc_cirf),
    calc_kirf = isTRUE(calc_kirf),
    calc_sirf = isTRUE(calc_sirf),
    calc_wirf = isTRUE(calc_wirf),
    cores = bootstrap_cores,
    progress = bootstrap_progress,
    chunk_size = bootstrap_chunk_size
  )

  out$bootstrap_irf <- bootstrap_irf$paths
  out$ci <- compute_irf_bootstrap_ci(bootstrap_irf$paths, ci_probs)
  out$bootstrap_info <- list(
    ci_probs = ci_probs,
    ci_level = if (length(ci_level) == 1L) ci_level else NA_real_,
    requested_replications = bootstrap_irf$requested_replications,
    candidate_replications = bootstrap_irf$candidate_replications,
    used_replications = bootstrap_irf$used_replications,
    used_indices = bootstrap_irf$used_indices,
    skipped_indices = bootstrap_irf$skipped_indices,
    errors = bootstrap_irf$errors,
    common_random_numbers = TRUE,
    cores = bootstrap_irf$cores,
    chunk_size = bootstrap_irf$chunk_size
  )

  out
}
