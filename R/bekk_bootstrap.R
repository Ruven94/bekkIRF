#' Bootstrap BEKK model parameters
#'
#' Simulates bootstrap return paths from a fitted BEKK model by resampling the
#' fitted standardized innovations, refits the supplied BEKK specification to
#' each simulated path, and returns bootstrap parameter arrays.
#'
#' @param bekk_model A fitted BEKK model object with entries `H_t`, `data`,
#'   `C0`, and either matrix parameters `A`, `G`, optionally `B`, or
#'   scalar-BEKK parameters `a`, `g`, optionally `b`.
#' @param bekk_spec_model BEKK specification used for re-estimation. If `NULL`,
#'   `bekk_model$spec` is used.
#' @param data Numeric return matrix. If `NULL`, `bekk_model$data` is used.
#' @param bootsamp Number of bootstrap replications.
#' @param cores Number of parallel workers.
#' @param root_type Matrix root used to standardize and simulate innovations.
#'   Use `"spectral"`/`"spec"` or `"cholesky"`/`"chol"`.
#' @param seed Integer seed. Replication `b` uses `seed + b`.
#' @param max_iter Maximum number of BEKK optimizer iterations passed to
#'   [BEKKs::bekk_fit()].
#' @param crit Convergence criterion passed to [BEKKs::bekk_fit()].
#' @param center Logical. If `TRUE`, each bootstrap return sample is centered
#'   column-wise before re-estimation.
#' @param xi_outlier_threshold Numeric threshold used to remove standardized
#'   innovations before resampling. Rows with at least one absolute component
#'   larger than the threshold are skipped. Set to `NULL` to disable filtering.
#' @param progress Logical. If `TRUE`, print a text progress bar.
#' @param chunk_size Number of bootstrap replications submitted per parallel
#'   batch. The default uses `cores`.
#'
#' @returns An object of class `"bekkBootstrap"` containing bootstrap parameter
#'   arrays `C0`, `A`, `G`, optionally `B`, diagnostic vectors `converged`
#'   and `error`, and a `settings` list with bootstrap choices and
#'   standardized-innovation filtering diagnostics. `converged` is a
#'   conservative acceptance flag: estimates must be finite, `BEKK_valid`
#'   must be true, and the optimizer must stop before `max_iter`. BEKKs does
#'   not expose a separate convergence flag, so fits reaching the iteration
#'   limit are excluded even if they might have converged on that iteration.
#'   The `fit_diagnostics` data frame retains iteration counts and validity.
#'
#' @examples
#' \dontrun{
#' data(gold_msci_returns)
#' spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = FALSE))
#' fit <- BEKKs::bekk_fit(spec, gold_msci_returns, max_iter = 50)
#'
#' boot <- bekk_bootstrap(
#'   fit,
#'   bekk_spec_model = spec,
#'   bootsamp = 99,
#'   cores = 2
#' )
#'
#' irf <- compute_irf(
#'   fit,
#'   shock_type = "historical",
#'   time = 444,
#'   bekk_bootstrap = boot
#' )
#' }
#' @export
bekk_bootstrap <- function(bekk_model,
                           bekk_spec_model = NULL,
                           data = NULL,
                           bootsamp = 999,
                           cores = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
                           root_type = c("spectral", "cholesky", "spec", "chol"),
                           seed = 123L,
                           max_iter = 50,
                           crit = 1e-9,
                           center = TRUE,
                           xi_outlier_threshold = 5,
                           progress = TRUE,
                           chunk_size = NULL) {
  if (!requireNamespace("BEKKs", quietly = TRUE)) {
    stop("Package `BEKKs` is required for `bekk_bootstrap()`.")
  }

  root_type <- match.arg(root_type)
  root_type <- switch(
    root_type,
    spectral = "spec",
    cholesky = "chol",
    spec = "spec",
    chol = "chol"
  )

  required_fields <- c("H_t", "C0")
  missing_fields <- setdiff(required_fields, names(bekk_model))
  if (length(missing_fields) > 0) {
    stop(
      "`bekk_model` is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      "."
    )
  }

  if (is.null(data)) {
    if (is.null(bekk_model$data)) {
      stop("`data` must be supplied if `bekk_model$data` is missing.")
    }
    data <- bekk_model$data
  }
  data <- as.matrix(data)

  data_mean <- colMeans(data)
  if (any(abs(data_mean) > 1e-10)) {
    warning(
      "`data` is not centered; column means are ",
      paste(signif(data_mean, 4), collapse = ", "),
      ". Bootstrap re-estimation samples are centered only if `center = TRUE`.",
      call. = FALSE
    )
  }

  if (is.null(bekk_spec_model)) {
    if (is.null(bekk_model$spec)) {
      stop("`bekk_spec_model` must be supplied if `bekk_model$spec` is missing.")
    }
    bekk_spec_model <- bekk_model$spec
  }

  H_t <- as.matrix(bekk_model$H_t)
  N <- nrow(data)
  K <- ncol(data)

  params <- bekk_extract_parameters(bekk_model, K)
  C0 <- params$C0
  A <- params$A
  G <- params$G
  B <- params$B
  asym <- params$asymmetric
  signs <- params$signs

  if (N < 2L) {
    stop("`data` must contain at least two rows.")
  }

  if (nrow(H_t) != N || ncol(H_t) != K^2) {
    stop("`bekk_model$H_t` must have dimensions `nrow(data) x ncol(data)^2`.")
  }

  if (length(bootsamp) != 1L || is.na(bootsamp) || bootsamp <= 0) {
    stop("`bootsamp` must be a positive integer.")
  }
  bootsamp <- as.integer(bootsamp)

  if (length(cores) != 1L || is.na(cores) || cores <= 0) {
    stop("`cores` must be a positive integer.")
  }
  cores <- min(as.integer(cores), bootsamp)

  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a single integer.")
  }
  seed <- as.integer(seed)

  if (length(center) != 1L || is.na(center)) {
    stop("`center` must be TRUE or FALSE.")
  }
  center <- isTRUE(center)

  if (!is.null(xi_outlier_threshold)) {
    if (length(xi_outlier_threshold) != 1L ||
        is.na(xi_outlier_threshold) ||
        xi_outlier_threshold <= 0) {
      stop("`xi_outlier_threshold` must be a positive number or NULL.")
    }
    xi_outlier_threshold <- as.numeric(xi_outlier_threshold)
  }

  if (is.null(chunk_size)) {
    chunk_size <- cores
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size <= 0) {
    stop("`chunk_size` must be a positive integer.")
  }
  chunk_size <- as.integer(chunk_size)

  xi <- compute_xi(H_t, data, root_type = root_type)
  xi_original_n <- nrow(xi)

  if (!is.null(xi_outlier_threshold)) {
    keep_xi <- apply(abs(xi) <= xi_outlier_threshold, 1L, all)
    xi <- xi[keep_xi, , drop = FALSE]

    if (nrow(xi) == 0L) {
      stop("No standardized innovations remain after applying `xi_outlier_threshold`.")
    }
  }

  xi_resample_n <- nrow(xi)
  H_start <- matrix(H_t[1L, ], nrow = K, ncol = K)

  worker_fun <- function(b) {
    set.seed(seed + b)
    idx <- sample.int(xi_resample_n, size = N, replace = TRUE)
    xi_sim <- xi[idx, , drop = FALSE]

    H_prev <- H_start
    ret_sim <- matrix(0, nrow = N, ncol = K)
    ret_sim[1L, ] <- drop(matroot(H_prev, type = root_type) %*% xi_sim[1L, ])

    for (i in 2:N) {
      e_prev <- ret_sim[i - 1L, ]
      H_new <- C0 %*% t(C0) +
        t(A) %*% tcrossprod(e_prev) %*% A +
        t(G) %*% H_prev %*% G

      if (asym) {
        eta <- compute_eta(e_prev, signs)
        H_new <- H_new + t(B) %*% tcrossprod(eta) %*% B
      }

      H_new <- (H_new + t(H_new)) / 2
      ret_sim[i, ] <- drop(matroot(H_new, type = root_type) %*% xi_sim[i, ])
      H_prev <- H_new
    }

    if (center) {
      ret_sim <- scale(ret_sim, center = TRUE, scale = FALSE)
    }

    fit <- tryCatch(
      BEKKs::bekk_fit(bekk_spec_model, ret_sim, max_iter = max_iter, crit = crit),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(
        C0 = matrix(NA_real_, K, K),
        A = matrix(NA_real_, K, K),
        G = matrix(NA_real_, K, K),
        B = if (asym) matrix(NA_real_, K, K) else NULL,
        converged = FALSE,
        error = conditionMessage(fit)
      ))
    }

    fit_params <- tryCatch(
      bekk_extract_parameters(fit, K),
      error = function(e) e
    )

    if (inherits(fit_params, "error")) {
      return(list(
        C0 = matrix(NA_real_, K, K),
        A = matrix(NA_real_, K, K),
        G = matrix(NA_real_, K, K),
        B = if (asym) matrix(NA_real_, K, K) else NULL,
        converged = FALSE,
        error = conditionMessage(fit_params)
      ))
    }

    diagnostics <- bekk_bootstrap_fit_diagnostics(fit, fit_params, max_iter)

    list(
      C0 = fit_params$C0,
      A = fit_params$A,
      G = fit_params$G,
      B = if (asym) fit_params$B else NULL,
      converged = diagnostics$accepted,
      error = diagnostics$error,
      iter = diagnostics$iter,
      valid = diagnostics$valid
    )
  }

  C0_boot <- array(NA_real_, dim = c(K, K, bootsamp))
  A_boot <- array(NA_real_, dim = c(K, K, bootsamp))
  G_boot <- array(NA_real_, dim = c(K, K, bootsamp))
  B_boot <- if (asym) array(NA_real_, dim = c(K, K, bootsamp)) else NULL
  converged <- rep(FALSE, bootsamp)
  error <- rep(NA_character_, bootsamp)
  fit_iter <- rep(NA_real_, bootsamp)
  fit_valid <- rep(NA, bootsamp)

  chunks <- split(seq_len(bootsamp), ceiling(seq_len(bootsamp) / chunk_size))

  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = bootsamp, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  done <- 0L
  cl <- NULL
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "worker_fun", "seed", "N", "K", "xi", "xi_resample_n", "H_start", "root_type",
        "C0", "A", "G", "B", "asym", "signs", "matroot", "compute_eta",
        "bekk_spec_model", "max_iter", "crit", "center",
        "bekk_extract_parameters", "bekk_parameter_matrix", "bekk_model_type",
        "bekk_bootstrap_fit_diagnostics"
      ),
      envir = environment()
    )
  }

  run_chunk <- function(chunk) {
    if (cores == 1L) {
      lapply(chunk, worker_fun)
    } else {
      parallel::parLapplyLB(cl, chunk, worker_fun)
    }
  }

  for (chunk in chunks) {
    chunk_res <- run_chunk(chunk)

    for (j in seq_along(chunk)) {
      b <- chunk[j]
      res <- chunk_res[[j]]
      C0_boot[, , b] <- res$C0
      A_boot[, , b] <- res$A
      G_boot[, , b] <- res$G
      if (asym) {
        B_boot[, , b] <- res$B
      }
      converged[b] <- isTRUE(res$converged)
      error[b] <- res$error
      if (!is.null(res$iter)) fit_iter[b] <- res$iter
      if (!is.null(res$valid)) fit_valid[b] <- res$valid
    }

    done <- done + length(chunk)
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, done)
    }
  }

  if (any(!converged)) {
    warning(
      sum(!converged), " of ", bootsamp,
      " bootstrap fits were excluded; inspect `error` and `fit_diagnostics`.",
      call. = FALSE
    )
  }

  out <- list(
    C0 = C0_boot,
    A = A_boot,
    G = G_boot,
    B = B_boot,
    converged = converged,
    error = error,
    fit_diagnostics = data.frame(iter = fit_iter, valid = fit_valid,
                                 accepted = converged),
    settings = list(
      bootsamp = bootsamp,
      cores = cores,
      root_type = root_type,
      seed = seed,
      max_iter = max_iter,
      crit = crit,
      center = center,
      xi_outlier_threshold = xi_outlier_threshold,
      xi_original_n = xi_original_n,
      xi_resample_n = xi_resample_n,
      xi_removed_n = xi_original_n - xi_resample_n,
      asymmetric = asym,
      model_type = params$model_type
    )
  )

  class(out) <- "bekkBootstrap"
  out
}

bekk_bootstrap_fit_diagnostics <- function(fit, fit_params, max_iter) {
  iter <- if (is.numeric(fit$iter) && length(fit$iter) == 1L &&
              is.finite(fit$iter)) as.numeric(fit$iter) else NA_real_
  valid <- isTRUE(fit$BEKK_valid)
  finite <- all(is.finite(unlist(fit_params[c("C0", "A", "G", "B")],
                                 use.names = FALSE))) &&
    is.numeric(fit$log_likelihood) && length(fit$log_likelihood) == 1L &&
    is.finite(fit$log_likelihood)
  problems <- character()
  if (!valid) problems <- c(problems, "BEKKs did not report a valid model")
  if (!finite) problems <- c(problems, "Non-finite parameter estimates or log likelihood")
  if (is.na(iter)) {
    problems <- c(problems, "Optimizer iteration count is unavailable")
  } else if (iter >= max_iter) {
    problems <- c(problems, "Optimizer reached max_iter; convergence is not confirmed")
  }
  list(accepted = length(problems) == 0L,
       error = if (length(problems)) paste(problems, collapse = "; ") else NA_character_,
       iter = iter, valid = valid)
}
