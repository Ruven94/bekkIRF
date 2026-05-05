#' Bootstrap BEKK model parameters
#'
#' Simulates bootstrap return paths from a fitted BEKK model by resampling the
#' fitted standardized innovations, refits the supplied BEKK specification to
#' each simulated path, and returns bootstrap parameter arrays.
#'
#' @param bekk_model A fitted BEKK model object with entries `H_t`, `data`,
#'   `C0`, `A`, `G`, and optionally `B`, `asymmetric`, `signs`, and `spec`.
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
#' @param progress Logical. If `TRUE`, print a text progress bar.
#' @param chunk_size Number of bootstrap replications submitted per parallel
#'   batch. The default uses `cores`.
#'
#' @returns An object of class `"bekkBootstrap"` containing bootstrap parameter
#'   arrays `C0`, `A`, `G`, optionally `B`, and diagnostic vectors `converged`
#'   and `error`.
#' @export
bekk_bootstrap <- function(bekk_model,
                           bekk_spec_model = NULL,
                           data = NULL,
                           bootsamp = 999,
                           cores = max(1L, parallel::detectCores() - 1L),
                           root_type = c("spectral", "cholesky", "spec", "chol"),
                           seed = 123L,
                           max_iter = 50,
                           crit = 1e-9,
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

  required_fields <- c("H_t", "C0", "A", "G")
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

  if (is.null(bekk_spec_model)) {
    if (is.null(bekk_model$spec)) {
      stop("`bekk_spec_model` must be supplied if `bekk_model$spec` is missing.")
    }
    bekk_spec_model <- bekk_model$spec
  }

  H_t <- as.matrix(bekk_model$H_t)
  C0 <- as.matrix(bekk_model$C0)
  A <- as.matrix(bekk_model$A)
  G <- as.matrix(bekk_model$G)

  N <- nrow(data)
  K <- ncol(data)

  if (N < 2L) {
    stop("`data` must contain at least two rows.")
  }

  if (nrow(H_t) != N || ncol(H_t) != K^2) {
    stop("`bekk_model$H_t` must have dimensions `nrow(data) x ncol(data)^2`.")
  }

  if (nrow(C0) != K || ncol(C0) != K ||
      nrow(A) != K || ncol(A) != K ||
      nrow(G) != K || ncol(G) != K) {
    stop("`C0`, `A`, and `G` must be square matrices compatible with `data`.")
  }

  asym <- isTRUE(bekk_model$asymmetric) || !is.null(bekk_model$B)
  if (asym) {
    if (is.null(bekk_model$B)) {
      stop("`bekk_model$B` is required for asymmetric BEKK bootstrap.")
    }
    B <- as.matrix(bekk_model$B)
    if (nrow(B) != K || ncol(B) != K) {
      stop("`B` must be a square matrix compatible with `data`.")
    }
  } else {
    B <- matrix(0, K, K)
  }

  signs <- bekk_model$signs
  if (is.null(signs)) {
    signs <- rep(-1, K)
  }
  signs <- as.numeric(signs)
  if (length(signs) != K || !all(signs %in% c(-1, 1))) {
    stop("`signs` must contain one -1 or 1 value per series.")
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

  if (is.null(chunk_size)) {
    chunk_size <- cores
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size <= 0) {
    stop("`chunk_size` must be a positive integer.")
  }
  chunk_size <- as.integer(chunk_size)

  xi <- compute_xi(H_t, data, root_type = root_type)
  H_start <- matrix(H_t[1L, ], nrow = K, ncol = K)

  worker_fun <- function(b) {
    set.seed(seed + b)
    idx <- sample.int(N, size = N, replace = TRUE)
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

    list(
      C0 = fit$C0,
      A = fit$A,
      G = fit$G,
      B = if (asym) fit$B else NULL,
      converged = TRUE,
      error = NA_character_
    )
  }

  C0_boot <- array(NA_real_, dim = c(K, K, bootsamp))
  A_boot <- array(NA_real_, dim = c(K, K, bootsamp))
  G_boot <- array(NA_real_, dim = c(K, K, bootsamp))
  B_boot <- if (asym) array(NA_real_, dim = c(K, K, bootsamp)) else NULL
  converged <- rep(FALSE, bootsamp)
  error <- rep(NA_character_, bootsamp)

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
        "worker_fun", "seed", "N", "K", "xi", "H_start", "root_type",
        "C0", "A", "G", "B", "asym", "signs", "matroot", "compute_eta",
        "bekk_spec_model", "max_iter", "crit"
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
    }

    done <- done + length(chunk)
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, done)
    }
  }

  out <- list(
    C0 = C0_boot,
    A = A_boot,
    G = G_boot,
    B = B_boot,
    converged = converged,
    error = error,
    settings = list(
      bootsamp = bootsamp,
      cores = cores,
      root_type = root_type,
      seed = seed,
      max_iter = max_iter,
      crit = crit,
      asymmetric = asym
    )
  )

  class(out) <- "bekkBootstrap"
  out
}
