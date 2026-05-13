compute_irf_ci_probs <- function(ci_level) {
  if (!is.numeric(ci_level) || anyNA(ci_level)) {
    stop("`ci_level` must be numeric.")
  }

  if (length(ci_level) == 1L) {
    if (ci_level <= 0 || ci_level >= 1) {
      stop("`ci_level` must be between 0 and 1.")
    }
    alpha <- 1 - ci_level
    return(c(alpha / 2, 1 - alpha / 2))
  }

  if (length(ci_level) == 2L) {
    ci_probs <- sort(ci_level)
    if (ci_probs[1L] < 0 || ci_probs[2L] > 1 || ci_probs[1L] >= ci_probs[2L]) {
      stop("Quantile probabilities in `ci_level` must satisfy 0 <= lower < upper <= 1.")
    }
    return(ci_probs)
  }

  stop("`ci_level` must be a confidence level or a length-two probability vector.")
}

compute_irf_bootstrap_paths <- function(bekk_bootstrap,
                                        main_out,
                                        H_0,
                                        shock,
                                        xi,
                                        signs,
                                        psi_kurt,
                                        psi_skew,
                                        root_type_id,
                                        asym,
                                        simsamp,
                                        n_ahead,
                                        seed,
                                        calc_virf,
                                        calc_cirf,
                                        calc_kirf,
                                        calc_sirf,
                                        calc_wirf,
                                        cores,
                                        progress,
                                        chunk_size) {
  if (!inherits(bekk_bootstrap, "bekkBootstrap")) {
    stop("`bekk_bootstrap` must be an object of class \"bekkBootstrap\".")
  }

  required_fields <- c("C0", "A", "G")
  missing_fields <- setdiff(required_fields, names(bekk_bootstrap))
  if (length(missing_fields) > 0) {
    stop(
      "`bekk_bootstrap` is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      "."
    )
  }

  if (length(dim(bekk_bootstrap$C0)) != 3L ||
      length(dim(bekk_bootstrap$A)) != 3L ||
      length(dim(bekk_bootstrap$G)) != 3L) {
    stop("`bekk_bootstrap$C0`, `$A`, and `$G` must be three-dimensional arrays.")
  }

  K <- nrow(H_0)
  if (is.null(K) || K != ncol(H_0)) {
    stop("`H_0` must be a square covariance matrix.")
  }

  boot_dim <- dim(bekk_bootstrap$C0)
  R <- boot_dim[3L]
  if (!identical(dim(bekk_bootstrap$A), boot_dim) ||
      !identical(dim(bekk_bootstrap$G), boot_dim)) {
    stop("`bekk_bootstrap$C0`, `$A`, and `$G` must have identical dimensions.")
  }

  if (!identical(boot_dim[1:2], c(K, K))) {
    stop("Bootstrap parameter dimensions are not compatible with `bekk_model`.")
  }

  if (asym) {
    if (is.null(bekk_bootstrap$B) || length(dim(bekk_bootstrap$B)) != 3L) {
      stop("`bekk_bootstrap$B` is required for asymmetric BEKK IRF bootstrap.")
    }
    if (!identical(dim(bekk_bootstrap$B), boot_dim)) {
      stop("`bekk_bootstrap$B` must have the same dimensions as `$C0`.")
    }
  }

  converged <- bekk_bootstrap$converged
  if (is.null(converged)) {
    converged <- rep(TRUE, R)
  }
  if (length(converged) != R) {
    stop("`bekk_bootstrap$converged` must have length equal to the number of bootstrap draws.")
  }

  param_finite <- vapply(seq_len(R), function(b) {
    ok <- all(is.finite(bekk_bootstrap$C0[, , b])) &&
      all(is.finite(bekk_bootstrap$A[, , b])) &&
      all(is.finite(bekk_bootstrap$G[, , b]))
    if (asym) {
      ok <- ok && all(is.finite(bekk_bootstrap$B[, , b]))
    }
    ok
  }, logical(1L))

  candidate_indices <- which(converged & param_finite)
  if (length(candidate_indices) == 0L) {
    stop("No converged finite bootstrap parameter draws are available.")
  }

  if (length(cores) != 1L || is.na(cores) || cores <= 0) {
    stop("`bootstrap_cores` must be a positive integer.")
  }
  cores <- min(as.integer(cores), length(candidate_indices))

  if (length(progress) != 1L || is.na(progress)) {
    stop("`bootstrap_progress` must be TRUE or FALSE.")
  }
  progress <- isTRUE(progress)

  if (is.null(chunk_size)) {
    chunk_size <- cores
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size <= 0) {
    stop("`bootstrap_chunk_size` must be a positive integer or NULL.")
  }
  chunk_size <- as.integer(chunk_size)

  selected <- list(
    VIRF = calc_virf,
    CIRF = calc_cirf,
    KIRF = calc_kirf,
    SIRF = calc_sirf,
    WIRF = calc_wirf
  )

  paths <- list()
  for (name in names(selected)[unlist(selected, use.names = FALSE)]) {
    mat <- main_out[[paste0(name, "_mean")]]
    if (!is.null(mat)) {
      paths[[name]] <- array(
        NA_real_,
        dim = c(nrow(mat), ncol(mat), length(candidate_indices)),
        dimnames = list(NULL, colnames(mat), NULL)
      )
    }
  }

  errors <- rep(NA_character_, R)
  used_indices <- integer(0L)
  local_pos <- 0L
  B_zero <- matrix(0, K, K)

  worker_fun <- function(b) {
    B_b <- if (asym) bekk_bootstrap$B[, , b] else B_zero

    boot_res <- tryCatch(
      compute_irf_core_cpp(
        H_0 = H_0,
        shock = shock,
        xi = xi,
        C = bekk_bootstrap$C0[, , b],
        A = bekk_bootstrap$A[, , b],
        G = bekk_bootstrap$G[, , b],
        B = B_b,
        signs = signs,
        psi_kurt = psi_kurt,
        psi_skew = psi_skew,
        root_type = root_type_id,
        asym = asym,
        simsamp = simsamp,
        n_ahead = n_ahead,
        seed = seed,
        calc_virf = calc_virf,
        calc_cirf = calc_cirf,
        calc_kirf = calc_kirf,
        calc_sirf = calc_sirf,
        calc_wirf = calc_wirf
      ),
      error = function(e) e
    )

    if (inherits(boot_res, "error")) {
      return(list(index = b, result = NULL, error = conditionMessage(boot_res)))
    }

    list(index = b, result = boot_res, error = NA_character_)
  }
  environment(worker_fun) <- list2env(
    list(
      asym = asym,
      bekk_bootstrap = bekk_bootstrap,
      B_zero = B_zero,
      H_0 = H_0,
      shock = shock,
      xi = xi,
      signs = signs,
      psi_kurt = psi_kurt,
      psi_skew = psi_skew,
      root_type_id = root_type_id,
      simsamp = simsamp,
      n_ahead = n_ahead,
      seed = seed,
      calc_virf = calc_virf,
      calc_cirf = calc_cirf,
      calc_kirf = calc_kirf,
      calc_sirf = calc_sirf,
      calc_wirf = calc_wirf,
      compute_irf_core_cpp = compute_irf_core_cpp
    ),
    parent = baseenv()
  )

  chunks <- split(candidate_indices, ceiling(seq_along(candidate_indices) / chunk_size))

  pb <- NULL
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = length(candidate_indices), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  cl <- NULL
  if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)

    lib_paths <- .libPaths()
    parallel::clusterExport(cl, varlist = "lib_paths", envir = environment())
    worker_ready <- tryCatch(
      parallel::clusterEvalQ(cl, {
        .libPaths(lib_paths)
        ns <- loadNamespace("bekkIRF")
        exists("compute_irf_core_cpp", envir = ns, inherits = FALSE)
      }),
      error = function(e) e
    )

    if (inherits(worker_ready, "error") || !all(unlist(worker_ready, use.names = FALSE))) {
      parallel::stopCluster(cl)
      cl <- NULL
      cores <- 1L
      warning(
        "Parallel bootstrap IRF workers could not load `bekkIRF`; ",
        "falling back to sequential computation.",
        call. = FALSE
      )
    } else {
      parallel::clusterExport(
        cl,
        varlist = c("worker_fun"),
        envir = environment()
      )
    }
  }

  run_chunk <- function(chunk) {
    if (cores == 1L) {
      lapply(chunk, worker_fun)
    } else {
      parallel::parLapplyLB(cl, chunk, worker_fun)
    }
  }

  done <- 0L
  for (chunk in chunks) {
    chunk_res <- run_chunk(chunk)

    for (res in chunk_res) {
      b <- res$index

      if (!is.na(res$error)) {
        errors[b] <- res$error
        next
      }

      finite_result <- TRUE
      boot_res <- res$result

      for (name in names(paths)) {
        mat <- boot_res[[paste0(name, "_mean")]]
        if (is.null(mat) || any(!is.finite(mat))) {
          finite_result <- FALSE
        }
      }

      if (!finite_result) {
        errors[b] <- "Non-finite bootstrap IRF result."
        next
      }

      local_pos <- local_pos + 1L
      used_indices <- c(used_indices, b)
      for (name in names(paths)) {
        paths[[name]][, , local_pos] <- boot_res[[paste0(name, "_mean")]]
      }
    }

    done <- done + length(chunk)
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, done)
    }
  }

  if (length(used_indices) == 0L) {
    stop("All bootstrap IRF computations failed.")
  }

  for (name in names(paths)) {
    paths[[name]] <- paths[[name]][, , seq_len(length(used_indices)), drop = FALSE]
  }

  list(
    paths = paths,
    requested_replications = R,
    candidate_replications = length(candidate_indices),
    used_replications = length(used_indices),
    used_indices = used_indices,
    skipped_indices = setdiff(seq_len(R), used_indices),
    errors = errors,
    cores = cores,
    chunk_size = chunk_size
  )
}

compute_irf_bootstrap_ci <- function(paths, ci_probs) {
  ci <- list()

  for (name in names(paths)) {
    arr <- paths[[name]]
    lower <- apply(arr, c(1, 2), stats::quantile, probs = ci_probs[1L], na.rm = TRUE)
    upper <- apply(arr, c(1, 2), stats::quantile, probs = ci_probs[2L], na.rm = TRUE)
    ci[[name]] <- list(
      lower = lower,
      upper = upper,
      probs = ci_probs
    )
  }

  ci
}
