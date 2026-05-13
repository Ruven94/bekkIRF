test_that("compute_irf supports spectral and cholesky roots", {
  set.seed(1)

  K <- 2
  N <- 30
  data <- matrix(rnorm(N * K), nrow = N, ncol = K)
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  out_spec <- compute_irf(
    bekk_model,
    root_type = "spectral",
    shock = c(1, 0),
    time = 10,
    simsamp = 5,
    n.ahead = 3
  )

  out_chol <- compute_irf(
    bekk_model,
    root_type = "cholesky",
    shock_type = "empirical",
    time = 10,
    simsamp = 5,
    n.ahead = 3
  )

  expect_equal(dim(out_spec$VIRF_mean), c(3, 3))
  expect_equal(dim(out_spec$CIRF_mean), c(3, 1))
  expect_equal(dim(out_chol$VIRF_mean), c(3, 3))
  expect_equal(dim(out_chol$CIRF_mean), c(3, 1))
  expect_s3_class(out_spec, "bekkIRF")
  expect_s3_class(out_chol, "bekkIRF")
  expect_null(out_spec$bootstrap_irf)
  expect_null(out_chol$bootstrap_irf)
})

test_that("compute_irf seed is reproducible", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  out_1 <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 5,
    n.ahead = 3,
    seed = 42
  )

  out_2 <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 5,
    n.ahead = 3,
    seed = 42
  )

  expect_equal(out_1$VIRF_mean, out_2$VIRF_mean)
  expect_equal(out_1$CIRF_mean, out_2$CIRF_mean)
})

test_that("compute_irf computes bootstrap IRF paths and confidence intervals", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  bekk_boot <- list(
    C0 = array(rep(diag(0.1, K), 3), dim = c(K, K, 3)),
    A = array(c(
      diag(0.08, K),
      diag(0.10, K),
      diag(0.12, K)
    ), dim = c(K, K, 3)),
    G = array(rep(diag(0.8, K), 3), dim = c(K, K, 3)),
    B = NULL,
    converged = c(TRUE, TRUE, FALSE)
  )
  class(bekk_boot) <- "bekkBootstrap"

  out <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 5,
    n.ahead = 3,
    seed = 42,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_wirf = FALSE,
    bekk_bootstrap = bekk_boot,
    ci_level = 0.8,
    bootstrap_cores = 1,
    bootstrap_progress = FALSE
  )

  expect_s3_class(out, "bekkIRF")
  expect_equal(out$bootstrap_info$ci_probs, c(0.1, 0.9))
  expect_equal(out$bootstrap_info$requested_replications, 3)
  expect_equal(out$bootstrap_info$used_replications, 2)
  expect_equal(out$bootstrap_info$used_indices, c(1L, 2L))
  expect_equal(out$bootstrap_info$cores, 1L)
  expect_equal(out$bootstrap_info$chunk_size, 1L)
  expect_equal(dim(out$bootstrap_irf$VIRF), c(3, 3, 2))
  expect_equal(dim(out$bootstrap_irf$CIRF), c(3, 1, 2))
  expect_equal(dim(out$ci$VIRF$lower), c(3, 3))
  expect_equal(dim(out$ci$CIRF$upper), c(3, 1))
})

test_that("compute_irf accepts explicit bootstrap quantile probabilities", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  bekk_boot <- list(
    C0 = array(rep(diag(0.1, K), 2), dim = c(K, K, 2)),
    A = array(rep(diag(0.1, K), 2), dim = c(K, K, 2)),
    G = array(rep(diag(0.8, K), 2), dim = c(K, K, 2)),
    B = NULL,
    converged = c(TRUE, TRUE)
  )
  class(bekk_boot) <- "bekkBootstrap"

  out <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 3,
    n.ahead = 2,
    calc_virf = TRUE,
    calc_cirf = FALSE,
    bekk_bootstrap = bekk_boot,
    ci_level = c(0.05, 0.95),
    bootstrap_cores = 1,
    bootstrap_progress = FALSE
  )

  expect_equal(out$bootstrap_info$ci_probs, c(0.05, 0.95))
  expect_null(out$bootstrap_irf$CIRF)
  expect_equal(dim(out$bootstrap_irf$VIRF), c(2, 3, 2))
})

test_that("compute_irf bootstrap paths are invariant to parallel execution", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  bekk_boot <- list(
    C0 = array(rep(diag(0.1, K), 3), dim = c(K, K, 3)),
    A = array(c(
      diag(0.08, K),
      diag(0.10, K),
      diag(0.12, K)
    ), dim = c(K, K, 3)),
    G = array(rep(diag(0.8, K), 3), dim = c(K, K, 3)),
    B = NULL,
    converged = c(TRUE, TRUE, TRUE)
  )
  class(bekk_boot) <- "bekkBootstrap"

  out_seq <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 4,
    n.ahead = 3,
    seed = 42,
    bekk_bootstrap = bekk_boot,
    ci_level = 0.8,
    bootstrap_cores = 1,
    bootstrap_progress = FALSE
  )

  out_par <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 4,
    n.ahead = 3,
    seed = 42,
    bekk_bootstrap = bekk_boot,
    ci_level = 0.8,
    bootstrap_cores = 2,
    bootstrap_progress = FALSE,
    bootstrap_chunk_size = 1
  )

  expect_equal(out_par$bootstrap_info$cores, 2L)
  expect_equal(out_par$bootstrap_info$chunk_size, 1L)
  expect_equal(out_par$bootstrap_irf, out_seq$bootstrap_irf)
  expect_equal(out_par$ci, out_seq$ci)
})

test_that("compute_irf validates bootstrap parallel control arguments", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  bekk_boot <- list(
    C0 = array(rep(diag(0.1, K), 2), dim = c(K, K, 2)),
    A = array(rep(diag(0.1, K), 2), dim = c(K, K, 2)),
    G = array(rep(diag(0.8, K), 2), dim = c(K, K, 2)),
    B = NULL,
    converged = c(TRUE, TRUE)
  )
  class(bekk_boot) <- "bekkBootstrap"

  expect_error(
    compute_irf(
      bekk_model,
      shock = c(1, 0),
      time = 10,
      simsamp = 3,
      n.ahead = 2,
      bekk_bootstrap = bekk_boot,
      bootstrap_cores = 0,
      bootstrap_progress = FALSE
    ),
    "`bootstrap_cores` must be a positive integer.",
    fixed = TRUE
  )

  expect_error(
    compute_irf(
      bekk_model,
      shock = c(1, 0),
      time = 10,
      simsamp = 3,
      n.ahead = 2,
      bekk_bootstrap = bekk_boot,
      bootstrap_cores = 1,
      bootstrap_progress = NA
    ),
    "`bootstrap_progress` must be TRUE or FALSE.",
    fixed = TRUE
  )

  expect_error(
    compute_irf(
      bekk_model,
      shock = c(1, 0),
      time = 10,
      simsamp = 3,
      n.ahead = 2,
      bekk_bootstrap = bekk_boot,
      bootstrap_cores = 1,
      bootstrap_progress = FALSE,
      bootstrap_chunk_size = 0
    ),
    "`bootstrap_chunk_size` must be a positive integer or NULL.",
    fixed = TRUE
  )
})
