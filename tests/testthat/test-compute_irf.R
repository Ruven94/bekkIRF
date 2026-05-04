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
