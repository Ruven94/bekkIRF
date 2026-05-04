test_that("compute_theoretical_virf returns vech output for structural shocks", {
  K <- 2
  N <- 10
  H <- diag(c(2, 3))

  bekk_model <- list(
    H_t = matrix(rep(as.vector(H), times = N), nrow = N, byrow = TRUE),
    data = matrix(0, nrow = N, ncol = K),
    A = diag(c(0.2, 0.3)),
    G = diag(c(0.6, 0.7)),
    asymmetric = FALSE
  )

  out <- compute_theoretical_virf(
    bekk_model,
    shock = c(1, 0),
    time = 3,
    n.ahead = 4
  )

  expect_equal(dim(out), c(4, 3))
  expect_equal(colnames(out), c("H11", "H21", "H22"))
})

test_that("compute_theoretical_virf first horizon follows HH2006 closed form", {
  K <- 2
  N <- 10
  H <- matrix(c(2, 0.4, 0.4, 3), 2, 2)
  A <- matrix(c(0.2, 0.1, 0.05, 0.3), 2, 2)
  G <- matrix(c(0.6, 0.05, 0.02, 0.7), 2, 2)
  shock <- c(1.2, -0.7)

  bekk_model <- list(
    H_t = matrix(rep(as.vector(H), times = N), nrow = N, byrow = TRUE),
    data = matrix(0, nrow = N, ncol = K),
    A = A,
    G = G,
    asymmetric = FALSE
  )

  out <- compute_theoretical_virf(
    bekk_model,
    root_type = "spectral",
    shock = shock,
    time = 4,
    n.ahead = 2,
    format = "array"
  )

  Q <- matroot(H, type = "spec")
  expected_first <- t(A) %*% (tcrossprod(Q %*% shock) - H) %*% A
  expected_second <- t(A) %*% expected_first %*% A +
    t(G) %*% expected_first %*% G

  expect_equal(out[, , 1], expected_first, tolerance = 1e-12)
  expect_equal(out[, , 2], expected_second, tolerance = 1e-12)
})

test_that("compute_theoretical_virf empirical shock matches equivalent structural shock", {
  K <- 2
  N <- 12
  H <- diag(c(1.5, 2.5))
  data <- matrix(seq_len(N * K) / 10, nrow = N, ncol = K)

  bekk_model <- list(
    H_t = matrix(rep(as.vector(H), times = N), nrow = N, byrow = TRUE),
    data = data,
    A = diag(c(0.2, 0.25)),
    G = diag(c(0.6, 0.65)),
    asymmetric = FALSE
  )

  xi <- compute_xi(bekk_model$H_t, data, root_type = "spec")

  out_empirical <- compute_theoretical_virf(
    bekk_model,
    shock_type = "empirical",
    time = 5,
    n.ahead = 3
  )

  out_structural <- compute_theoretical_virf(
    bekk_model,
    shock_type = "structural",
    shock = xi[5, ],
    time = 5,
    n.ahead = 3
  )

  expect_equal(out_empirical, out_structural, tolerance = 1e-12)
})

test_that("compute_theoretical_virf supports root aliases", {
  K <- 2
  N <- 10

  bekk_model <- list(
    H_t = matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE),
    data = matrix(0, nrow = N, ncol = K),
    A = diag(0.2, K),
    G = diag(0.7, K),
    asymmetric = FALSE
  )

  out_spectral <- compute_theoretical_virf(
    bekk_model,
    root_type = "spectral",
    shock = c(1, 0),
    time = 2,
    n.ahead = 3
  )

  out_spec <- compute_theoretical_virf(
    bekk_model,
    root_type = "spec",
    shock = c(1, 0),
    time = 2,
    n.ahead = 3
  )

  out_cholesky <- compute_theoretical_virf(
    bekk_model,
    root_type = "cholesky",
    shock = c(1, 0),
    time = 2,
    n.ahead = 3
  )

  out_chol <- compute_theoretical_virf(
    bekk_model,
    root_type = "chol",
    shock = c(1, 0),
    time = 2,
    n.ahead = 3
  )

  expect_equal(out_spectral, out_spec)
  expect_equal(out_cholesky, out_chol)
})

test_that("compute_theoretical_virf rejects asymmetric models", {
  K <- 2
  N <- 10

  bekk_model <- list(
    H_t = matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE),
    data = matrix(0, nrow = N, ncol = K),
    A = diag(0.2, K),
    G = diag(0.7, K),
    asymmetric = TRUE
  )

  expect_error(
    compute_theoretical_virf(bekk_model, shock = c(1, 0), time = 2),
    "only for symmetric BEKK models",
    fixed = TRUE
  )
})

test_that("compute_theoretical_virf validates shock and time inputs", {
  K <- 2
  N <- 10

  bekk_model <- list(
    H_t = matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE),
    data = matrix(0, nrow = N, ncol = K),
    A = diag(0.2, K),
    G = diag(0.7, K),
    asymmetric = FALSE
  )

  expect_error(
    compute_theoretical_virf(bekk_model, time = 2),
    "`shock` must be supplied",
    fixed = TRUE
  )

  expect_error(
    compute_theoretical_virf(bekk_model, shock = c(1, 0), time = 0),
    "`time` must be between 1",
    fixed = TRUE
  )

  expect_error(
    compute_theoretical_virf(bekk_model, shock = c(1, 0, 0), time = 2),
    "`shock` must have length",
    fixed = TRUE
  )
})
