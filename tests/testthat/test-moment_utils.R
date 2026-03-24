test_that("compute_psi_kurtosis returns matrix with correct dimensions", {
  xi <- matrix(rnorm(20), ncol = 2)

  Psi <- compute_psi_kurtosis(xi)

  expect_true(is.matrix(Psi))
  expect_equal(dim(Psi), c(4, 4))
})

test_that("compute_psi_skewness returns matrix with correct dimensions", {
  xi <- matrix(rnorm(20), ncol = 2)

  Psi <- compute_psi_skewness(xi)

  expect_true(is.matrix(Psi))
  expect_equal(dim(Psi), c(4, 2))
})

test_that("compute_psi_kurtosis matches manual calculation for simple input", {
  xi <- matrix(c(1, 2,
                 3, 4), ncol = 2, byrow = TRUE)

  A1 <- tcrossprod(xi[1, ])
  A2 <- tcrossprod(xi[2, ])
  Psi_manual <- (kronecker(A1, A1) + kronecker(A2, A2)) / 2

  Psi <- compute_psi_kurtosis(xi)

  expect_equal(Psi, Psi_manual, tolerance = 1e-8)
})

test_that("compute_psi_skewness matches manual calculation for simple input", {
  xi <- matrix(c(1, 2,
                 3, 4), ncol = 2, byrow = TRUE)

  A1 <- tcrossprod(xi[1, ])
  A2 <- tcrossprod(xi[2, ])
  B1 <- xi[1, ]
  B2 <- xi[2, ]

  Psi_manual <- (kronecker(A1, B1) + kronecker(A2, B2)) / 2

  Psi <- compute_psi_skewness(xi)

  expect_equal(Psi, Psi_manual, tolerance = 1e-8)
})

test_that("compute_psi_kurtosis works for one-row input", {
  xi <- matrix(c(1, -2), nrow = 1)

  A <- tcrossprod(xi[1, ])
  Psi_manual <- kronecker(A, A)

  Psi <- compute_psi_kurtosis(xi)

  expect_equal(Psi, Psi_manual, tolerance = 1e-8)
})

test_that("compute_psi_skewness works for one-row input", {
  xi <- matrix(c(1, -2), nrow = 1)

  A <- tcrossprod(xi[1, ])
  B <- xi[1, ]
  Psi_manual <- kronecker(A, B)

  Psi <- compute_psi_skewness(xi)

  expect_equal(Psi, Psi_manual, tolerance = 1e-8)
})

test_that("compute_psi_kurtosis fails if xi is not a matrix", {
  expect_error(
    compute_psi_kurtosis(1:10),
    "`xi` must be a matrix.",
    fixed = TRUE
  )
})

test_that("compute_psi_skewness fails if xi is not a matrix", {
  expect_error(
    compute_psi_skewness(1:10),
    "`xi` must be a matrix.",
    fixed = TRUE
  )
})

test_that("compute_psi_kurtosis fails if xi is not numeric", {
  xi <- matrix(c("a", "b", "c", "d"), ncol = 2)

  expect_error(
    compute_psi_kurtosis(xi),
    "`xi` must be numeric.",
    fixed = TRUE
  )
})

test_that("compute_psi_skewness fails if xi is not numeric", {
  xi <- matrix(c("a", "b", "c", "d"), ncol = 2)

  expect_error(
    compute_psi_skewness(xi),
    "`xi` must be numeric.",
    fixed = TRUE
  )
})

test_that("compute_psi_kurtosis fails for empty input", {
  xi <- matrix(numeric(0), ncol = 2)

  expect_error(
    compute_psi_kurtosis(xi),
    "`xi` must have at least one row.",
    fixed = TRUE
  )
})

test_that("compute_psi_skewness fails for empty input", {
  xi <- matrix(numeric(0), ncol = 2)

  expect_error(
    compute_psi_skewness(xi),
    "`xi` must have at least one row.",
    fixed = TRUE
  )
})

test_that("compute_psi_kurtosis returns symmetric matrix", {
  xi <- matrix(rnorm(20), ncol = 2)

  Psi <- compute_psi_kurtosis(xi)

  expect_equal(Psi, t(Psi), tolerance = 1e-8)
})


