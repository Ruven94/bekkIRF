test_that("spectral root reconstructs a positive definite matrix", {
  M <- matrix(c(4, 1,
                1, 3), 2, 2)

  R <- matroot(M, type = "spec")

  expect_true(is.matrix(R))
  expect_equal(R %*% R, M, tolerance = 1e-8)
})

test_that("spectral root reconstructs a positive semi-definite matrix", {
  M <- matrix(c(1, 2,
                2, 4), 2, 2)

  R <- matroot(M, type = "spec")

  expect_true(is.matrix(R))
  expect_equal(R %*% R, M, tolerance = 1e-8)
})

test_that("cholesky root reconstructs a positive definite matrix", {
  M <- matrix(c(4, 1,
                1, 3), 2, 2)

  R <- matroot(M, type = "chol")

  expect_true(is.matrix(R))
  expect_equal(R %*% t(R), M, tolerance = 1e-8)
})

test_that("cholesky root fails for positive semi-definite but not positive definite matrix", {
  M <- matrix(c(1, 2,
                2, 4), 2, 2)

  expect_error(
    matroot(M, type = "chol"),
    "positive definite|leading minor"
  )
})

test_that("matroot fails for non-matrix input", {
  x <- c(1, 2, 3)

  expect_error(
    matroot(x),
    "matrix"
  )
})

test_that("matroot fails for non-square matrix", {
  M <- matrix(1:6, nrow = 2, ncol = 3)

  expect_error(
    matroot(M),
    "square"
  )
})

test_that("matroot fails for non-symmetric matrix", {
  M <- matrix(c(1, 2,
                3, 4), 2, 2)

  expect_error(
    matroot(M),
    "symmetric"
  )
})

test_that("spectral root fails for matrix with truly negative eigenvalue", {
  M <- matrix(c(1, 2,
                2, 1), 2, 2)

  expect_error(
    matroot(M, type = "spec"),
    "positive semi-definite"
  )
})

test_that("type argument is matched correctly", {
  M <- matrix(c(4, 1,
                1, 3), 2, 2)

  expect_error(
    matroot(M, type = "banana"),
    "arg should be one of|should be one of"
  )
})

test_that("spectral root remains symmetric", {
  M <- matrix(c(4, 1,
                1, 3), 2, 2)

  R <- matroot(M, type = "spec")

  expect_equal(R, t(R), tolerance = 1e-8)
})
