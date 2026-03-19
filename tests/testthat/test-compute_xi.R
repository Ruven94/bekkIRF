test_that("compute_xi returns matrix with correct dimensions", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H_t <- rbind(
    c(1, 0, 0, 1),
    c(1, 0, 0, 1)
  )

  xi <- compute_xi(H_t, data, root_type = "spec")

  expect_true(is.matrix(xi))
  expect_equal(dim(xi), dim(data))
})

test_that("compute_xi reproduces data for identity covariance matrices", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H_t <- rbind(
    c(1, 0, 0, 1),
    c(1, 0, 0, 1)
  )

  xi <- compute_xi(H_t, data, root_type = "spec")

  expect_equal(xi, data, tolerance = 1e-8)
})

test_that("compute_xi works with spectral root for positive definite matrices", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H1 <- matrix(c(2, 0.5,
                 0.5, 1.5), 2, 2)

  H_t <- rbind(
    as.vector(H1),
    as.vector(H1)
  )

  xi <- compute_xi(H_t, data, root_type = "spec")

  expect_true(is.matrix(xi))
  expect_equal(dim(xi), dim(data))

  Q <- matroot(H1, type = "spec")
  expect_equal(drop(Q %*% xi[1, ]), data[1, ], tolerance = 1e-8)
  expect_equal(drop(Q %*% xi[2, ]), data[2, ], tolerance = 1e-8)
})

test_that("compute_xi works with cholesky root for positive definite matrices", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H1 <- matrix(c(2, 0.5,
                 0.5, 1.5), 2, 2)

  H_t <- rbind(
    as.vector(H1),
    as.vector(H1)
  )

  xi <- compute_xi(H_t, data, root_type = "chol")

  expect_true(is.matrix(xi))
  expect_equal(dim(xi), dim(data))

  Q <- matroot(H1, type = "chol")
  expect_equal(drop(Q %*% xi[1, ]), data[1, ], tolerance = 1e-8)
  expect_equal(drop(Q %*% xi[2, ]), data[2, ], tolerance = 1e-8)
})

test_that("compute_xi fails if H_t is not a matrix", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H_t <- c(1, 0, 0, 1)

  expect_error(
    compute_xi(H_t, data),
    "`H_t` must be a matrix.",
    fixed = TRUE
  )
})

test_that("compute_xi fails if data is not a matrix", {
  data <- c(1, 2, 3, 4)

  H_t <- rbind(
    c(1, 0, 0, 1),
    c(1, 0, 0, 1)
  )

  expect_error(
    compute_xi(H_t, data),
    "`data` must be a matrix.",
    fixed = TRUE
  )
})

test_that("compute_xi fails if H_t and data have different numbers of rows", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H_t <- matrix(c(1, 0, 0, 1), nrow = 1)

  expect_error(
    compute_xi(H_t, data),
    "same number of rows"
  )
})

test_that("compute_xi fails if H_t has wrong number of columns", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H_t <- matrix(1, nrow = 2, ncol = 3)

  expect_error(
    compute_xi(H_t, data),
    "columns"
  )
})

test_that("compute_xi fails for cholesky root if covariance matrix is not positive definite", {
  data <- matrix(c(1, 2,
                   3, 4), ncol = 2, byrow = TRUE)

  H1 <- matrix(c(1, 2,
                 2, 4), 2, 2)

  H_t <- rbind(
    as.vector(H1),
    as.vector(H1)
  )

  expect_error(
    compute_xi(H_t, data, root_type = "chol")
  )
})

test_that("compute_xi matches manual computation for identity covariance", {
  data <- matrix(c(5, -1,
                   2,  7), ncol = 2, byrow = TRUE)

  H_t <- rbind(
    c(1, 0, 0, 1),
    c(1, 0, 0, 1)
  )

  xi <- compute_xi(H_t, data)

  expect_equal(xi[1, ], data[1, ], tolerance = 1e-8)
  expect_equal(xi[2, ], data[2, ], tolerance = 1e-8)
})

