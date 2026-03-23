test_that("compute_eta returns e_t if joint sign condition is satisfied", {
  e_t <- c(-2, -1)
  signs <- c(-1, -1)

  expect_equal(compute_eta(e_t, signs), e_t)
})

test_that("compute_eta returns zero vector if joint sign condition is not satisfied", {
  e_t <- c(-2, 1)
  signs <- c(-1, -1)

  expect_equal(compute_eta(e_t, signs), c(0, 0))
})

test_that("compute_eta works for mixed sign pattern", {
  e_t <- c(-2, 3)
  signs <- c(-1, 1)

  expect_equal(compute_eta(e_t, signs), e_t)
})

test_that("compute_eta returns zero for partially matching mixed sign pattern", {
  e_t <- c(-2, -3)
  signs <- c(-1, 1)

  expect_equal(compute_eta(e_t, signs), c(0, 0))
})

test_that("compute_eta fails if e_t is not numeric", {
  expect_error(
    compute_eta(e_t = c("a", "b"), signs = c(-1, -1)),
    "`e_t` must be numeric.",
    fixed = TRUE
  )
})

test_that("compute_eta fails if signs is not numeric", {
  expect_error(
    compute_eta(e_t = c(-1, -2), signs = c("a", "b")),
    "`signs` must be numeric.",
    fixed = TRUE
  )
})

test_that("compute_eta fails if signs has wrong length", {
  expect_error(
    compute_eta(e_t = c(-1, -2), signs = c(-1)),
    "`signs` must have the same length as `e_t`.",
    fixed = TRUE
  )
})

test_that("compute_eta fails if signs contains values other than -1 and 1", {
  expect_error(
    compute_eta(e_t = c(-1, -2), signs = c(-1, 0)),
    "`signs` must contain only -1 and 1.",
    fixed = TRUE
  )
})

test_that("compute_eta handles all four two-dimensional sign combinations", {
  signs_mat <- matrix(c(
    1,  1,
    -1,  1,
    1, -1,
    -1, -1
  ), 4, 2, byrow = TRUE)

  e_list <- list(
    c( 2,  3),
    c(-2,  3),
    c( 2, -3),
    c(-2, -3)
  )

  for (i in 1:4) {
    expect_equal(compute_eta(e_list[[i]], signs_mat[i, ]), e_list[[i]])
  }
})

test_that("compute_eta is consistent with asymmetric BEKK recursion", {
  # skip_if_not_installed("BEKKs")

  set.seed(123)

  # Small artificial bivariate dataset with all sign combinations represented
  data <- matrix(rnorm(400), ncol = 2)

  # Ensure centering, as typically done in the BEKK workflow
  data <- scale(data, center = TRUE, scale = FALSE)

  signs_mat <- matrix(c(
    1,  1,
    -1,  1,
    1, -1,
    -1, -1
  ), 4, 2, byrow = TRUE)

  K <- ncol(data)

  for (i in 1:4) {
    bekk_asym <- BEKKs::bekk_spec(
      model = list(type = "bekk", asymmetric = TRUE),
      signs = signs_mat[i, ]
    )

    bekk_asym_model <- BEKKs::bekk_fit(bekk_asym, data)

    for (j in 1:4) {
      sign_search <- signs_mat[j, ]

      idx_candidates <- which(apply(sweep(data, 2, sign_search, `*`) > 0, 1, all))

      # We need at least one valid index and also idx + 1 to exist
      idx_candidates <- idx_candidates[idx_candidates < nrow(data)]

      expect_true(length(idx_candidates) > 0)

      idx <- idx_candidates[1]

      e_t <- data[idx, ]
      H_t <- matrix(bekk_asym_model$H_t[idx, ], K, K)
      H_tp1_obj <- matrix(bekk_asym_model$H_t[idx + 1, ], K, K)

      eta_t <- compute_eta(e_t, signs_mat[i, ])

      H_tp1_manual <- bekk_asym_model$C0 %*% t(bekk_asym_model$C0) +
        t(bekk_asym_model$A) %*% (e_t %*% t(e_t)) %*% bekk_asym_model$A +
        t(bekk_asym_model$G) %*% H_t %*% bekk_asym_model$G +
        t(bekk_asym_model$B) %*% (eta_t %*% t(eta_t)) %*% bekk_asym_model$B

      expect_equal(H_tp1_manual, H_tp1_obj, tolerance = 1e-8)
    }
  }
})
