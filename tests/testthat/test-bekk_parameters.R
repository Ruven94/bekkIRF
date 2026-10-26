test_that("scalar BEKK responses use linear BEKKs coefficients", {
  # With zero future innovations, the recursion is deterministic. This tests
  # the actual response, including persistence, rather than just matrix shape.
  K <- 2L
  N <- 4L
  fit <- list(
    H_t = matrix(rep(as.vector(diag(K)), N), nrow = N, byrow = TRUE),
    data = matrix(0, N, K),
    C0 = diag(0.2, K), a = 0.09, g = 0.7, asymmetric = FALSE
  )
  shock <- c(-2, -1)
  for (asymmetric in c(FALSE, TRUE)) {
    fit$asymmetric <- asymmetric
    fit$b <- if (asymmetric) 0.04 else NULL
    irf <- compute_irf(fit, shock = shock, time = 1, simsamp = 3,
                       n.ahead = 4, calc_wirf = TRUE)
    base <- shocked <- diag(K)
    expected <- matrix(0, 4, 3)
    expected_cor <- matrix(0, 4, 1)
    expected_weights <- matrix(0, 4, 2)
    for (h in seq_len(4)) {
      base <- tcrossprod(fit$C0) + fit$g * base
      shocked <- tcrossprod(fit$C0) + fit$g * shocked
      if (h == 1) {
        shocked <- shocked + (fit$a + if (asymmetric) fit$b else 0) *
          tcrossprod(shock)
      }
      expected[h, ] <- (shocked - base)[lower.tri(base, diag = TRUE)]
      expected_cor[h, ] <- cov2cor(shocked)[2, 1] - cov2cor(base)[2, 1]
      ws <- solve(shocked, rep(1, K))
      wb <- solve(base, rep(1, K))
      expected_weights[h, ] <- ws / sum(ws) - wb / sum(wb)
    }
    expect_equal(unname(irf$VIRF_mean), expected, tolerance = 1e-12)
    expect_equal(unname(irf$CIRF_mean), expected_cor, tolerance = 1e-12)
    expect_equal(unname(irf$WIRF_mean), expected_weights, tolerance = 1e-12)
  }
})

test_that("scalar dynamic coefficients must be finite and non-negative", {
  for (field in c("a", "b", "g")) {
    for (value in c(NA_real_, Inf, -Inf, -0.1)) {
      model <- setNames(list(value), field)
      expect_error(bekk_parameter_matrix(model, toupper(field), field, 2),
                   "finite|non-negative")
    }
  }
})
