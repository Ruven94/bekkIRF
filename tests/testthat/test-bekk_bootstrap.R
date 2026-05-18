test_that("bekk_bootstrap supports BEKKs bekk, dbekk, and sbekk model outputs", {
  skip_if_not_installed("BEKKs")

  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- as.matrix(gold_msci_returns[1:80, ])
  model_data <- scale(model_data, center = TRUE, scale = FALSE)
  K <- ncol(model_data)

  for (model_type in c("bekk", "dbekk", "sbekk")) {
    spec <- BEKKs::bekk_spec(model = list(type = model_type, asymmetric = FALSE))
    fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 1, crit = 1e-2)

    boot <- bekk_bootstrap(
      fit,
      bekk_spec_model = spec,
      bootsamp = 2,
      cores = 1,
      max_iter = 1,
      crit = 1e-2,
      progress = FALSE,
      xi_outlier_threshold = NULL
    )

    expect_s3_class(boot, "bekkBootstrap")
    expect_equal(boot$settings$model_type, model_type)
    expect_equal(dim(boot$C0), c(K, K, 2))
    expect_equal(dim(boot$A), c(K, K, 2))
    expect_equal(dim(boot$G), c(K, K, 2))
    expect_null(boot$B)
    expect_equal(length(boot$converged), 2)
    expect_true(all(boot$converged))
    expect_true(all(is.finite(boot$A)))
    expect_true(all(is.finite(boot$G)))
  }
})

test_that("bekk_bootstrap stores scalar BEKK bootstrap estimates as matrix arrays", {
  skip_if_not_installed("BEKKs")

  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- as.matrix(gold_msci_returns[1:80, ])
  model_data <- scale(model_data, center = TRUE, scale = FALSE)
  K <- ncol(model_data)

  spec <- BEKKs::bekk_spec(model = list(type = "sbekk", asymmetric = FALSE))
  fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 1, crit = 1e-2)

  boot <- bekk_bootstrap(
    fit,
    bekk_spec_model = spec,
    bootsamp = 1,
    cores = 1,
    max_iter = 1,
    crit = 1e-2,
    progress = FALSE,
    xi_outlier_threshold = NULL
  )

  expect_equal(dim(boot$A), c(K, K, 1))
  expect_equal(dim(boot$G), c(K, K, 1))
  expect_true(isTRUE(all.equal(boot$A[, , 1], diag(diag(boot$A[, , 1]), K))))
  expect_true(isTRUE(all.equal(boot$G[, , 1], diag(diag(boot$G[, , 1]), K))))
})

test_that("bekk_bootstrap supports asymmetric scalar BEKKs output", {
  skip_if_not_installed("BEKKs")

  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- as.matrix(gold_msci_returns[1:80, ])
  model_data <- scale(model_data, center = TRUE, scale = FALSE)
  K <- ncol(model_data)

  spec <- BEKKs::bekk_spec(model = list(type = "sbekk", asymmetric = TRUE))
  fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 1, crit = 1e-2)

  boot <- bekk_bootstrap(
    fit,
    bekk_spec_model = spec,
    bootsamp = 1,
    cores = 1,
    max_iter = 1,
    crit = 1e-2,
    progress = FALSE,
    xi_outlier_threshold = NULL
  )

  expect_s3_class(boot, "bekkBootstrap")
  expect_true(boot$settings$asymmetric)
  expect_equal(boot$settings$model_type, "sbekk")
  expect_equal(dim(boot$A), c(K, K, 1))
  expect_equal(dim(boot$B), c(K, K, 1))
  expect_equal(dim(boot$G), c(K, K, 1))
  expect_true(all(boot$converged))
  expect_true(all(is.finite(boot$B)))
})
