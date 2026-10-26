test_that("bekk_bootstrap supports BEKKs bekk, dbekk, and sbekk model outputs", {
  skip_if_not_installed("BEKKs")

  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- as.matrix(gold_msci_returns[1:80, ])
  model_data <- scale(model_data, center = TRUE, scale = FALSE)
  K <- ncol(model_data)

  for (model_type in c("bekk", "dbekk", "sbekk")) {
    spec <- BEKKs::bekk_spec(model = list(type = model_type, asymmetric = FALSE))
    fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 1, crit = 1e-2)

    expect_warning(boot <- bekk_bootstrap(
      fit,
      bekk_spec_model = spec,
      bootsamp = 2,
      cores = 1,
      max_iter = 1,
      crit = 1e-2,
      progress = FALSE,
      xi_outlier_threshold = NULL
    ), "bootstrap fits were excluded")

    expect_s3_class(boot, "bekkBootstrap")
    expect_equal(boot$settings$model_type, model_type)
    expect_equal(dim(boot$C0), c(K, K, 2))
    expect_equal(dim(boot$A), c(K, K, 2))
    expect_equal(dim(boot$G), c(K, K, 2))
    expect_null(boot$B)
    expect_equal(length(boot$converged), 2)
    expect_false(any(boot$converged))
    expect_true(all(boot$fit_diagnostics$iter == 1))
    expect_true(all(grepl("max_iter", boot$error, fixed = TRUE)))
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

  expect_warning(boot <- bekk_bootstrap(
    fit,
    bekk_spec_model = spec,
    bootsamp = 1,
    cores = 1,
    max_iter = 1,
    crit = 1e-2,
    progress = FALSE,
    xi_outlier_threshold = NULL
  ), "bootstrap fits were excluded")

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

  expect_warning(boot <- bekk_bootstrap(
    fit,
    bekk_spec_model = spec,
    bootsamp = 1,
    cores = 1,
    max_iter = 1,
    crit = 1e-2,
    progress = FALSE,
    xi_outlier_threshold = NULL
  ), "bootstrap fits were excluded")

  expect_s3_class(boot, "bekkBootstrap")
  expect_true(boot$settings$asymmetric)
  expect_equal(boot$settings$model_type, "sbekk")
  expect_equal(dim(boot$A), c(K, K, 1))
  expect_equal(dim(boot$B), c(K, K, 1))
  expect_equal(dim(boot$G), c(K, K, 1))
  expect_false(any(boot$converged))
  expect_true(all(boot$fit_diagnostics$iter == 1))
  expect_true(all(is.finite(boot$B)))
})

test_that("bootstrap acceptance requires valid finite fits below the iteration limit", {
  params <- list(C0 = diag(2), A = diag(0.1, 2),
                 G = diag(0.8, 2), B = matrix(0, 2, 2))
  fit <- list(iter = 4, BEKK_valid = TRUE, log_likelihood = -100)
  expect_true(bekk_bootstrap_fit_diagnostics(fit, params, 50)$accepted)
  fit$BEKK_valid <- FALSE
  expect_false(bekk_bootstrap_fit_diagnostics(fit, params, 50)$accepted)
  fit$BEKK_valid <- TRUE
  fit$iter <- 50
  expect_match(bekk_bootstrap_fit_diagnostics(fit, params, 50)$error, "max_iter")
  fit$iter <- NULL
  expect_false(bekk_bootstrap_fit_diagnostics(fit, params, 50)$accepted)
  fit$iter <- 4
  params$A[1, 1] <- Inf
  expect_false(bekk_bootstrap_fit_diagnostics(fit, params, 50)$accepted)
})

test_that("bootstrap accepts actual fits that finish below the iteration limit", {
  skip_if_not_installed("BEKKs")
  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- scale(as.matrix(gold_msci_returns[1:80, ]),
                      center = TRUE, scale = FALSE)
  spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = FALSE))
  fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 100, crit = 1e-5)
  boot <- bekk_bootstrap(fit, spec, bootsamp = 2, cores = 1,
                         max_iter = 100, crit = 1e-5, progress = FALSE,
                         xi_outlier_threshold = NULL)
  expect_true(all(boot$converged))
  expect_true(all(boot$fit_diagnostics$valid))
  expect_true(all(boot$fit_diagnostics$iter < 100))
  expect_true(all(is.na(boot$error)))
})
