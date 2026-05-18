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

test_that("compute_irf supports BEKKs bekk, dbekk, and sbekk model outputs", {
  skip_if_not_installed("BEKKs")

  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- as.matrix(gold_msci_returns[1:80, ])
  model_data <- scale(model_data, center = TRUE, scale = FALSE)

  for (model_type in c("bekk", "dbekk", "sbekk")) {
    spec <- BEKKs::bekk_spec(model = list(type = model_type, asymmetric = FALSE))
    fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 1, crit = 1e-2)

    out <- compute_irf(
      fit,
      shock = c(1, 0),
      time = 10,
      simsamp = 3,
      n.ahead = 2,
      calc_virf = TRUE,
      calc_cirf = TRUE,
      calc_kirf = FALSE,
      calc_sirf = FALSE,
      calc_wirf = FALSE
    )

    expect_s3_class(out, "bekkIRF")
    expect_equal(out$settings$model_type, model_type)
    expect_equal(dim(out$VIRF_mean), c(2, 3))
    expect_equal(dim(out$CIRF_mean), c(2, 1))
    expect_true(all(is.finite(out$VIRF_mean)))
    expect_true(all(is.finite(out$CIRF_mean)))
  }
})

test_that("compute_irf supports asymmetric scalar BEKKs output", {
  skip_if_not_installed("BEKKs")

  data("gold_msci_returns", package = "bekkIRF", envir = environment())
  model_data <- as.matrix(gold_msci_returns[1:80, ])
  model_data <- scale(model_data, center = TRUE, scale = FALSE)

  spec <- BEKKs::bekk_spec(model = list(type = "sbekk", asymmetric = TRUE))
  fit <- BEKKs::bekk_fit(spec, model_data, max_iter = 1, crit = 1e-2)

  out <- compute_irf(
    fit,
    shock = c(1, 0),
    time = 10,
    simsamp = 3,
    n.ahead = 2,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_kirf = TRUE,
    calc_sirf = TRUE,
    calc_wirf = TRUE
  )

  expect_s3_class(out, "bekkIRF")
  expect_true(out$settings$asymmetric)
  expect_equal(out$settings$model_type, "sbekk")
  expect_equal(dim(out$VIRF_mean), c(2, 3))
  expect_equal(dim(out$KIRF_mean), c(2, 3))
  expect_equal(dim(out$SIRF_mean), c(2, 2))
  expect_equal(dim(out$WIRF_mean), c(2, 2))
  expect_true(all(is.finite(out$VIRF_mean)))
})

test_that("compute_irf seed is reproducible", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  colnames(data) <- c("msci.Ret", "gold.Ret")
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
  colnames(data) <- c("msci.Ret", "gold.Ret")
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

test_that("compute_irf computes bootstrap confidence intervals for K greater than 2", {
  K <- 3
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  colnames(data) <- c("MSCI", "Gold", "Oil")
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

  out <- compute_irf(
    bekk_model,
    shock = c(1, 0, 0),
    time = 10,
    simsamp = 4,
    n.ahead = 3,
    seed = 42,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_kirf = TRUE,
    calc_sirf = TRUE,
    calc_wirf = TRUE,
    bekk_bootstrap = bekk_boot,
    ci_level = 0.8,
    bootstrap_cores = 1,
    bootstrap_progress = FALSE
  )

  expect_s3_class(out, "bekkIRF")
  expect_equal(out$settings$series_names, c("MSCI", "Gold", "Oil"))
  expect_equal(out$bootstrap_info$used_replications, 3)
  expect_equal(dim(out$VIRF_mean), c(3, 6))
  expect_equal(dim(out$CIRF_mean), c(3, 3))
  expect_equal(dim(out$SIRF_mean), c(3, 3))
  expect_equal(dim(out$KIRF_mean), c(3, 6))
  expect_equal(dim(out$WIRF_mean), c(3, 3))
  expect_equal(dim(out$bootstrap_irf$VIRF), c(3, 6, 3))
  expect_equal(dim(out$bootstrap_irf$CIRF), c(3, 3, 3))
  expect_equal(dim(out$bootstrap_irf$KIRF), c(3, 6, 3))
  expect_equal(dim(out$ci$VIRF$lower), c(3, 6))
  expect_equal(dim(out$ci$CIRF$upper), c(3, 3))
  expect_true(all(is.finite(out$ci$VIRF$lower)))
  expect_true(all(is.finite(out$ci$CIRF$upper)))
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

test_that("bekkIRF print, summary, and plot methods work", {
  K <- 2
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  colnames(data) <- c("msci.Ret", "gold.Ret")
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  out <- compute_irf(
    bekk_model,
    shock = c(1, 0),
    time = 10,
    simsamp = 3,
    n.ahead = 2,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_kirf = FALSE,
    calc_sirf = FALSE,
    calc_wirf = FALSE
  )

  printed <- capture.output(print(out))
  expect_true(any(grepl("bekkIRF object", printed, fixed = TRUE)))
  expect_true(any(grepl("Bootstrap: none", printed, fixed = TRUE)))

  out_summary <- summary(out)
  expect_s3_class(out_summary, "summary_bekkIRF")
  expect_equal(out_summary$irf$type, c("VIRF", "CIRF"))
  expect_equal(out_summary$settings$n.ahead, 2L)
  expect_true(isFALSE(out_summary$settings$bootstrap))

  summary_printed <- capture.output(print(out_summary))
  expect_true(any(grepl("Summary of bekkIRF object", summary_printed, fixed = TRUE)))

  p_virf <- plot(out, type = "VIRF")
  p_all <- plot(out)
  expect_s3_class(p_virf, "ggplot")
  expect_s3_class(p_all, "ggplot")
  expect_null(p_virf$labels$subtitle)
  expect_equal(unique(as.character(p_virf$data$component)), c(
    "VIRF_{msci.Ret}",
    "VIRF_{msci.Ret, gold.Ret}",
    "VIRF_{gold.Ret}"
  ))

  p_pair <- plot(
    out,
    type = "virf",
    components = c(1, 2),
    title = "Custom title",
    subtitle = "Custom subtitle",
    xlab = "h",
    ylab = "response",
    xlim = c(1, 2),
    ylim = c(-1, 1),
    type_label = "Variance response",
    component_labels = "MSCI-Gold covariance",
    line_color = "red",
    ci_fill = "blue",
    ci_alpha = 0.1,
    line_width = 1.2,
    zero_line = FALSE
  )
  expect_s3_class(p_pair, "ggplot")
  expect_equal(unique(as.character(p_pair$data$component)), "VIRF_{msci.Ret, gold.Ret}")
  expect_equal(p_pair$labels$title, "Custom title")
  expect_equal(p_pair$labels$subtitle, "Custom subtitle")
  expect_equal(p_pair$labels$x, "h")
  expect_equal(p_pair$labels$y, "response")
  expect_equal(unique(as.character(p_pair$data$type_label)), "Variance response")
  expect_equal(unique(as.character(p_pair$data$component_label)), "MSCI-Gold covariance")
  expect_equal(p_pair$coordinates$limits$x, c(1, 2))
  expect_equal(p_pair$coordinates$limits$y, c(-1, 1))

  expect_error(
    plot(out, type = "WIRF"),
    "Requested IRF type is not available: WIRF.",
    fixed = TRUE
  )
})

test_that("bekkIRF plot component labels and filters work for K greater than 2", {
  K <- 3
  N <- 30
  data <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
  colnames(data) <- c("MSCI", "Gold", "Oil")
  H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)

  bekk_model <- list(
    H_t = H_t,
    data = data,
    C0 = diag(0.1, K),
    A = diag(0.1, K),
    G = diag(0.8, K),
    asymmetric = FALSE
  )

  out <- compute_irf(
    bekk_model,
    shock = c(1, 0, 0),
    time = 10,
    simsamp = 3,
    n.ahead = 2,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_kirf = FALSE,
    calc_sirf = TRUE,
    calc_wirf = TRUE
  )

  expect_equal(out$settings$series_names, c("MSCI", "Gold", "Oil"))

  p_virf <- plot(out, type = "virf")
  expect_equal(unique(as.character(p_virf$data$component)), c(
    "VIRF_{MSCI}",
    "VIRF_{MSCI, Gold}",
    "VIRF_{MSCI, Oil}",
    "VIRF_{Gold}",
    "VIRF_{Gold, Oil}",
    "VIRF_{Oil}"
  ))

  p_pair <- plot(out, type = "VIRF", components = c(1, 3))
  expect_equal(unique(as.character(p_pair$data$component)), "VIRF_{MSCI, Oil}")

  p_own <- plot(out, type = "sirf", components = 3)
  expect_equal(unique(as.character(p_own$data$component)), "SIRF_{Oil}")

  p_all_pair <- plot(out, type = "all", components = c(2, 3))
  expect_equal(unique(as.character(p_all_pair$data$type)), c("VIRF", "CIRF"))
  expect_equal(unique(as.character(p_all_pair$data$component)), c(
    "VIRF_{Gold, Oil}",
    "CIRF_{Gold, Oil}"
  ))

  p_custom_labels <- plot(
    out,
    type = "all",
    components = c(2, 3),
    type_label = c(VIRF = "Variance", CIRF = "Correlation"),
    component_labels = c(
      "VIRF_{Gold, Oil}" = "Gold-Oil covariance",
      "CIRF_{Gold, Oil}" = "Gold-Oil correlation"
    )
  )
  expect_equal(unique(as.character(p_custom_labels$data$type_label)), c("Variance", "Correlation"))
  expect_equal(unique(as.character(p_custom_labels$data$component_label)), c(
    "Gold-Oil covariance",
    "Gold-Oil correlation"
  ))

  expect_error(
    plot(out, type = "CIRF", components = 3),
    "Requested component is not available for the selected IRF type(s).",
    fixed = TRUE
  )
})
