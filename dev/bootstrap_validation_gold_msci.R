# Internal validation script:
# Run a small BEKK parameter bootstrap on the package data and check the
# function behavior without adding a slow stochastic test to testthat.

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_file, value = TRUE)
root <- if (length(file_arg) == 1L) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."))
} else {
  normalizePath(".")
}
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the package root or via Rscript dev/bootstrap_validation_gold_msci.R.")
}

source(file.path(root, "R", "matrix_utils.R"))
source(file.path(root, "R", "compute_xi.R"))
source(file.path(root, "R", "compute_eta.R"))
source(file.path(root, "R", "bekk_bootstrap.R"))

load(file.path(root, "data", "gold_msci_returns.rda"))

dir.create(file.path(root, "dev", "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "dev", "results"), recursive = TRUE, showWarnings = FALSE)

spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = FALSE))
fit_file <- file.path(root, "dev", "cache", "bekk_symmetric_gold_msci.rds")

if (file.exists(fit_file)) {
  bekk_model <- readRDS(fit_file)
} else {
  bekk_model <- BEKKs::bekk_fit(spec, gold_msci_returns, max_iter = 50)
  saveRDS(bekk_model, fit_file)
}

expect_true <- function(x, label) {
  if (!isTRUE(x)) {
    stop("Validation failed: ", label, call. = FALSE)
  }
}

summarise_bootstrap <- function(x, label) {
  finite_params <- c(as.vector(x$C0), as.vector(x$A), as.vector(x$G))
  if (!is.null(x$B)) {
    finite_params <- c(finite_params, as.vector(x$B))
  }

  data.frame(
    label = label,
    bootsamp = x$settings$bootsamp,
    cores = x$settings$cores,
    root_type = x$settings$root_type,
    max_iter = x$settings$max_iter,
    converged = sum(x$converged),
    failed = sum(!x$converged),
    error_count = sum(!is.na(x$error)),
    any_nonfinite_parameter = any(!is.finite(finite_params), na.rm = TRUE),
    c0_dim = paste(dim(x$C0), collapse = "x"),
    a_dim = paste(dim(x$A), collapse = "x"),
    g_dim = paste(dim(x$G), collapse = "x"),
    b_is_null = is.null(x$B),
    stringsAsFactors = FALSE
  )
}

K <- ncol(gold_msci_returns)
main_bootsamp <- 50L
detected_cores <- parallel::detectCores()
if (is.na(detected_cores) || detected_cores < 2L) {
  main_cores <- 1L
} else {
  main_cores <- min(2L, detected_cores - 1L)
}

cat("Running main bootstrap validation with", main_bootsamp, "replications...\n")
main_boot <- bekk_bootstrap(
  bekk_model,
  bekk_spec_model = spec,
  data = gold_msci_returns,
  bootsamp = main_bootsamp,
  cores = main_cores,
  root_type = "spectral",
  seed = 20260505L,
  max_iter = 5,
  progress = TRUE,
  chunk_size = main_cores
)

expect_true(inherits(main_boot, "bekkBootstrap"), "main output class")
expect_true(identical(dim(main_boot$C0), c(K, K, main_bootsamp)), "main C0 dimensions")
expect_true(identical(dim(main_boot$A), c(K, K, main_bootsamp)), "main A dimensions")
expect_true(identical(dim(main_boot$G), c(K, K, main_bootsamp)), "main G dimensions")
expect_true(is.null(main_boot$B), "main symmetric B is NULL")
expect_true(length(main_boot$converged) == main_bootsamp, "main converged length")
expect_true(length(main_boot$error) == main_bootsamp, "main error length")
expect_true(any(main_boot$converged), "at least one main bootstrap fit succeeded")

cat("\nRunning serial argument check...\n")
serial_boot <- bekk_bootstrap(
  bekk_model,
  bootsamp = 3L,
  cores = 1L,
  root_type = "spec",
  seed = 101L,
  max_iter = 1,
  progress = FALSE
)
expect_true(identical(dim(serial_boot$C0), c(K, K, 3L)), "serial C0 dimensions")
expect_true(serial_boot$settings$cores == 1L, "serial cores setting")

cat("Running parallel argument check...\n")
parallel_boot <- bekk_bootstrap(
  bekk_model,
  bekk_spec_model = spec,
  data = gold_msci_returns,
  bootsamp = 2L,
  cores = 2L,
  root_type = "spectral",
  seed = 202L,
  max_iter = 1,
  progress = FALSE,
  chunk_size = 1L
)
expect_true(identical(dim(parallel_boot$C0), c(K, K, 2L)), "parallel C0 dimensions")
expect_true(parallel_boot$settings$cores == 2L, "parallel cores setting")

cat("Running reproducibility check...\n")
repro_1 <- bekk_bootstrap(
  bekk_model,
  bootsamp = 2L,
  cores = 1L,
  root_type = "spectral",
  seed = 909L,
  max_iter = 1,
  progress = FALSE
)
repro_2 <- bekk_bootstrap(
  bekk_model,
  bootsamp = 2L,
  cores = 1L,
  root_type = "spectral",
  seed = 909L,
  max_iter = 1,
  progress = FALSE
)
expect_true(isTRUE(all.equal(repro_1$C0, repro_2$C0, tolerance = 0)), "C0 reproducibility")
expect_true(isTRUE(all.equal(repro_1$A, repro_2$A, tolerance = 0)), "A reproducibility")
expect_true(isTRUE(all.equal(repro_1$G, repro_2$G, tolerance = 0)), "G reproducibility")

cat("Running cholesky root alias check...\n")
chol_boot <- bekk_bootstrap(
  bekk_model,
  bekk_spec_model = spec,
  data = gold_msci_returns,
  bootsamp = 2L,
  cores = 1L,
  root_type = "cholesky",
  seed = 303L,
  max_iter = 1,
  progress = FALSE
)
expect_true(chol_boot$settings$root_type == "chol", "cholesky alias maps to chol")
expect_true(identical(dim(chol_boot$G), c(K, K, 2L)), "cholesky G dimensions")

cat("Running input validation checks...\n")
err_bootsamp <- tryCatch(
  bekk_bootstrap(bekk_model, bootsamp = 0L, progress = FALSE),
  error = function(e) e
)
err_cores <- tryCatch(
  bekk_bootstrap(bekk_model, bootsamp = 1L, cores = 0L, progress = FALSE),
  error = function(e) e
)
expect_true(inherits(err_bootsamp, "error"), "invalid bootsamp errors")
expect_true(inherits(err_cores, "error"), "invalid cores errors")

summary <- rbind(
  summarise_bootstrap(main_boot, "main_50_spectral_parallel"),
  summarise_bootstrap(serial_boot, "serial_3_spec"),
  summarise_bootstrap(parallel_boot, "parallel_2_spectral"),
  summarise_bootstrap(repro_1, "repro_1"),
  summarise_bootstrap(repro_2, "repro_2"),
  summarise_bootstrap(chol_boot, "chol_2")
)

parameter_summary <- data.frame(
  parameter = c(
    paste0("C0[", row(main_boot$C0[, , 1L]), ",", col(main_boot$C0[, , 1L]), "]"),
    paste0("A[", row(main_boot$A[, , 1L]), ",", col(main_boot$A[, , 1L]), "]"),
    paste0("G[", row(main_boot$G[, , 1L]), ",", col(main_boot$G[, , 1L]), "]")
  ),
  mean = c(
    apply(main_boot$C0, c(1, 2), mean, na.rm = TRUE),
    apply(main_boot$A, c(1, 2), mean, na.rm = TRUE),
    apply(main_boot$G, c(1, 2), mean, na.rm = TRUE)
  ),
  sd = c(
    apply(main_boot$C0, c(1, 2), stats::sd, na.rm = TRUE),
    apply(main_boot$A, c(1, 2), stats::sd, na.rm = TRUE),
    apply(main_boot$G, c(1, 2), stats::sd, na.rm = TRUE)
  ),
  q025 = c(
    apply(main_boot$C0, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE),
    apply(main_boot$A, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE),
    apply(main_boot$G, c(1, 2), stats::quantile, probs = 0.025, na.rm = TRUE)
  ),
  q975 = c(
    apply(main_boot$C0, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE),
    apply(main_boot$A, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE),
    apply(main_boot$G, c(1, 2), stats::quantile, probs = 0.975, na.rm = TRUE)
  )
)

summary_file <- file.path(root, "dev", "results", "bootstrap_validation_summary.csv")
parameter_file <- file.path(root, "dev", "results", "bootstrap_validation_parameter_summary.csv")
object_file <- file.path(root, "dev", "results", "bootstrap_validation_main_bootstrap.rds")

utils::write.csv(summary, summary_file, row.names = FALSE)
utils::write.csv(parameter_summary, parameter_file, row.names = FALSE)
saveRDS(main_boot, object_file)

print(summary)
cat("Parameter summary written to:", parameter_file, "\n")
cat("Validation summary written to:", summary_file, "\n")
cat("Main bootstrap object written to:", object_file, "\n")
