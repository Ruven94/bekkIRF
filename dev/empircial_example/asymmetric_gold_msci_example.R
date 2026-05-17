# Empirical asymmetric BEKK IRF example.
#
# This is an internal development script, not a package test. It estimates an
# asymmetric BEKK model on the package data, computes all available IRFs, runs a
# small parameter bootstrap, and checks bootstrap-based IRF confidence bands.

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_file, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  file.path(normalizePath("."), "dev", "empircial_example")
}
root <- normalizePath(file.path(script_dir, "..", ".."))
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the package root or via Rscript dev/empircial_example/asymmetric_gold_msci_example.R.")
}

if (!requireNamespace("BEKKs", quietly = TRUE)) {
  stop("Package `BEKKs` is required.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package `devtools` is required for loading the package source.")
}

devtools::load_all(root, quiet = TRUE)
load(file.path(root, "data", "gold_msci_returns.rda"))

example_dir <- file.path(root, "dev", "empircial_example")
cache_dir <- file.path(example_dir, "cache")
figures_dir <- file.path(example_dir, "figures")
results_dir <- file.path(example_dir, "results")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Edit these for the full empirical run.
root_type <- "spectral"
shock_type <- "empirical"
time <- 444L
n.ahead <- 100L
simsamp <- 100000L
bootsamp <- 100L
seed <- 20260513L
ci_level <- 0.95
bekk_max_iter <- 50L
bootstrap_max_iter <- 5L

cores <- parallel::detectCores()

cat("Empirical asymmetric BEKK IRF example\n")
cat("Data rows:", nrow(gold_msci_returns), "Series:", ncol(gold_msci_returns), "\n")
cat("simsamp:", simsamp, "n.ahead:", n.ahead, "bootsamp:", bootsamp, "cores:", cores, "\n")

spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = TRUE))

fit_file <- file.path(cache_dir, "bekk_asymmetric_gold_msci.rds")
irf_file <- file.path(cache_dir, "irf_asymmetric_gold_msci_no_bootstrap.rds")
boot_file <- file.path(cache_dir, "bekk_bootstrap_asymmetric_gold_msci.rds")
boot_irf_file <- file.path(cache_dir, "irf_asymmetric_gold_msci_with_bootstrap.rds")
plot_data_file <- file.path(results_dir, "irf_asymmetric_gold_msci_plot_data.csv")
summary_file <- file.path(results_dir, "irf_asymmetric_gold_msci_summary.csv")
plot_file <- file.path(figures_dir, "irf_asymmetric_gold_msci_all_irfs.jpeg")
ci_plot_file <- file.path(figures_dir, "irf_asymmetric_gold_msci_all_irfs_with_ci.jpeg")

if (file.exists(fit_file)) {
  cat("Loading cached asymmetric BEKK fit:", fit_file, "\n")
  bekk_model <- readRDS(fit_file)
} else {
  cat("Estimating asymmetric BEKK model...\n")
  bekk_model <- BEKKs::bekk_fit(
    spec,
    gold_msci_returns,
    max_iter = bekk_max_iter
  )
  saveRDS(bekk_model, fit_file)
}



cat("Computing main IRFs without bootstrap...\n")
t <- Sys.time()
irf_main <- compute_irf(
  bekk_model,
  root_type = root_type,
  shock_type = shock_type,
  time = time,
  simsamp = simsamp,
  n.ahead = n.ahead,
  seed = seed,
  calc_virf = TRUE,
  calc_cirf = TRUE,
  calc_kirf = TRUE,
  calc_sirf = TRUE,
  calc_wirf = TRUE
)
print(Sys.time() - t)

saveRDS(irf_main, irf_file)

irf_matrix_to_long <- function(x, irf_type, value_name = "value") {
  if (is.null(x)) {
    return(NULL)
  }
  component_names <- colnames(x)
  if (is.null(component_names)) {
    component_names <- paste0(irf_type, "_", seq_len(ncol(x)))
  }
  out <- data.frame(
    horizon = rep(seq_len(nrow(x)), times = ncol(x)),
    irf_type = irf_type,
    component = rep(component_names, each = nrow(x)),
    value = as.vector(x),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "value"] <- value_name
  out
}

irf_long <- do.call(
  rbind,
  Filter(
    Negate(is.null),
    list(
      irf_matrix_to_long(irf_main$VIRF_mean, "VIRF"),
      irf_matrix_to_long(irf_main$CIRF_mean, "CIRF"),
      irf_matrix_to_long(irf_main$KIRF_mean, "KIRF"),
      irf_matrix_to_long(irf_main$SIRF_mean, "SIRF"),
      irf_matrix_to_long(irf_main$WIRF_mean, "WIRF")
    )
  )
)

p <- ggplot2::ggplot(irf_long, ggplot2::aes(x = horizon, y = value)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  ggplot2::geom_line(linewidth = 0.7, color = "#1f6f8b") +
  ggplot2::facet_wrap(irf_type ~ component, scales = "free_y") +
  ggplot2::labs(
    title = "Asymmetric BEKK impulse response functions",
    subtitle = paste0(
      "Gold/MSCI returns, ", shock_type, " shock at t = ", time,
      ", simsamp = ", format(simsamp, big.mark = ",")
    ),
    x = "Horizon",
    y = "IRF"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(plot_file, p, width = 13, height = 9, dpi = 160)
utils::write.csv(irf_long, plot_data_file, row.names = FALSE)

cat("Running asymmetric BEKK parameter bootstrap...\n")
bekk_boot <- bekk_bootstrap(
  bekk_model,
  bekk_spec_model = spec,
  data = gold_msci_returns,
  bootsamp = bootsamp,
  cores = cores,
  root_type = root_type,
  seed = seed + 1000L,
  max_iter = bootstrap_max_iter,
  center = TRUE,
  xi_outlier_threshold = 5,
  progress = TRUE,
  chunk_size = cores
)
class(bekk_boot)
saveRDS(bekk_boot, boot_file)

cat("Computing bootstrap-based IRF confidence intervals...\n")
irf_boot <- compute_irf(
  bekk_model,
  root_type = root_type,
  shock_type = shock_type,
  time = time,
  simsamp = simsamp,
  n.ahead = n.ahead,
  seed = seed,
  calc_virf = TRUE,
  calc_cirf = TRUE,
  calc_kirf = TRUE,
  calc_sirf = TRUE,
  calc_wirf = TRUE,
  bekk_bootstrap = bekk_boot,
  ci_level = ci_level,
  bootstrap_cores = cores,
  bootstrap_progress = TRUE,
  bootstrap_chunk_size = cores
)
class(irf_boot)
names(irf_boot)
saveRDS(irf_boot, boot_irf_file)

ci_to_long <- function(irf_obj, irf_type) {
  mean_mat <- irf_obj[[paste0(irf_type, "_mean")]]
  ci_obj <- irf_obj$ci[[irf_type]]
  if (is.null(mean_mat) || is.null(ci_obj)) {
    return(NULL)
  }
  center <- irf_matrix_to_long(mean_mat, irf_type, "mean")
  lower <- irf_matrix_to_long(ci_obj$lower, irf_type, "lower")
  upper <- irf_matrix_to_long(ci_obj$upper, irf_type, "upper")
  Reduce(
    function(x, y) merge(x, y, by = c("horizon", "irf_type", "component"), all = TRUE),
    list(center, lower, upper)
  )
}

ci_long <- do.call(
  rbind,
  Filter(
    Negate(is.null),
    lapply(c("VIRF", "CIRF", "KIRF", "SIRF", "WIRF"), function(type) ci_to_long(irf_boot, type))
  )
)

p_ci <- ggplot2::ggplot(ci_long, ggplot2::aes(x = horizon, y = mean)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower, ymax = upper),
    fill = "#8ecae6",
    alpha = 0.35
  ) +
  ggplot2::geom_line(linewidth = 0.7, color = "#023047") +
  ggplot2::facet_wrap(irf_type ~ component, scales = "free_y") +
  ggplot2::labs(
    title = "Asymmetric BEKK impulse response functions with bootstrap intervals",
    subtitle = paste0(
      "Gold/MSCI returns, ", bootsamp, " bootstrap draws, ",
      round(100 * ci_level), "% pointwise CI"
    ),
    x = "Horizon",
    y = "IRF"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(ci_plot_file, p_ci, width = 13, height = 9, dpi = 160)

summary <- data.frame(
  root_type = root_type,
  shock_type = shock_type,
  time = time,
  n.ahead = n.ahead,
  simsamp = simsamp,
  bootsamp = bootsamp,
  bootstrap_converged = sum(bekk_boot$converged),
  bootstrap_failed = sum(!bekk_boot$converged),
  irf_bootstrap_used = irf_boot$bootstrap_info$used_replications,
  irf_bootstrap_candidates = irf_boot$bootstrap_info$candidate_replications,
  irf_bootstrap_cores = irf_boot$bootstrap_info$cores,
  ci_lower = irf_boot$bootstrap_info$ci_probs[1L],
  ci_upper = irf_boot$bootstrap_info$ci_probs[2L]
)
utils::write.csv(summary, summary_file, row.names = FALSE)
print(summary)

cat("Cached asymmetric BEKK fit:", fit_file, "\n")
cat("Cached main IRF object:", irf_file, "\n")
cat("Cached bootstrap object:", boot_file, "\n")
cat("Cached bootstrap IRF object:", boot_irf_file, "\n")
cat("IRF plot:", plot_file, "\n")
cat("Bootstrap CI plot:", ci_plot_file, "\n")
cat("Summary:", summary_file, "\n")


# ---------------------------------------------------------------------------
# Load cached outputs
# ---------------------------------------------------------------------------
# Run this block on its own if you only want to inspect previously computed
# results without re-estimating the model or recomputing IRFs.


bekk_model_cached <- readRDS(fit_file)
irf_main_cached <- readRDS(irf_file)
bekk_boot_cached <- readRDS(boot_file)
irf_boot_cached <- readRDS(boot_irf_file)
summary_cached <- utils::read.csv(summary_file)
plot_data_cached <- utils::read.csv(plot_data_file)

irf_boot_cached$settings$series_names
irf_boot_cached$settings$series_names <- colnames(bekk_model_cached$data)
irf_main_cached$settings$series_names <- colnames(bekk_model_cached$data)

print(irf_boot_cached)
summary(irf_boot_cached)

# ---------------------------------------------------------------------------
# Plot function examples
# ---------------------------------------------------------------------------
# These examples cover all arguments of plot.bekkIRF(). Run them after loading
# the cached objects above.

plot(irf_boot_cached)

plot(
  irf_boot_cached,
  type = "virf",                 # case-insensitive: "virf", "Virf", "VIRF"
  ci = TRUE,                     # show bootstrap confidence interval ribbons
  components = c(1, 2),          # own component: 1; pair component: c(1, 2)
  title = "Gold/MSCI VIRF",
  subtitle = "Empirical shock at t = 444",
  xlab = "Horizon",
  ylab = "VIRF",
  xlim = c(1, 100),
  ylim = c(-1, 2),
  type_label = "Variance impulse response",
  component_labels = "MSCI-Gold covariance",
  line_color = "black",
  ci_fill = "blue",
  ci_alpha = 0.25,
  line_width = 0.9,
  zero_line = TRUE
)

plot(
  irf_boot_cached,
  type = "CIRF",
  ci = TRUE,
  title = "Gold/MSCI CIRF",
  ylab = "Correlation IRF",
  ylim = c(-0.5, 1),
  line_color = "#023047",
  ci_fill = "#8ecae6",
  ci_alpha = 0.35,
  zero_line = TRUE
)

plot(
  irf_boot_cached,
  type = "all",
  ci = FALSE,
  title = "All BEKK impulse response functions",
  type_label = c(
    VIRF = "Variance",
    CIRF = "Correlation",
    SIRF = "Skewness",
    KIRF = "Kurtosis",
    WIRF = "GMV weights"
  ),
  line_color = "black",
  line_width = 0.7,
  zero_line = TRUE
)

plot(
  irf_boot_cached,
  type = "VIRF",
  component_labels = c(
    "VIRF_{msci.Ret}" = "MSCI variance",
    "VIRF_{msci.Ret, gold.Ret}" = "MSCI-Gold covariance",
    "VIRF_{gold.Ret}" = "Gold variance"
  )
)

plot(
  irf_main_cached,
  type = "VIRF",
  component_labels = c(
    "VIRF_{msci.Ret}" = "MSCI variance",
    "VIRF_{msci.Ret, gold.Ret}" = "MSCI-Gold covariance",
    "VIRF_{gold.Ret}" = "Gold variance"
  )
)

