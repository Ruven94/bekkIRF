# Empirical Gold / MSCI World Developed Markets BEKK IRF example
#
# Internal development script, not a package test.
# The script estimates asymmetric and symmetric BEKK models on the package data,
# computes simulation-based IRFs, bootstrap confidence intervals, and a
# theoretical-vs-simulation VIRF comparison for the symmetric model.

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

required_packages <- c("BEKKs", "devtools", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), ".")
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

# ---------------------------------------------------------------------------
# User parameters
# ---------------------------------------------------------------------------
# Increase these for the paper run, e.g. simsamp = 100000 and bootsamp = 1000.
use_cache <- FALSE
root_type <- "spectral"
shock_type <- "empirical"
shock_time <- 444L
n.ahead <- 100
simsamp <- 100000
bootsamp <- 1000
ci_level <- 0.95
seed <- 20260518L
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
cores <- max(1L, detected_cores - 1L)
bootstrap_chunk_size <- cores

bekk_max_iter_asym <- 50L
bekk_max_iter_sym <- 50L
bootstrap_max_iter <- 10L
bootstrap_center <- TRUE
xi_outlier_threshold <- 5

# ---------------------------------------------------------------------------
# File locations
# ---------------------------------------------------------------------------
data_head_file <- file.path(results_dir, "gold_msci_returns_head.csv")
run_summary_file <- file.path(results_dir, "empirical_example_run_summary.csv")
asym_fit_file <- file.path(cache_dir, "bekk_asymmetric_gold_msci.rds")
asym_irf_file <- file.path(cache_dir, "irf_asymmetric_gold_msci_no_bootstrap.rds")
asym_boot_file <- file.path(cache_dir, "bekk_bootstrap_asymmetric_gold_msci.rds")
asym_boot_irf_file <- file.path(cache_dir, "irf_asymmetric_gold_msci_with_bootstrap.rds")
sym_fit_file <- file.path(cache_dir, "bekk_symmetric_gold_msci.rds")
sym_irf_file <- file.path(cache_dir, "irf_symmetric_gold_msci_simulated_virf.rds")
sym_theoretical_file <- file.path(cache_dir, "irf_symmetric_gold_msci_theoretical_virf.rds")
sym_compare_file <- file.path(results_dir, "symmetric_virf_theoretical_vs_simulated.csv")

plot_asym_all_file <- file.path(figures_dir, "asymmetric_irf_all_no_ci.jpeg")
plot_asym_all_ci_file <- file.path(figures_dir, "asymmetric_irf_all_with_ci.jpeg")
plot_asym_type_file <- function(type) {
  file.path(figures_dir, paste0("asymmetric_irf_", tolower(type), "_with_ci.jpeg"))
}
plot_sym_virf_compare_file <- file.path(figures_dir, "symmetric_virf_theoretical_vs_simulated.jpeg")

# ---------------------------------------------------------------------------
# Small timing and caching helpers
# ---------------------------------------------------------------------------
timings <- data.frame(step = character(), seconds = numeric(), stringsAsFactors = FALSE)

time_step <- function(label, expr) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("=", 72), "\n", sep = "")
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  timings <<- rbind(timings, data.frame(step = label, seconds = elapsed))
  cat(sprintf("Runtime: %.2f seconds\n", elapsed))
  value
}

load_or_compute <- function(file, label, expr) {
  if (isTRUE(use_cache) && file.exists(file)) {
    cat("Loading cached file:", file, "\n")
    return(readRDS(file))
  }

  value <- time_step(label, expr)
  saveRDS(value, file)
  cat("Saved:", file, "\n")
  value
}

save_plot <- function(plot, file, width = 13, height = 9, dpi = 170) {
  ggplot2::ggsave(file, plot, width = width, height = height, dpi = dpi)
  cat("Saved plot:", file, "\n")
  invisible(file)
}

cat("Empirical Gold / MSCI World Developed Markets BEKK IRF example\n")
cat("Rows:", nrow(gold_msci_returns), "Columns:", ncol(gold_msci_returns), "\n")
cat("Series:", paste(colnames(gold_msci_returns), collapse = ", "), "\n")
cat("Settings: simsamp =", simsamp, ", n.ahead =", n.ahead,
    ", bootsamp =", bootsamp, ", cores =", cores, "\n")

# ---------------------------------------------------------------------------
# Data head
# ---------------------------------------------------------------------------
data_head <- utils::head(gold_msci_returns, 10L)
print(data_head)
utils::write.csv(data_head, data_head_file)
cat("Saved data head:", data_head_file, "\n")

# ---------------------------------------------------------------------------
# Asymmetric BEKK model, simulation IRFs, bootstrap, bootstrap IRFs
# ---------------------------------------------------------------------------
asym_spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = TRUE))

asym_fit <- load_or_compute(
  asym_fit_file,
  "Estimate asymmetric BEKK model",
  BEKKs::bekk_fit(
    asym_spec,
    gold_msci_returns,
    max_iter = bekk_max_iter_asym
  )
)

asym_irf <- load_or_compute(
  asym_irf_file,
  "Compute asymmetric simulation IRFs without bootstrap",
  compute_irf(
    asym_fit,
    root_type = root_type,
    shock_type = shock_type,
    time = shock_time,
    simsamp = simsamp,
    n.ahead = n.ahead,
    seed = seed,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_kirf = TRUE,
    calc_sirf = TRUE,
    calc_wirf = TRUE
  )
)

asym_boot <- load_or_compute(
  asym_boot_file,
  "Run asymmetric BEKK parameter bootstrap",
  bekk_bootstrap(
    asym_fit,
    bekk_spec_model = asym_spec,
    data = gold_msci_returns,
    bootsamp = bootsamp,
    cores = cores,
    root_type = root_type,
    seed = seed + 1000L,
    max_iter = bootstrap_max_iter,
    center = bootstrap_center,
    xi_outlier_threshold = xi_outlier_threshold,
    progress = TRUE,
    chunk_size = bootstrap_chunk_size
  )
)

asym_boot_irf <- load_or_compute(
  asym_boot_irf_file,
  "Compute asymmetric bootstrap IRF confidence intervals",
  compute_irf(
    asym_fit,
    root_type = root_type,
    shock_type = shock_type,
    time = shock_time,
    simsamp = simsamp,
    n.ahead = n.ahead,
    seed = seed,
    calc_virf = TRUE,
    calc_cirf = TRUE,
    calc_kirf = TRUE,
    calc_sirf = TRUE,
    calc_wirf = TRUE,
    bekk_bootstrap = asym_boot,
    ci_level = ci_level,
    bootstrap_cores = cores,
    bootstrap_progress = TRUE,
    bootstrap_chunk_size = bootstrap_chunk_size
  )
)

# ---------------------------------------------------------------------------
# Symmetric BEKK model: theoretical VIRF vs simulation-based VIRF
# ---------------------------------------------------------------------------
sym_spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = FALSE))

sym_fit <- load_or_compute(
  sym_fit_file,
  "Estimate symmetric BEKK model",
  BEKKs::bekk_fit(
    sym_spec,
    gold_msci_returns,
    max_iter = bekk_max_iter_sym
  )
)

sym_irf <- load_or_compute(
  sym_irf_file,
  "Compute symmetric simulation-based VIRF",
  compute_irf(
    sym_fit,
    root_type = root_type,
    shock_type = shock_type,
    time = shock_time,
    simsamp = simsamp,
    n.ahead = n.ahead,
    seed = seed,
    calc_virf = TRUE,
    calc_cirf = FALSE,
    calc_kirf = FALSE,
    calc_sirf = FALSE,
    calc_wirf = FALSE
  )
)

sym_theoretical_virf <- load_or_compute(
  sym_theoretical_file,
  "Compute symmetric theoretical VIRF",
  compute_theoretical_virf(
    sym_fit,
    root_type = root_type,
    shock_type = shock_type,
    time = shock_time,
    n.ahead = n.ahead,
    format = "vech"
  )
)

virf_component_names <- colnames(sym_irf$VIRF_mean)
if (is.null(virf_component_names)) {
  virf_component_names <- paste0("VIRF_", seq_len(ncol(sym_irf$VIRF_mean)))
}
colnames(sym_theoretical_virf) <- virf_component_names

sym_compare <- rbind(
  data.frame(
    horizon = rep(seq_len(nrow(sym_irf$VIRF_mean)), times = ncol(sym_irf$VIRF_mean)),
    component = rep(virf_component_names, each = nrow(sym_irf$VIRF_mean)),
    method = "Simulation-based",
    value = as.vector(sym_irf$VIRF_mean),
    stringsAsFactors = FALSE
  ),
  data.frame(
    horizon = rep(seq_len(nrow(sym_theoretical_virf)), times = ncol(sym_theoretical_virf)),
    component = rep(virf_component_names, each = nrow(sym_theoretical_virf)),
    method = "Theoretical",
    value = as.vector(sym_theoretical_virf),
    stringsAsFactors = FALSE
  )
)
utils::write.csv(sym_compare, sym_compare_file, row.names = FALSE)
cat("Saved symmetric VIRF comparison data:", sym_compare_file, "\n")

# ---------------------------------------------------------------------------
# Reload cached outputs
# ---------------------------------------------------------------------------
# Run this block by itself when the files already exist and you only want to
# inspect or plot the stored results without re-estimating/recomputing.

cat("\n", strrep("-", 72), "\n", sep = "")
cat("Reload cached outputs\n")
cat(strrep("-", 72), "\n", sep = "")

gold_msci_head_cached <- utils::read.csv(data_head_file, row.names = 1L, check.names = FALSE)
asym_fit_cached <- readRDS(asym_fit_file)
asym_irf_cached <- readRDS(asym_irf_file)
asym_boot_cached <- readRDS(asym_boot_file)
asym_boot_irf_cached <- readRDS(asym_boot_irf_file)
sym_fit_cached <- readRDS(sym_fit_file)
sym_irf_cached <- readRDS(sym_irf_file)
sym_theoretical_virf_cached <- readRDS(sym_theoretical_file)
sym_compare_cached <- utils::read.csv(sym_compare_file)

print(asym_irf_cached)
summary(asym_boot_irf_cached)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
plot_asym_all <- plot(
  asym_irf_cached,
  type = "all",
  ci = FALSE,
  title = "Asymmetric BEKK impulse response functions",
  subtitle = paste0("Empirical shock at t = ", shock_time),
  line_color = "black",
  zero_line = TRUE
)
save_plot(plot_asym_all, plot_asym_all_file)

plot_asym_all_ci <- plot(
  asym_boot_irf_cached,
  type = "all",
  ci = TRUE,
  title = "Asymmetric BEKK impulse response functions",
  subtitle = paste0(round(100 * ci_level), "% bootstrap confidence intervals"),
  line_color = "black",
  ci_fill = "blue",
  ci_alpha = 0.22,
  zero_line = TRUE
)
save_plot(plot_asym_all_ci, plot_asym_all_ci_file)

for (irf_type in c("VIRF", "CIRF", "SIRF", "KIRF", "WIRF")) {
  p_type <- plot(
    asym_boot_irf_cached,
    type = irf_type,
    ci = TRUE,
    title = paste(irf_type, "impulse response function"),
    subtitle = paste0(round(100 * ci_level), "% bootstrap confidence intervals"),
    line_color = "black",
    ci_fill = "blue",
    ci_alpha = 0.22,
    zero_line = TRUE
  )
  save_plot(p_type, plot_asym_type_file(irf_type), width = 12, height = 7)
}

plot_sym_compare <- ggplot2::ggplot(
  sym_compare_cached,
  ggplot2::aes(x = horizon, y = value, color = method, linetype = method)
) +
  ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.45) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::facet_wrap(~ component, scales = "free_y") +
  ggplot2::scale_color_manual(values = c("Simulation-based" = "black", "Theoretical" = "#0072B2")) +
  ggplot2::labs(
    title = "Symmetric BEKK VIRF: theoretical vs simulation-based",
    subtitle = paste0("Empirical shock at t = ", shock_time, ", simsamp = ", simsamp),
    x = "Horizon",
    y = "VIRF",
    color = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )
save_plot(plot_sym_compare, plot_sym_virf_compare_file, width = 12, height = 7)

# ---------------------------------------------------------------------------
# Run summary
# ---------------------------------------------------------------------------
run_summary <- data.frame(
  root_type = root_type,
  shock_type = shock_type,
  shock_time = shock_time,
  n.ahead = n.ahead,
  simsamp = simsamp,
  bootsamp = bootsamp,
  ci_level = ci_level,
  cores = cores,
  asym_bootstrap_converged = sum(asym_boot_cached$converged),
  asym_bootstrap_failed = sum(!asym_boot_cached$converged),
  asym_irf_bootstrap_used = asym_boot_irf_cached$bootstrap_info$used_replications,
  asym_irf_bootstrap_requested = asym_boot_irf_cached$bootstrap_info$requested_replications,
  stringsAsFactors = FALSE
)
utils::write.csv(run_summary, run_summary_file, row.names = FALSE)
print(run_summary)

cat("\nCached files:\n")
cat("  data head:", data_head_file, "\n")
cat("  asymmetric BEKK fit:", asym_fit_file, "\n")
cat("  asymmetric IRF:", asym_irf_file, "\n")
cat("  asymmetric bootstrap:", asym_boot_file, "\n")
cat("  asymmetric bootstrap IRF:", asym_boot_irf_file, "\n")
cat("  symmetric BEKK fit:", sym_fit_file, "\n")
cat("  symmetric simulated VIRF:", sym_irf_file, "\n")
cat("  symmetric theoretical VIRF:", sym_theoretical_file, "\n")
cat("  symmetric VIRF comparison:", sym_compare_file, "\n")

cat("\nFigures:\n")
cat("  all asymmetric IRFs:", plot_asym_all_file, "\n")
cat("  all asymmetric IRFs with CI:", plot_asym_all_ci_file, "\n")
cat("  symmetric VIRF comparison:", plot_sym_virf_compare_file, "\n")

cat("\nTimings:\n")
print(timings, row.names = FALSE)


