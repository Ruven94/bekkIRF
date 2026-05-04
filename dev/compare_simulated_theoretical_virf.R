# Internal validation script:
# Compare simulation-based VIRFs from compute_irf() with closed-form symmetric
# VIRFs from compute_theoretical_virf() for the package data.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

script_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
root <- if (!is.na(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."))
} else {
  normalizePath(".")
}
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  root <- normalizePath(".")
}

source(file.path(root, "R", "matrix_utils.R"))
source(file.path(root, "R", "compute_xi.R"))
source(file.path(root, "R", "moment_utils.R"))
source(file.path(root, "R", "compute_irf.R"))
source(file.path(root, "R", "compute_theoretical_virf.R"))
Rcpp::sourceCpp(file.path(root, "src", "compute_irf.cpp"))

load(file.path(root, "data", "gold_msci_returns.rda"))

dir.create(file.path(root, "dev", "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "dev", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "dev", "results"), recursive = TRUE, showWarnings = FALSE)

root_type <- "spectral"
shock_type <- "empirical"
time <- 444L
n.ahead <- 25L
simsamp <- 100000L
seed <- 123L

fit_file <- file.path(root, "dev", "cache", "bekk_symmetric_gold_msci.rds")

if (file.exists(fit_file)) {
  bekk_model <- readRDS(fit_file)
} else {
  spec <- BEKKs::bekk_spec(model = list(type = "bekk", asymmetric = FALSE))
  bekk_model <- BEKKs::bekk_fit(spec, gold_msci_returns, max_iter = 50)
  saveRDS(bekk_model, fit_file)
}

simulated <- compute_irf(
  bekk_model,
  root_type = root_type,
  shock_type = shock_type,
  time = time,
  simsamp = simsamp,
  n.ahead = n.ahead,
  seed = seed,
  calc_virf = TRUE,
  calc_cirf = FALSE,
  calc_kirf = FALSE,
  calc_sirf = FALSE,
  calc_wirf = FALSE
)$VIRF_mean

theoretical <- compute_theoretical_virf(
  bekk_model,
  root_type = root_type,
  shock_type = shock_type,
  time = time,
  n.ahead = n.ahead,
  format = "vech"
)

diff <- simulated - theoretical
comparison_summary <- data.frame(
  root_type = root_type,
  shock_type = shock_type,
  time = time,
  n.ahead = n.ahead,
  simsamp = simsamp,
  seed = seed,
  max_abs_diff = max(abs(diff)),
  mean_abs_diff = mean(abs(diff)),
  rmse = sqrt(mean(diff^2)),
  abs_equal_1e_3 = all(abs(diff) <= 1e-3),
  abs_equal_1e_2 = all(abs(diff) <= 1e-2),
  abs_equal_5e_2 = all(abs(diff) <= 5e-2)
)

print(comparison_summary)

virf_names <- colnames(theoretical)
if (is.null(virf_names)) {
  virf_names <- paste0("VIRF_", seq_len(ncol(theoretical)))
}

plot_data <- rbind(
  data.frame(
    horizon = rep(seq_len(n.ahead), times = ncol(simulated)),
    component = rep(virf_names, each = n.ahead),
    value = as.vector(simulated),
    type = "simulation"
  ),
  data.frame(
    horizon = rep(seq_len(n.ahead), times = ncol(theoretical)),
    component = rep(virf_names, each = n.ahead),
    value = as.vector(theoretical),
    type = "theoretical"
  )
)

summary_file <- file.path(root, "dev", "results", "simulation_vs_theoretical_virf_summary.csv")
values_file <- file.path(root, "dev", "results", "simulation_vs_theoretical_virf_values.csv")
utils::write.csv(comparison_summary, summary_file, row.names = FALSE)
utils::write.csv(plot_data, values_file, row.names = FALSE)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = horizon, y = value, color = type)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::facet_wrap(~ component, scales = "free_y") +
  ggplot2::labs(
    title = "Simulation vs theoretical symmetric VIRF",
    subtitle = paste0(
      "root = ", root_type,
      ", empirical shock at t = ", time,
      ", simsamp = ", format(simsamp, big.mark = ",")
    ),
    x = "Horizon",
    y = "VIRF",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")

plot_file <- file.path(root, "dev", "figures", "simulation_vs_theoretical_virf_spectral_t444.png")
ggplot2::ggsave(plot_file, p, width = 9, height = 5, dpi = 160)

if (interactive()) {
  print(p)
}

cat("Plot written to:", plot_file, "\n")
cat("Summary written to:", summary_file, "\n")
cat("Values written to:", values_file, "\n")
