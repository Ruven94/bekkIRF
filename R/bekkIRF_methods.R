bekk_irf_type_order <- c("VIRF", "CIRF", "SIRF", "KIRF", "WIRF")

bekk_irf_available_types <- function(x) {
  types <- bekk_irf_type_order[
    vapply(bekk_irf_type_order, function(type) {
      !is.null(x[[paste0(type, "_mean")]])
    }, logical(1L))
  ]
  types
}

bekk_irf_infer_k <- function(n_components, type) {
  if (type %in% c("VIRF", "KIRF")) {
    k <- (sqrt(8 * n_components + 1) - 1) / 2
  } else if (type == "CIRF") {
    k <- (1 + sqrt(1 + 8 * n_components)) / 2
  } else {
    k <- n_components
  }

  if (!is.finite(k) || abs(k - round(k)) > sqrt(.Machine$double.eps)) {
    return(n_components)
  }
  as.integer(round(k))
}

bekk_irf_series_names <- function(x, mat = NULL, type = NULL) {
  series_names <- x$settings$series_names
  if (!is.null(series_names)) {
    return(as.character(series_names))
  }

  if (!is.null(mat) && !is.null(type)) {
    return(paste0("Series ", seq_len(bekk_irf_infer_k(ncol(mat), type))))
  }

  "Series 1"
}

bekk_irf_component_label <- function(type, series_names, i, j = i) {
  if (i == j) {
    paste0(type, "_{", series_names[i], "}")
  } else {
    paste0(type, "_{", series_names[i], ", ", series_names[j], "}")
  }
}

bekk_irf_component_meta <- function(mat, type, series_names) {
  n_components <- ncol(mat)
  K <- length(series_names)

  if (type %in% c("VIRF", "KIRF")) {
    meta <- data.frame(
      component_id = seq_len(n_components),
      i = integer(n_components),
      j = integer(n_components),
      stringsAsFactors = FALSE
    )

    idx <- 1L
    for (col in seq_len(K)) {
      for (row in col:K) {
        if (idx <= n_components) {
          meta$i[idx] <- col
          meta$j[idx] <- row
          idx <- idx + 1L
        }
      }
    }
  } else if (type == "CIRF") {
    meta <- data.frame(
      component_id = seq_len(n_components),
      i = integer(n_components),
      j = integer(n_components),
      stringsAsFactors = FALSE
    )

    idx <- 1L
    for (i in seq_len(max(K - 1L, 0L))) {
      for (j in (i + 1L):K) {
        if (idx <= n_components) {
          meta$i[idx] <- i
          meta$j[idx] <- j
          idx <- idx + 1L
        }
      }
    }
  } else {
    meta <- data.frame(
      component_id = seq_len(n_components),
      i = seq_len(n_components),
      j = seq_len(n_components),
      stringsAsFactors = FALSE
    )
  }

  component_names <- colnames(mat)
  if (is.null(component_names)) {
    component_names <- vapply(
      seq_len(nrow(meta)),
      function(k) bekk_irf_component_label(type, series_names, meta$i[k], meta$j[k]),
      character(1L)
    )
  }

  meta$component <- component_names
  meta
}

bekk_irf_matrix_to_long <- function(mat, type, series_names, value_name = "value") {
  if (is.null(mat)) {
    return(NULL)
  }
  meta <- bekk_irf_component_meta(mat, type, series_names)

  out <- data.frame(
    horizon = rep(seq_len(nrow(mat)), times = ncol(mat)),
    type = type,
    component_id = rep(meta$component_id, each = nrow(mat)),
    i = rep(meta$i, each = nrow(mat)),
    j = rep(meta$j, each = nrow(mat)),
    component = rep(meta$component, each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "value"] <- value_name
  out
}

bekk_irf_normalize_type <- function(type) {
  if (length(type) != 1L || is.na(type)) {
    stop("`type` must be a single IRF type.")
  }
  type_upper <- toupper(type)
  if (type_upper == "ALL") {
    return("all")
  }
  valid <- bekk_irf_type_order
  match_idx <- match(type_upper, valid)
  if (is.na(match_idx)) {
    stop("`type` must be one of: all, ", paste(valid, collapse = ", "), ".")
  }
  valid[match_idx]
}

bekk_irf_normalize_components <- function(components) {
  if (is.null(components)) {
    return(NULL)
  }
  if (!is.numeric(components) || anyNA(components) || length(components) > 2L || length(components) < 1L) {
    stop("`components` must be NULL, a single component index, or a length-two numeric pair.")
  }
  if (any(components < 1) || any(components != as.integer(components))) {
    stop("`components` must contain positive integer indices.")
  }
  as.integer(components)
}

bekk_irf_filter_components <- function(plot_data, components) {
  if (is.null(components)) {
    return(plot_data)
  }

  if (length(components) == 1L) {
    i <- components[1L]
    keep <- plot_data$i == i & plot_data$j == i
  } else {
    i <- min(components)
    j <- max(components)
    keep <- plot_data$i == i & plot_data$j == j
  }

  plot_data[keep, , drop = FALSE]
}

bekk_irf_validate_limits <- function(limits, name) {
  if (is.null(limits)) {
    return(NULL)
  }
  if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) || limits[1L] >= limits[2L]) {
    stop("`", name, "` must be NULL or a numeric vector with lower < upper.")
  }
  limits
}

bekk_irf_apply_labels <- function(values, labels, name) {
  values <- as.character(values)
  unique_values <- unique(values)

  if (is.null(labels)) {
    return(factor(values, levels = unique_values))
  }
  if (!is.character(labels) || anyNA(labels)) {
    stop("`", name, "` must be NULL or a character vector.")
  }

  if (!is.null(names(labels)) && any(nzchar(names(labels)))) {
    mapped <- labels[values]
    missing <- unique(values[is.na(mapped)])
    if (length(missing) > 0L) {
      stop(
        "`", name, "` is missing labels for: ",
        paste(missing, collapse = ", "),
        "."
      )
    }
    return(factor(unname(mapped), levels = unique(unname(labels[unique_values]))))
  }

  if (length(labels) == 1L) {
    return(factor(rep(labels, length(values)), levels = labels))
  }
  if (length(labels) != length(unique_values)) {
    stop(
      "`", name, "` must have length 1, length equal to the number of plotted panels, ",
      "or be a named character vector."
    )
  }

  mapped <- stats::setNames(labels, unique_values)[values]
  factor(unname(mapped), levels = labels)
}

bekk_irf_plot_data <- function(x, types, ci = TRUE, components = NULL) {
  plot_data <- do.call(
    rbind,
    lapply(types, function(type) {
      mat <- x[[paste0(type, "_mean")]]
      center <- bekk_irf_matrix_to_long(
        mat,
        type,
        bekk_irf_series_names(x, mat, type),
        "mean"
      )
      if (is.null(center)) {
        return(NULL)
      }

      if (isTRUE(ci) && !is.null(x$ci[[type]])) {
        lower <- bekk_irf_matrix_to_long(x$ci[[type]]$lower, type, bekk_irf_series_names(x, mat, type), "lower")
        upper <- bekk_irf_matrix_to_long(x$ci[[type]]$upper, type, bekk_irf_series_names(x, mat, type), "upper")
        by_cols <- c("horizon", "type", "component_id", "i", "j", "component")
        center <- merge(center, lower, by = by_cols, all.x = TRUE)
        center <- merge(center, upper, by = by_cols, all.x = TRUE)
      } else {
        center$lower <- NA_real_
        center$upper <- NA_real_
      }

      center
    })
  )

  plot_data <- bekk_irf_filter_components(plot_data, components)
  if (nrow(plot_data) == 0L) {
    stop("Requested component is not available for the selected IRF type(s).")
  }

  plot_data$type <- factor(plot_data$type, levels = bekk_irf_type_order)
  plot_data$component <- factor(plot_data$component, levels = unique(plot_data$component))
  plot_data
}

bekk_irf_setting_value <- function(settings, names, default = "not stored") {
  for (name in names) {
    value <- settings[[name]]
    if (!is.null(value) && length(value) > 0L) {
      return(value)
    }
  }
  default
}

bekk_irf_format_scalar <- function(value, big_mark = FALSE) {
  if (is.null(value) || length(value) == 0L) {
    return("not stored")
  }
  if (is.numeric(value) || is.integer(value)) {
    value <- value[1L]
    if (is.na(value)) {
      return("not stored")
    }
    if (big_mark) {
      return(format(value, big.mark = ",", scientific = FALSE, trim = TRUE))
    }
    return(as.character(value))
  }
  value <- as.character(value[1L])
  if (is.na(value) || !nzchar(value)) {
    return("not stored")
  }
  value
}

bekk_irf_format_vector <- function(value, digits = 6L) {
  if (is.null(value) || length(value) == 0L) {
    return(NULL)
  }
  value <- as.numeric(value)
  if (anyNA(value)) {
    return(NULL)
  }
  paste0("c(", paste(signif(value, digits), collapse = ", "), ")")
}

bekk_irf_settings_line <- function(x) {
  settings <- x$settings
  if (is.null(settings) || length(settings) == 0L) {
    return(NULL)
  }

  asym_text <- NULL
  if (!is.null(settings$asymmetric)) {
    asym_text <- if (isTRUE(settings$asymmetric)) "yes" else "no"
  }

  line <- paste0(
    "root_type = ", bekk_irf_format_scalar(bekk_irf_setting_value(settings, c("root_type", "root"))),
    ", shock_type = ", bekk_irf_format_scalar(bekk_irf_setting_value(settings, c("shock_type"))),
    ", time = ", bekk_irf_format_scalar(bekk_irf_setting_value(settings, c("time"))),
    ", simsamp = ", bekk_irf_format_scalar(bekk_irf_setting_value(settings, c("simsamp")), big_mark = TRUE),
    ", n.ahead = ", bekk_irf_format_scalar(bekk_irf_setting_value(settings, c("n.ahead", "n_ahead")))
  )
  if (!is.null(asym_text)) {
    line <- paste0(line, ", asymmetric = ", asym_text)
  }
  line
}

bekk_irf_shock_line <- function(x) {
  settings <- x$settings
  if (is.null(settings) || length(settings) == 0L) {
    return(NULL)
  }

  shock <- bekk_irf_format_vector(bekk_irf_setting_value(settings, c("shock", "shock_vector"), default = NULL))
  if (is.null(shock)) {
    return(NULL)
  }

  paste0("shock = ", shock)
}

#' Print a BEKK IRF object
#'
#' @param x Object of class `"bekkIRF"`.
#' @param ... Additional arguments ignored.
#' @export
print.bekkIRF <- function(x, ...) {
  types <- bekk_irf_available_types(x)

  cat("bekkIRF object\n")
  settings_line <- bekk_irf_settings_line(x)
  if (!is.null(settings_line)) {
    cat(settings_line, "\n", sep = "")
  }
  shock_line <- bekk_irf_shock_line(x)
  if (!is.null(shock_line)) {
    cat(shock_line, "\n", sep = "")
  }

  if (length(types) == 0L) {
    cat("IRFs: none\n")
  } else {
    cat("IRFs: ", paste(types, collapse = ", "), "\n", sep = "")
  }

  if (!is.null(x$bootstrap_info)) {
    cat(
      "Bootstrap: ",
      x$bootstrap_info$used_replications,
      " used of ",
      x$bootstrap_info$requested_replications,
      " draws",
      sep = ""
    )
    if (!is.null(x$bootstrap_info$ci_probs)) {
      cat(
        ", CI = [",
        x$bootstrap_info$ci_probs[1L],
        ", ",
        x$bootstrap_info$ci_probs[2L],
        "]",
        sep = ""
      )
    }
    cat("\n")
  } else {
    cat("Bootstrap: none\n")
  }

  cat("Use summary() for details and plot() for IRF plots.\n")

  invisible(x)
}

#' Summarise a BEKK IRF object
#'
#' @param object Object of class `"bekkIRF"`.
#' @param ... Additional arguments ignored.
#'
#' @details
#' The summary table reports simple diagnostics for each selected IRF type:
#' `max_abs` is the largest absolute response across all horizons and
#' components, while `final_mean_abs` is the mean absolute response at the last
#' simulated horizon across components.
#' @export
summary.bekkIRF <- function(object, ...) {
  types <- bekk_irf_available_types(object)

  irf_summary <- do.call(
    rbind,
    lapply(types, function(type) {
      mat <- object[[paste0(type, "_mean")]]
      data.frame(
        type = type,
        horizons = nrow(mat),
        components = ncol(mat),
        min = min(mat, na.rm = TRUE),
        max = max(mat, na.rm = TRUE),
        max_abs = max(abs(mat), na.rm = TRUE),
        final_mean_abs = mean(abs(mat[nrow(mat), ]), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  out <- list(
    settings = object$settings,
    irf = irf_summary,
    bootstrap = object$bootstrap_info
  )
  class(out) <- c("summary_bekkIRF", "list")
  out
}

#' Print a BEKK IRF summary
#'
#' @param x Object returned by [summary.bekkIRF()].
#' @param ... Additional arguments ignored.
#' @export
print.summary_bekkIRF <- function(x, ...) {
  cat("Summary of bekkIRF object\n")

  if (!is.null(x$settings)) {
    cat(bekk_irf_settings_line(x), "\n", sep = "")
    shock_line <- bekk_irf_shock_line(x)
    if (!is.null(shock_line)) {
      cat(shock_line, "\n", sep = "")
    }
  }

  if (!is.null(x$irf) && nrow(x$irf) > 0L) {
    print(x$irf, row.names = FALSE)
  } else {
    cat("No IRF matrices available.\n")
  }

  if (!is.null(x$bootstrap)) {
    cat(
      "Bootstrap: ",
      x$bootstrap$used_replications,
      " used of ",
      x$bootstrap$requested_replications,
      " draws",
      sep = ""
    )
    if (!is.null(x$bootstrap$cores)) {
      cat(", cores = ", x$bootstrap$cores, sep = "")
    }
    cat("\n")
  } else {
    cat("Bootstrap: none\n")
  }

  invisible(x)
}

#' Plot BEKK impulse response functions
#'
#' Plots one or more impulse response functions stored in a `"bekkIRF"` object.
#' Bootstrap confidence intervals are shown automatically when available and
#' `ci = TRUE`.
#'
#' @param x Object of class `"bekkIRF"`.
#' @param y Ignored.
#' @param type IRF type to plot. Use `"all"` for all available IRFs, or one of
#'   `"VIRF"` (variance impulse response function), `"CIRF"` (correlation
#'   impulse response function), `"SIRF"` (skewness impulse response function),
#'   `"KIRF"` (kurtosis impulse response function), or `"WIRF"` (weights
#'   impulse response function for optimal portfolio weights). Matching is
#'   case-insensitive.
#' @param ci Logical. If `TRUE`, plot bootstrap confidence intervals when
#'   available.
#' @param components Optional component selector. Use a single integer for an
#'   own-series component, e.g. `1`, or a length-two pair such as `c(1, 2)`.
#' @param title,subtitle Plot title and subtitle. If `title = NULL`, an
#'   informative default title is used. If `subtitle = NULL`, no subtitle is
#'   shown.
#' @param xlab,ylab Axis labels.
#' @param xlim,ylim Optional axis limits passed to `ggplot2::coord_cartesian()`.
#' @param type_label Optional custom label(s) for the IRF type strip. Use a
#'   single character value, a vector in plotted type order, or a named vector.
#' @param component_labels Optional custom label(s) for component strips. Use a
#'   single character value, a vector in plotted component order, or a named
#'   vector mapping default component labels to custom labels.
#' @param line_color Line color.
#' @param ci_fill Fill color for bootstrap confidence ribbons.
#' @param ci_alpha Alpha transparency for bootstrap confidence ribbons.
#' @param line_width Line width.
#' @param zero_line Logical. If `TRUE`, draw a horizontal zero line.
#' @param ... Additional arguments ignored.
#'
#' @returns A `ggplot` object.
#'
#' @examples
#' K <- 2
#' N <- 30
#' x <- matrix(seq_len(N * K) / 100, nrow = N, ncol = K)
#' colnames(x) <- c("series1", "series2")
#' H_t <- matrix(rep(as.vector(diag(K)), times = N), nrow = N, byrow = TRUE)
#'
#' fit <- list(
#'   H_t = H_t,
#'   data = x,
#'   C0 = diag(0.1, K),
#'   A = diag(0.1, K),
#'   G = diag(0.8, K),
#'   asymmetric = FALSE
#' )
#'
#' irf <- compute_irf(fit, shock = c(1, 0), time = 10, simsamp = 10, n.ahead = 5)
#' plot(irf, type = "VIRF")
#' @export
plot.bekkIRF <- function(x,
                         y = NULL,
                         type = "all",
                         ci = TRUE,
                         components = NULL,
                         title = NULL,
                         subtitle = NULL,
                         xlab = "Horizon",
                         ylab = "IRF",
                         xlim = NULL,
                         ylim = NULL,
                         type_label = NULL,
                         component_labels = NULL,
                         line_color = "#023047",
                         ci_fill = "#8ecae6",
                         ci_alpha = 0.35,
                         line_width = 0.75,
                         zero_line = TRUE,
                         ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.")
  }

  type <- bekk_irf_normalize_type(type)
  components <- bekk_irf_normalize_components(components)
  available_types <- bekk_irf_available_types(x)
  if (length(available_types) == 0L) {
    stop("No IRF matrices are available to plot.")
  }

  types <- if (type == "all") available_types else type
  missing_types <- setdiff(types, available_types)
  if (length(missing_types) > 0L) {
    stop("Requested IRF type is not available: ", paste(missing_types, collapse = ", "), ".")
  }

  if (length(ci) != 1L || is.na(ci)) {
    stop("`ci` must be TRUE or FALSE.")
  }
  if (length(zero_line) != 1L || is.na(zero_line)) {
    stop("`zero_line` must be TRUE or FALSE.")
  }
  xlim <- bekk_irf_validate_limits(xlim, "xlim")
  ylim <- bekk_irf_validate_limits(ylim, "ylim")

  plot_data <- bekk_irf_plot_data(x, types, ci = isTRUE(ci), components = components)
  plot_data$type_label <- bekk_irf_apply_labels(plot_data$type, type_label, "type_label")
  plot_data$component_label <- bekk_irf_apply_labels(plot_data$component, component_labels, "component_labels")
  has_ci <- isTRUE(ci) && any(is.finite(plot_data$lower)) && any(is.finite(plot_data$upper))

  if (is.null(title)) {
    title <- if (type == "all") {
      "BEKK impulse response functions"
    } else {
      paste(type, "impulse response function")
    }
  }
  if (is.null(subtitle)) {
    subtitle <- NULL
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = horizon, y = mean))

  if (isTRUE(zero_line)) {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.45)
  }

  if (has_ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = ci_fill,
        alpha = ci_alpha,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_line(color = line_color, linewidth = line_width, na.rm = TRUE) +
    ggplot2::facet_wrap(type_label ~ component_label, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }

  p
}
