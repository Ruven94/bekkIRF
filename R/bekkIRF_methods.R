bekk_irf_type_order <- c("VIRF", "CIRF", "SIRF", "KIRF", "WIRF")

bekk_irf_available_types <- function(x) {
  types <- bekk_irf_type_order[
    vapply(bekk_irf_type_order, function(type) {
      !is.null(x[[paste0(type, "_mean")]])
    }, logical(1L))
  ]
  types
}

bekk_irf_component_names <- function(mat, type) {
  component_names <- colnames(mat)
  if (is.null(component_names)) {
    component_names <- paste0(type, "_", seq_len(ncol(mat)))
  }
  component_names
}

bekk_irf_matrix_to_long <- function(mat, type, value_name = "value") {
  if (is.null(mat)) {
    return(NULL)
  }

  out <- data.frame(
    horizon = rep(seq_len(nrow(mat)), times = ncol(mat)),
    type = type,
    component = rep(bekk_irf_component_names(mat, type), each = nrow(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "value"] <- value_name
  out
}

bekk_irf_plot_data <- function(x, types, ci = TRUE) {
  plot_data <- do.call(
    rbind,
    lapply(types, function(type) {
      center <- bekk_irf_matrix_to_long(x[[paste0(type, "_mean")]], type, "mean")
      if (is.null(center)) {
        return(NULL)
      }

      if (isTRUE(ci) && !is.null(x$ci[[type]])) {
        lower <- bekk_irf_matrix_to_long(x$ci[[type]]$lower, type, "lower")
        upper <- bekk_irf_matrix_to_long(x$ci[[type]]$upper, type, "upper")
        center <- merge(center, lower, by = c("horizon", "type", "component"), all.x = TRUE)
        center <- merge(center, upper, by = c("horizon", "type", "component"), all.x = TRUE)
      } else {
        center$lower <- NA_real_
        center$upper <- NA_real_
      }

      center
    })
  )

  plot_data$type <- factor(plot_data$type, levels = bekk_irf_type_order)
  plot_data
}

bekk_irf_settings_line <- function(x) {
  settings <- x$settings
  if (is.null(settings)) {
    return(NULL)
  }

  paste0(
    "root = ", settings$root_type,
    ", shock = ", settings$shock_type,
    ", time = ", settings$time,
    ", simsamp = ", format(settings$simsamp, big.mark = ","),
    ", n.ahead = ", settings$n.ahead
  )
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

  if (length(types) == 0L) {
    cat("IRFs: none\n")
  } else {
    dims <- vapply(types, function(type) {
      paste(dim(x[[paste0(type, "_mean")]]), collapse = " x ")
    }, character(1L))
    cat("IRFs: ", paste(paste0(types, " (", dims, ")"), collapse = ", "), "\n", sep = "")
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

  invisible(x)
}

#' Summarise a BEKK IRF object
#'
#' @param object Object of class `"bekkIRF"`.
#' @param ... Additional arguments ignored.
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
        has_ci = !is.null(object$ci[[type]]),
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
    cat(
      "root = ", x$settings$root_type,
      ", shock = ", x$settings$shock_type,
      ", time = ", x$settings$time,
      ", simsamp = ", format(x$settings$simsamp, big.mark = ","),
      ", n.ahead = ", x$settings$n.ahead,
      "\n",
      sep = ""
    )
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
  }

  invisible(x)
}

#' Plot BEKK impulse response functions
#'
#' @param x Object of class `"bekkIRF"`.
#' @param y Ignored.
#' @param type IRF type to plot. Use `"all"` for all available IRFs, or one of
#'   `"VIRF"`, `"CIRF"`, `"SIRF"`, `"KIRF"`, `"WIRF"`.
#' @param ci Logical. If `TRUE`, plot bootstrap confidence intervals when
#'   available.
#' @param ... Additional arguments ignored.
#'
#' @returns A `ggplot` object.
#' @export
plot.bekkIRF <- function(x,
                         y = NULL,
                         type = c("all", "VIRF", "CIRF", "SIRF", "KIRF", "WIRF"),
                         ci = TRUE,
                         ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.")
  }

  type <- match.arg(type)
  available_types <- bekk_irf_available_types(x)
  if (length(available_types) == 0L) {
    stop("No IRF matrices are available to plot.")
  }

  types <- if (type == "all") available_types else type
  missing_types <- setdiff(types, available_types)
  if (length(missing_types) > 0L) {
    stop("Requested IRF type is not available: ", paste(missing_types, collapse = ", "), ".")
  }

  plot_data <- bekk_irf_plot_data(x, types, ci = ci)
  has_ci <- isTRUE(ci) && any(is.finite(plot_data$lower)) && any(is.finite(plot_data$upper))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = horizon, y = mean)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey78", linewidth = 0.3)

  if (has_ci) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = "#8ecae6",
        alpha = 0.35,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_line(color = "#023047", linewidth = 0.75, na.rm = TRUE) +
    ggplot2::facet_wrap(type ~ component, scales = "free_y") +
    ggplot2::labs(
      title = if (type == "all") "BEKK impulse response functions" else paste(type, "impulse response functions"),
      subtitle = bekk_irf_settings_line(x),
      x = "Horizon",
      y = "IRF"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}
