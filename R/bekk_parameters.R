bekk_model_type <- function(bekk_model) {
  type <- NULL

  if (!is.null(bekk_model$spec$model$type)) {
    type <- bekk_model$spec$model$type
  }

  if (is.null(type)) {
    classes <- class(bekk_model)
    if (any(classes %in% c("sbekk", "sbekka"))) {
      type <- "sbekk"
    } else if (any(classes %in% c("dbekk", "dbekka"))) {
      type <- "dbekk"
    } else if (any(classes %in% c("bekk", "bekka"))) {
      type <- "bekk"
    }
  }

  if (is.null(type) || length(type) == 0L || is.na(type[1L])) {
    return(NA_character_)
  }

  as.character(type[1L])
}

bekk_parameter_matrix <- function(bekk_model,
                                  matrix_name,
                                  scalar_name,
                                  K,
                                  required = TRUE) {
  if (!is.null(bekk_model[[matrix_name]])) {
    mat <- as.matrix(bekk_model[[matrix_name]])
  } else if (!is.null(bekk_model[[scalar_name]])) {
    scalar <- as.numeric(bekk_model[[scalar_name]])
    if (length(scalar) != 1L || is.na(scalar)) {
      stop("`bekk_model$", scalar_name, "` must be a single numeric scalar.")
    }
    mat <- diag(scalar, K)
  } else if (required) {
    stop(
      "`bekk_model` must contain either `",
      matrix_name,
      "` or scalar `",
      scalar_name,
      "`."
    )
  } else {
    return(NULL)
  }

  if (nrow(mat) != K || ncol(mat) != K) {
    stop("`", matrix_name, "` must be a square matrix compatible with `data`.")
  }

  mat
}

bekk_extract_parameters <- function(bekk_model, K) {
  C0 <- bekk_parameter_matrix(
    bekk_model = bekk_model,
    matrix_name = "C0",
    scalar_name = "c0",
    K = K
  )
  A <- bekk_parameter_matrix(
    bekk_model = bekk_model,
    matrix_name = "A",
    scalar_name = "a",
    K = K
  )
  G <- bekk_parameter_matrix(
    bekk_model = bekk_model,
    matrix_name = "G",
    scalar_name = "g",
    K = K
  )

  asym <- isTRUE(bekk_model$asymmetric) ||
    !is.null(bekk_model$B) ||
    !is.null(bekk_model$b)

  if (asym) {
    B <- bekk_parameter_matrix(
      bekk_model = bekk_model,
      matrix_name = "B",
      scalar_name = "b",
      K = K
    )
  } else {
    B <- matrix(0, K, K)
  }

  signs <- bekk_model$signs
  if (is.null(signs)) {
    signs <- rep(-1, K)
  }
  signs <- as.numeric(signs)
  if (length(signs) != K || !all(signs %in% c(-1, 1))) {
    stop("`signs` must contain one -1 or 1 value per series.")
  }

  list(
    C0 = C0,
    A = A,
    G = G,
    B = B,
    asymmetric = asym,
    signs = signs,
    model_type = bekk_model_type(bekk_model)
  )
}
