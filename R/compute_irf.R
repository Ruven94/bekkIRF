

compute_irf <- function(bekk_model,
                        # bootstrap_bekk_model,
                        root_type = c("spec", "chol"),
                        shock_type = c("structural", "empirical"),
                        shock = NULL,
                        time = NULL,
                        simsamp = 10000,
                        n.ahead = 100,
                        seed = 123){

  xi <- compute_xi(bekk_model$H_t, bekk_model$data, root_type = root_type)
  psi_skew <- compute_psi_skewness(xi)
  psi_kurt <- compute_psi_kurtosis(xi)

  # Parameter
  A <- bekk_model$A
  if (bekk_model$asymmetric){
    B <- bekk_model$B
  }
  C <- bekk_model$C0
  G <- bekk_model$G
  K <- ncol(bekk_model$data)
  N <- nrow(bekk_model$data)
  signs <- bekk_model$signs

  H_0 <- matrix(bekk_model$H_t[time,], K, K)

  if(shock_type == "structural"){
    shock <- shock
  }else if(shock_type == "empirical"){
    shock <- xi[time,]
  }

  e <- matroot(H_0, type = root_type) %*% shock

  H_1 <- C %*% t(C) +
    t(A) %*% e %*% t(e) %*% A +
    t(G) %*% H_0 %*% G +
    if (bekk_model$asymmetric){
    t(B) %*% calculate_eta(e, signs = signs) %*% t(calculate_eta(e, signs = signs)) %*% B
    }

  # H_1 muss in die C++ function, da es hier auch bootstrapped parameter nehmen soll, nur shock bleibt gleich

  # if (bootstrap_bekk_model == NULL){C++}
  # if (bootstrap_bekk_model != NULL){parallelsierung einbauen}
  AVIRF_core_cpp(H_t  = H_t,
                 xi   = xi,
                 C    = C,
                 A    = A,
                 G    = G,
                 B    = B_mat,
                 Psi  = Psi,
                 psi  = psi,
                 timeforVIRF = timeforVIRF,
                 n_ahead     = n.ahead,
                 asym        = asym,
                 simsamp     = simsamp,
                 seed        = seed)
}


## --- Paramter ----
n.ahead <- 100
timeforVIRF <- 444

## ---- IRF parallel function ----
IRF_parallel <- function(data,                 # returns
                         C, A, G,              # BEKK-Parameter
                         B = NULL,             # only for asymmetric model
                         timeforVIRF,          # start index t
                         n.ahead,              # h-ahead forecast
                         asym = FALSE,         # Using sym / asym?
                         simsamp = NULL,       # simulation sample size
                         C_boot = NULL, A_boot = NULL,
                         G_boot = NULL, B_boot = NULL,
                         cores = parallel::detectCores(), # aktuell in C++ nicht mehr genutzt
                         seed = 123) {

  ## -- 1  Checks --------------------------------------------------------------
  if (is.null(B) && asym) stop("Asym-Modell benötigt B!")
  if (is.null(simsamp)) simsamp <- 2000 * n.ahead
  if (abs(mean(data[, 1])) > 1e-15 || abs(mean(data[, 2])) > 1e-15) {
    data <- scale(data, center = TRUE, scale = FALSE)
    message("Data was centered")
  }

  ## -- 2  BEKK-Reconstruction--------------------------------------------------
  N <- nrow(data)
  K <- ncol(data)

  H_t <- matrix(0, N, K^2)
  xi  <- matrix(0, N, K)

  H_t[1, ] <- t(data) %*% data / N
  for (i in 2:N) {
    eta <- calculate_eta(data[i - 1, ])
    H_t[i, ] <- as.vector(C %*% t(C) +
                            t(A) %*% data[i - 1, ] %*% t(data[i - 1, ]) %*% A +
                            t(G) %*% matrix(H_t[i - 1, ], K, K) %*% G)
    if (asym) {
      H_t[i, ] <- H_t[i, ] + as.vector(t(B) %*% eta %*% t(eta) %*% B)
    }
    xi[i, ] <- chol2inv(chol(matroot(matrix(H_t[i, ], K, K)))) %*% data[i, ]
  }

  ## -- 3  Optional: Bootstrap-Parameter --------------------------------------
  if (!is.null(C_boot) || !is.null(A_boot) || !is.null(G_boot) || !is.null(B_boot)) {
    C <- C_boot; A <- A_boot; G <- G_boot
    if (!is.null(B_boot)) B <- B_boot
  }

  if (is.null(B)) {
    B_mat <- matrix(0, K, K)
  } else {
    B_mat <- B
  }

  ## -- 4  Vorbereitung (Psi/psi in R, Rest in C++) ---------------------------
  Psi <- Psi_Kurtosis_function(xi)
  psi <- Psi_Skewness_function(xi)

  # C++-Kern aufrufen
  res <- AVIRF_core_cpp(H_t  = H_t,
                        xi   = xi,
                        C    = C,
                        A    = A,
                        G    = G,
                        B    = B_mat,
                        Psi  = Psi,
                        psi  = psi,
                        timeforVIRF = timeforVIRF,
                        n_ahead     = n.ahead,
                        asym        = asym,
                        simsamp     = simsamp,
                        seed        = seed)

  res
}
