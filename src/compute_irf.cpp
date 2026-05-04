#include <RcppArmadillo.h>
#include <random>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// -----------------------------------------------------------------------------
// Helper: matrix root
// root_type = 0 -> spectral
// root_type = 1 -> cholesky
// -----------------------------------------------------------------------------
inline arma::mat matroot_cpp(const arma::mat& H, const int root_type, const double tol = 1e-12) {
  if (root_type == 0) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, H);

    for (arma::uword i = 0; i < eigval.n_elem; ++i) {
      if (eigval(i) < 0.0 && eigval(i) > -tol) {
        eigval(i) = 0.0;
      }
    }

    if (arma::any(eigval < 0.0)) {
      Rcpp::stop("`H` is not positive semi-definite.");
    }

    arma::mat L = arma::diagmat(arma::sqrt(eigval));
    return eigvec * L * eigvec.t();
  }

  if (root_type == 1) {
    arma::mat out;
    bool ok = arma::chol(out, H, "lower");
    if (!ok) {
      Rcpp::stop("`H` must be positive definite for `root_type = \"chol\"`.");
    }
    return out;
  }

  Rcpp::stop("Invalid `root_type`.");
  return arma::mat();
}

// -----------------------------------------------------------------------------
// Helper: compute eta
// Must match the R function compute_eta() and BEKKs-consistent joint sign logic
// -----------------------------------------------------------------------------
inline arma::vec compute_eta_cpp(const arma::vec& e, const arma::vec& signs) {
  if (e.n_elem != signs.n_elem) {
    Rcpp::stop("`e` and `signs` must have the same length.");
  }

  for (arma::uword i = 0; i < signs.n_elem; ++i) {
    if (signs(i) != -1.0 && signs(i) != 1.0) {
      Rcpp::stop("`signs` must contain only -1 and 1.");
    }
  }

  arma::vec prod = signs % e;

  for (arma::uword i = 0; i < prod.n_elem; ++i) {
    if (prod(i) <= 0.0) {
      return arma::zeros<arma::vec>(e.n_elem);
    }
  }

  return e;
}

// -----------------------------------------------------------------------------
// Helper: vectorize lower triangle (vech)
// -----------------------------------------------------------------------------
inline arma::vec vech_cpp(const arma::mat& M) {
  int K = M.n_rows;
  int p = K * (K + 1) / 2;
  arma::vec out(p);
  int idx = 0;

  for (int col = 0; col < K; ++col) {
    for (int row = col; row < K; ++row) {
      out(idx++) = M(row, col);
    }
  }

  return out;
}

// -----------------------------------------------------------------------------
// Helper: update BEKK covariance
// -----------------------------------------------------------------------------
inline arma::mat update_bekk_cpp(const arma::mat& H_prev,
                                 const arma::vec& e,
                                 const arma::mat& C,
                                 const arma::mat& A,
                                 const arma::mat& G,
                                 const arma::mat& B,
                                 const arma::vec& signs,
                                 const bool asym) {
  arma::mat H_new = C * C.t() +
    A.t() * e * e.t() * A +
    G.t() * H_prev * G;

  if (asym) {
    arma::vec eta = compute_eta_cpp(e, signs);
    H_new += B.t() * eta * eta.t() * B;
  }

  return H_new;
}

// -----------------------------------------------------------------------------
// Helper: correlation IRF for all pairwise correlations
// -----------------------------------------------------------------------------
inline arma::rowvec cirf_row_cpp(const arma::mat& H1, const arma::mat& H0) {
  int K = H1.n_rows;
  int n_cirf = K * (K - 1) / 2;
  arma::rowvec out(n_cirf, fill::zeros);

  int idx = 0;
  for (int m = 0; m < K - 1; ++m) {
    for (int k = m + 1; k < K; ++k) {
      double rho1 = H1(k, m) / std::sqrt(H1(k, k) * H1(m, m));
      double rho0 = H0(k, m) / std::sqrt(H0(k, k) * H0(m, m));
      out(idx++) = rho1 - rho0;
    }
  }

  return out;
}

// -----------------------------------------------------------------------------
// Helper: GMV weights IRF
// -----------------------------------------------------------------------------
inline arma::rowvec wirf_row_cpp(const arma::mat& H1, const arma::mat& H0) {
  int K = H1.n_rows;
  arma::vec one = arma::ones<arma::vec>(K);

  arma::mat invH1 = arma::inv_sympd(H1);
  arma::vec w1 = invH1 * one;
  w1 /= arma::as_scalar(one.t() * w1);

  arma::mat invH0 = arma::inv_sympd(H0);
  arma::vec w0 = invH0 * one;
  w0 /= arma::as_scalar(one.t() * w0);

  return (w1 - w0).t();
}

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List compute_irf_core_cpp(const arma::mat& H_0,
                                const arma::vec& shock,
                                const arma::mat& xi,
                                const arma::mat& C,
                                const arma::mat& A,
                                const arma::mat& G,
                                const arma::mat& B,
                                const arma::vec& signs,
                                const arma::mat& psi_kurt,
                                const arma::mat& psi_skew,
                                const int root_type,
                                const bool asym,
                                const int simsamp,
                                const int n_ahead,
                                const int seed,
                                const bool calc_virf,
                                const bool calc_cirf,
                                const bool calc_kirf,
                                const bool calc_sirf,
                                const bool calc_wirf) {

  int N = xi.n_rows;
  int K = xi.n_cols;
  int K2 = K * K;
  int n_sel = K * (K + 1) / 2;
  int n_cirf = K * (K - 1) / 2;

  if (N <= 0) {
    Rcpp::stop("`xi` must have at least one row.");
  }

  if (K <= 0) {
    Rcpp::stop("`xi` must have at least one column.");
  }

  if (simsamp <= 0) {
    Rcpp::stop("`simsamp` must be positive.");
  }

  if (n_ahead <= 0) {
    Rcpp::stop("`n_ahead` must be positive.");
  }

  if (H_0.n_rows != (arma::uword)K || H_0.n_cols != (arma::uword)K) {
    Rcpp::stop("`H_0` must be a square matrix with dimensions compatible with `xi`.");
  }

  if (C.n_rows != (arma::uword)K || C.n_cols != (arma::uword)K ||
      A.n_rows != (arma::uword)K || A.n_cols != (arma::uword)K ||
      G.n_rows != (arma::uword)K || G.n_cols != (arma::uword)K ||
      B.n_rows != (arma::uword)K || B.n_cols != (arma::uword)K) {
    Rcpp::stop("`C`, `A`, `G`, and `B` must be square matrices compatible with `xi`.");
  }

  if ((int)shock.n_elem != K) {
    Rcpp::stop("`shock` must have length equal to `ncol(xi)`.");
  }

  if ((int)signs.n_elem != K) {
    Rcpp::stop("`signs` must have length equal to `ncol(xi)`.");
  }

  if (psi_kurt.n_rows != (arma::uword)K2 || psi_kurt.n_cols != (arma::uword)K2) {
    Rcpp::stop("`psi_kurt` must have dimensions `K^2 x K^2`.");
  }

  if (psi_skew.n_rows != (arma::uword)K2 || psi_skew.n_cols != (arma::uword)K) {
    Rcpp::stop("`psi_skew` must have dimensions `K^2 x K`.");
  }

  // Result containers
  arma::cube VIRF_all;
  arma::cube CIRF_all;
  arma::cube KIRF_all;
  arma::cube SIRF_all;
  arma::cube WIRF_all;

  if (calc_virf) VIRF_all = arma::cube(n_ahead, n_sel, simsamp, fill::zeros);
  if (calc_cirf) CIRF_all = arma::cube(n_ahead, n_cirf, simsamp, fill::zeros);
  if (calc_kirf) KIRF_all = arma::cube(n_ahead, n_sel, simsamp, fill::zeros);
  if (calc_sirf) SIRF_all = arma::cube(n_ahead, K, simsamp, fill::zeros);
  if (calc_wirf) WIRF_all = arma::cube(n_ahead, K, simsamp, fill::zeros);

  std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));
  std::uniform_int_distribution<int> draw_idx(0, N - 1);

  for (int s = 0; s < simsamp; ++s) {
    arma::ivec idx(n_ahead);
    for (int i = 0; i < n_ahead; ++i) {
      idx(i) = draw_idx(rng);
    }

    arma::mat H_base = H_0;
    arma::mat H_shock = H_0;

    arma::mat VIRF_j;
    arma::mat CIRF_j;
    arma::mat KIRF_j;
    arma::mat SIRF_j;
    arma::mat WIRF_j;

    if (calc_virf) VIRF_j = arma::mat(n_ahead, n_sel, fill::zeros);
    if (calc_cirf) CIRF_j = arma::mat(n_ahead, n_cirf, fill::zeros);
    if (calc_kirf) KIRF_j = arma::mat(n_ahead, n_sel, fill::zeros);
    if (calc_sirf) SIRF_j = arma::mat(n_ahead, K, fill::zeros);
    if (calc_wirf) WIRF_j = arma::mat(n_ahead, K, fill::zeros);

    // -----------------------------------------------------------------------
    // Step 1: baseline path
    // -----------------------------------------------------------------------
    arma::vec x0 = xi.row(idx(0)).t();
    arma::mat Q0 = matroot_cpp(H_base, root_type);
    arma::vec e_base = Q0 * x0;
    arma::mat H_base_next = update_bekk_cpp(H_base, e_base, C, A, G, B, signs, asym);

    // -----------------------------------------------------------------------
    // Step 1: shock path
    // -----------------------------------------------------------------------
    arma::vec e_shock = Q0 * shock;
    arma::mat H_shock_next = update_bekk_cpp(H_shock, e_shock, C, A, G, B, signs, asym);

    // -----------------------------------------------------------------------
    // IRFs at horizon 1
    // -----------------------------------------------------------------------
    arma::mat diffH = H_shock_next - H_base_next;

    if (calc_virf) {
      VIRF_j.row(0) = vech_cpp(diffH).t();
    }

    if (calc_cirf) {
      CIRF_j.row(0) = cirf_row_cpp(H_shock_next, H_base_next);
    }

    if (calc_wirf) {
      WIRF_j.row(0) = wirf_row_cpp(H_shock_next, H_base_next);
    }

    if (calc_kirf || calc_sirf) {
      arma::mat R1 = matroot_cpp(H_shock_next, root_type);
      arma::mat R0 = matroot_cpp(H_base_next, root_type);
      arma::mat Kp1 = arma::kron(R1, R1);
      arma::mat Kp0 = arma::kron(R0, R0);

      if (calc_kirf) {
        int js = 0;
        for (int col = 0; col < K; ++col) {
          for (int row = col; row < K; ++row) {
            int j0 = row + col * K;

            arma::rowvec row1 = Kp1.row(j0);
            arma::colvec col1 = Kp1.col(j0);
            arma::rowvec row0 = Kp0.row(j0);
            arma::colvec col0 = Kp0.col(j0);

            double num1 = arma::as_scalar(row1 * psi_kurt * col1);
            double num0 = arma::as_scalar(row0 * psi_kurt * col0);

            double den1 = H_shock_next(row, row) * H_shock_next(col, col);
            double den0 = H_base_next(row, row) * H_base_next(col, col);

            KIRF_j(0, js++) = num1 / den1 - num0 / den0;
          }
        }
      }

      if (calc_sirf) {
        for (int j = 0; j < K; ++j) {
          int col_idx = j + j * K;

          arma::rowvec r1 = R1.row(j);
          arma::rowvec r0 = R0.row(j);
          arma::colvec v1 = Kp1.col(col_idx);
          arma::colvec v0 = Kp0.col(col_idx);

          double s1 = arma::as_scalar(r1 * psi_skew.t() * v1);
          double s0 = arma::as_scalar(r0 * psi_skew.t() * v0);

          double h1jj = H_shock_next(j, j);
          double h0jj = H_base_next(j, j);

          SIRF_j(0, j) = s1 / std::pow(h1jj, 1.5) - s0 / std::pow(h0jj, 1.5);
        }
      }
    }

    H_base = H_base_next;
    H_shock = H_shock_next;

    // -----------------------------------------------------------------------
    // Horizons 2 ... n_ahead
    // -----------------------------------------------------------------------
    for (int i = 1; i < n_ahead; ++i) {
      arma::vec x_future = xi.row(idx(i)).t();

      arma::mat Q_base = matroot_cpp(H_base, root_type);
      arma::vec e_base_i = Q_base * x_future;
      H_base = update_bekk_cpp(H_base, e_base_i, C, A, G, B, signs, asym);

      arma::mat Q_shock = matroot_cpp(H_shock, root_type);
      arma::vec e_shock_i = Q_shock * x_future;
      H_shock = update_bekk_cpp(H_shock, e_shock_i, C, A, G, B, signs, asym);

      diffH = H_shock - H_base;

      if (calc_virf) {
        VIRF_j.row(i) = vech_cpp(diffH).t();
      }

      if (calc_cirf) {
        CIRF_j.row(i) = cirf_row_cpp(H_shock, H_base);
      }

      if (calc_wirf) {
        WIRF_j.row(i) = wirf_row_cpp(H_shock, H_base);
      }

      if (calc_kirf || calc_sirf) {
        arma::mat R1 = matroot_cpp(H_shock, root_type);
        arma::mat R0 = matroot_cpp(H_base, root_type);
        arma::mat Kp1 = arma::kron(R1, R1);
        arma::mat Kp0 = arma::kron(R0, R0);

        if (calc_kirf) {
          int js = 0;
          for (int col = 0; col < K; ++col) {
            for (int row = col; row < K; ++row) {
              int j0 = row + col * K;

              arma::rowvec row1 = Kp1.row(j0);
              arma::colvec col1 = Kp1.col(j0);
              arma::rowvec row0 = Kp0.row(j0);
              arma::colvec col0 = Kp0.col(j0);

              double num1 = arma::as_scalar(row1 * psi_kurt * col1);
              double num0 = arma::as_scalar(row0 * psi_kurt * col0);

              double den1 = H_shock(row, row) * H_shock(col, col);
              double den0 = H_base(row, row) * H_base(col, col);

              KIRF_j(i, js++) = num1 / den1 - num0 / den0;
            }
          }
        }

        if (calc_sirf) {
          for (int j = 0; j < K; ++j) {
            int col_idx = j + j * K;

            arma::rowvec r1 = R1.row(j);
            arma::rowvec r0 = R0.row(j);
            arma::colvec v1 = Kp1.col(col_idx);
            arma::colvec v0 = Kp0.col(col_idx);

            double s1 = arma::as_scalar(r1 * psi_skew.t() * v1);
            double s0 = arma::as_scalar(r0 * psi_skew.t() * v0);

            double h1jj = H_shock(j, j);
            double h0jj = H_base(j, j);

            SIRF_j(i, j) = s1 / std::pow(h1jj, 1.5) - s0 / std::pow(h0jj, 1.5);
          }
        }
      }
    }

    if (calc_virf) VIRF_all.slice(s) = VIRF_j;
    if (calc_cirf) CIRF_all.slice(s) = CIRF_j;
    if (calc_kirf) KIRF_all.slice(s) = KIRF_j;
    if (calc_sirf) SIRF_all.slice(s) = SIRF_j;
    if (calc_wirf) WIRF_all.slice(s) = WIRF_j;
  }

  // Means over simulation runs
  Rcpp::List out;

  if (calc_virf) {
    arma::mat VIRF_mean = arma::mean(VIRF_all, 2);
    out["VIRF_mean"] = VIRF_mean;
  } else {
    out["VIRF_mean"] = R_NilValue;
  }

  if (calc_cirf) {
    arma::mat CIRF_mean = arma::mean(CIRF_all, 2);
    out["CIRF_mean"] = CIRF_mean;
  } else {
    out["CIRF_mean"] = R_NilValue;
  }

  if (calc_kirf) {
    arma::mat KIRF_mean = arma::mean(KIRF_all, 2);
    out["KIRF_mean"] = KIRF_mean;
  } else {
    out["KIRF_mean"] = R_NilValue;
  }

  if (calc_sirf) {
    arma::mat SIRF_mean = arma::mean(SIRF_all, 2);
    out["SIRF_mean"] = SIRF_mean;
  } else {
    out["SIRF_mean"] = R_NilValue;
  }

  if (calc_wirf) {
    arma::mat WIRF_mean = arma::mean(WIRF_all, 2);
    out["WIRF_mean"] = WIRF_mean;
  } else {
    out["WIRF_mean"] = R_NilValue;
  }

  return out;
}
