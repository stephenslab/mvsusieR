// C++ implementations of hot loops for the per-condition (3d) missing
// data methods (approximate and exact).
//
// Each function has a corresponding R implementation in R/missing_y_utils.R
// and R/individual_data_methods.R.  The R dispatchers fall back to R when
// the C++ versions are unavailable (e.g., during devtools::load_all()
// without compilation).
//
// Functions:
//   compute_VinvR_3d_rcpp   - per-pattern precision application (N x R)
//   compute_XtR_3d_rcpp     - V^{-1}-weighted X'R cross-product  (J x R)
//   compute_Xb_3d_rcpp      - per-condition X*b with optional Xbar correction
//   compute_betahat_3d_rcpp - batch matrix-vector: betahat[j] = svs[j] * XtR[j]
//   compute_svs_inv_3d_rcpp - assemble J svs_inv matrices from precomputed sums
//   compute_vbxxb_rcpp      - sum_j alpha_j tr(svs_inv_j * mu2_j) scalar

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// 1. compute_VinvR_3d_rcpp
//
// Apply per-pattern V_i^{-1} to each row of an N x R matrix.
// pattern_assign is 1-indexed (R convention).
//
// mat:            N x R input matrix
// Vinv_list:      list of K  R x R  precision matrices
// pattern_assign: N-vector of pattern indices (1-based)
//
// Returns N x R matrix.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat compute_VinvR_3d_rcpp(const arma::mat& mat,
                                const Rcpp::List& Vinv_list,
                                const arma::ivec& pattern_assign) {
  unsigned int N = mat.n_rows;
  unsigned int R = mat.n_cols;
  int K = Vinv_list.size();

  arma::mat result(N, R, fill::zeros);

  // Cache Vinv matrices in a vector for fast lookup
  std::vector<arma::mat> Vinv(K);
  for (int k = 0; k < K; k++) {
    Vinv[k] = Rcpp::as<arma::mat>(Vinv_list[k]);
  }

  // Single pass over rows: look up pattern, apply Vinv
  for (unsigned int i = 0; i < N; i++) {
    int k = pattern_assign(i) - 1;  // 0-based
    // result[i,:] = mat[i,:] * Vinv[k]'  (row-vector times transpose)
    // Equivalent to Vinv[k] * mat[i,:]' but stored as row
    result.row(i) = mat.row(i) * Vinv[k].t();
  }

  return result;
}


// ---------------------------------------------------------------------------
// 2. compute_XtR_3d_rcpp
//
// Compute V^{-1}-weighted cross-product X'R for approximate or exact method.
//
// Steps:
//   1. VinvR[i,:] = Vinv_{pattern(i)} * R_mat[i,:]   (per-pattern precision)
//   2. XtR[j,r] = sum_i X_3d[i,j,r] * VinvR[i,r]    (per-condition cross-product)
//   3. (exact only) XtR[j,:] -= Xbar[j,,] * colSums(VinvR)
//
// X_3d:           N x J x R array (already centered/scaled, NAs -> 0)
// R_mat:          N x R residual matrix
// Vinv_list:      list of K  R x R  precision matrices
// pattern_assign: N-vector (1-based)
// method:         1 = approximate, 2 = exact
// Xbar:           J x R x R array (exact only; ignored for approximate)
//
// Returns J x R matrix.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat compute_XtR_3d_rcpp(const arma::cube& X_3d,
                              const arma::mat& R_mat,
                              const Rcpp::List& Vinv_list,
                              const arma::ivec& pattern_assign,
                              int method,
                              const arma::cube& Xbar) {
  unsigned int N = X_3d.n_rows;
  unsigned int J = X_3d.n_cols;
  unsigned int R = X_3d.n_slices;
  int K = Vinv_list.size();

  // Cache Vinv
  std::vector<arma::mat> Vinv(K);
  for (int k = 0; k < K; k++) {
    Vinv[k] = Rcpp::as<arma::mat>(Vinv_list[k]);
  }

  // Step 1: VinvR (N x R)
  arma::mat VinvR(N, R, fill::zeros);
  for (unsigned int i = 0; i < N; i++) {
    int k = pattern_assign(i) - 1;
    VinvR.row(i) = R_mat.row(i) * Vinv[k].t();
  }

  // Step 2: Per-condition cross-product
  // XtR[j,r] = sum_i X_3d[i,j,r] * VinvR[i,r]
  arma::mat XtR(J, R, fill::zeros);
  for (unsigned int r = 0; r < R; r++) {
    // X_3d.slice(r) is N x J; VinvR.col(r) is N x 1
    // XtR[:,r] = X_3d[:,:,r]' * VinvR[:,r]
    XtR.col(r) = X_3d.slice(r).t() * VinvR.col(r);
  }

  // Step 3 (exact only): Xbar correction
  if (method == 2) {
    arma::rowvec colsum_VinvR = sum(VinvR, 0);  // 1 x R
    for (unsigned int j = 0; j < J; j++) {
      // Xbar(j, :, :) is R x R -- need to extract it
      arma::mat Xbar_j(R, R);
      for (unsigned int r1 = 0; r1 < R; r1++) {
        for (unsigned int r2 = 0; r2 < R; r2++) {
          Xbar_j(r1, r2) = Xbar(j, r1, r2);
        }
      }
      XtR.row(j) -= colsum_VinvR * Xbar_j;
    }
  }

  return XtR;
}


// ---------------------------------------------------------------------------
// 3. compute_Xb_3d_rcpp
//
// Per-condition X*b with optional Xbar correction (exact method).
//
// X_3d:   N x J x R array
// b:      J x R coefficient matrix
// method: 1 = approximate, 2 = exact
// Xbar:   J x R x R array (exact only; ignored for approximate)
//
// Returns N x R matrix.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat compute_Xb_3d_rcpp(const arma::cube& X_3d,
                             const arma::mat& b,
                             int method,
                             const arma::cube& Xbar) {
  unsigned int N = X_3d.n_rows;
  unsigned int J = X_3d.n_cols;
  unsigned int R = X_3d.n_slices;

  arma::mat Xb(N, R, fill::zeros);

  // Per-condition: Xb[:,r] = X_3d[:,:,r] * b[:,r]
  for (unsigned int r = 0; r < R; r++) {
    Xb.col(r) = X_3d.slice(r) * b.col(r);
  }

  // Exact: subtract Xbar correction
  // Xbarb = sum_j Xbar[j,,] %*% b[j,:]  (R-vector, broadcast to all N rows)
  if (method == 2) {
    arma::vec Xbarb(R, fill::zeros);
    for (unsigned int j = 0; j < J; j++) {
      arma::mat Xbar_j(R, R);
      for (unsigned int r1 = 0; r1 < R; r1++) {
        for (unsigned int r2 = 0; r2 < R; r2++) {
          Xbar_j(r1, r2) = Xbar(j, r1, r2);
        }
      }
      Xbarb += Xbar_j * b.row(j).t();
    }
    // Subtract from each row
    Xb.each_row() -= Xbarb.t();
  }

  return Xb;
}


// ---------------------------------------------------------------------------
// 4. compute_betahat_3d_rcpp
//
// Batch GLS: betahat[j,:] = svs[j] * XtR[j,:] for j = 1..J.
//
// svs_3d: R x R x J array of per-variable covariance matrices
// XtR:    J x R matrix
//
// Returns J x R matrix (NaN replaced with 0).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat compute_betahat_3d_rcpp(const arma::cube& svs_3d,
                                  const arma::mat& XtR) {
  unsigned int J = XtR.n_rows;
  unsigned int R = XtR.n_cols;

  arma::mat bhat(J, R);
  for (unsigned int j = 0; j < J; j++) {
    bhat.row(j) = (svs_3d.slice(j) * XtR.row(j).t()).t();
  }

  // Replace NaN with 0
  bhat.replace(datum::nan, 0.0);

  return bhat;
}


// ---------------------------------------------------------------------------
// 5. compute_svs_inv_3d_rcpp
//
// Assemble all J svs_inv matrices from precomputed K x J summary matrices.
//
// For the approximate method:
//   svs_inv_j[r,s] = sum_k Vinv_k[r,s] *
//     (raw_sq_sum[k,j] - (cm[j,r]+cm[j,s])*raw_sum[k,j] + n_k*cm[j,r]*cm[j,s])
//     / (csd[j,r] * csd[j,s])
//   where sum is only over (r,s) both in obs_k.
//
// For the exact method: same A1 plus A2/Xbar correction terms.
//
// Vinv_list:   list of K  R x R  precision matrices
// pattern:     K x R logical (observed conditions per pattern)
// raw_sq_sum:  K x J matrix
// raw_sum:     K x J matrix
// n_k:         K-vector
// cm:          J x R per-condition column means
// csd:         J x R per-condition column SDs
// method:      1 = approximate, 2 = exact
// Xbar:        J x R x R array (exact only)
// Vinvsum:     R x R (exact only)
//
// Returns R x R x J array.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::cube compute_svs_inv_3d_rcpp(const Rcpp::List& Vinv_list,
                                   const arma::imat& pattern,
                                   const arma::mat& raw_sq_sum,
                                   const arma::mat& raw_sum,
                                   const arma::ivec& n_k,
                                   const arma::mat& cm,
                                   const arma::mat& csd,
                                   int method,
                                   const arma::cube& Xbar,
                                   const arma::mat& Vinvsum) {
  unsigned int K = pattern.n_rows;
  unsigned int R = pattern.n_cols;
  unsigned int J = raw_sq_sum.n_cols;

  // Cache Vinv and per-pattern observed indices
  std::vector<arma::mat> Vinv(K);
  std::vector<std::vector<unsigned int>> obs_r_list(K);
  for (unsigned int k = 0; k < K; k++) {
    Vinv[k] = Rcpp::as<arma::mat>(Vinv_list[k]);
    for (unsigned int r = 0; r < R; r++) {
      if (pattern(k, r)) obs_r_list[k].push_back(r);
    }
  }

  arma::cube result(R, R, J, fill::zeros);

  for (unsigned int j = 0; j < J; j++) {
    arma::mat A1(R, R, fill::zeros);

    if (method == 1) {
      // Approximate method
      for (unsigned int k = 0; k < K; k++) {
        const auto& obs_r = obs_r_list[k];
        unsigned int n_obs = obs_r.size();
        if (n_obs == 0) continue;

        for (unsigned int ri = 0; ri < n_obs; ri++) {
          unsigned int r = obs_r[ri];
          for (unsigned int si = ri; si < n_obs; si++) {
            unsigned int s = obs_r[si];
            double cross_rs = (raw_sq_sum(k, j) -
              (cm(j, r) + cm(j, s)) * raw_sum(k, j) +
              n_k(k) * cm(j, r) * cm(j, s)) /
              (csd(j, r) * csd(j, s));
            double val = Vinv[k](r, s) * cross_rs;
            A1(r, s) += val;
            if (r != s) A1(s, r) += val;
          }
        }
      }
      result.slice(j) = A1;
    } else {
      // Exact method
      arma::mat A2(R, R, fill::zeros);

      for (unsigned int k = 0; k < K; k++) {
        const auto& obs_r = obs_r_list[k];
        unsigned int n_obs = obs_r.size();
        if (n_obs == 0) continue;

        // A1 contribution
        for (unsigned int ri = 0; ri < n_obs; ri++) {
          unsigned int r = obs_r[ri];
          for (unsigned int si = ri; si < n_obs; si++) {
            unsigned int s = obs_r[si];
            double cross_rs = (raw_sq_sum(k, j) -
              (cm(j, r) + cm(j, s)) * raw_sum(k, j) +
              n_k(k) * cm(j, r) * cm(j, s)) /
              (csd(j, r) * csd(j, s));
            double val = Vinv[k](r, s) * cross_rs;
            A1(r, s) += val;
            if (r != s) A1(s, r) += val;
          }
        }

        // A2 contribution
        for (unsigned int ri = 0; ri < n_obs; ri++) {
          unsigned int r = obs_r[ri];
          double xsum_r = (raw_sum(k, j) - n_k(k) * cm(j, r)) / csd(j, r);
          A2.col(r) += Vinv[k].col(r) * xsum_r;
        }
      }

      // Extract Xbar_j (R x R)
      arma::mat Xbar_j(R, R);
      for (unsigned int r1 = 0; r1 < R; r1++) {
        for (unsigned int r2 = 0; r2 < R; r2++) {
          Xbar_j(r1, r2) = Xbar(j, r1, r2);
        }
      }

      result.slice(j) = A1 - Xbar_j.t() * A2 - A2.t() * Xbar_j +
        Xbar_j.t() * Vinvsum * Xbar_j;
    }
  }

  return result;
}


// ---------------------------------------------------------------------------
// 6. compute_vbxxb_rcpp
//
// Compute the vbxxb scalar from SER_posterior_e_loglik:
//   vbxxb = sum_j alpha[j] * tr(svs_inv[j] * mu2[j])
//
// Also computes bxxb (R x R):
//   bxxb = sum_j d[j] * alpha[j] * mu2[j]
//
// alpha:      J-vector
// mu2_3d:     J x R x R array of posterior second moments
// svs_inv_3d: R x R x J array
// d:          J-vector (X'X diagonal)
//
// Returns list(bxxb = R x R matrix, vbxxb = scalar)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List compute_vbxxb_rcpp(const arma::vec& alpha,
                              const arma::cube& mu2_3d,
                              const arma::cube& svs_inv_3d,
                              const arma::vec& d) {
  unsigned int J = alpha.n_elem;
  unsigned int R = mu2_3d.n_cols;

  arma::mat bxxb(R, R, fill::zeros);
  double vbxxb = 0.0;

  for (unsigned int j = 0; j < J; j++) {
    // mu2_j is stored as mu2_3d(j, :, :) -- extract R x R
    arma::mat mu2_j(R, R);
    for (unsigned int r1 = 0; r1 < R; r1++) {
      for (unsigned int r2 = 0; r2 < R; r2++) {
        mu2_j(r1, r2) = mu2_3d(j, r1, r2);
      }
    }

    double a_j = alpha(j);
    arma::mat pb2_j = a_j * mu2_j;

    bxxb += d(j) * pb2_j;

    // tr(svs_inv[j] * pb2_j) = accu(svs_inv[j] % pb2_j)
    vbxxb += accu(svs_inv_3d.slice(j) % pb2_j);
  }

  return Rcpp::List::create(
    Rcpp::Named("bxxb") = bxxb,
    Rcpp::Named("vbxxb") = vbxxb
  );
}


// ---------------------------------------------------------------------------
// 7. compute_Xbar_from_sums_rcpp
//
// Compute the Xbar correction array for the exact missing-data method
// using precomputed K x J summary matrices instead of an O(J*N*R^2) loop.
//
// For each variable j:
//   A2[,r] = sum_k Vinv_k[,r] * xsum_kr   where xsum_kr is the
//   per-pattern centered/scaled column sum
//   Xbar[j,,] = Vinvsuminv %*% A2
//
// This is O(J * K * R^2 + J * R^3) instead of O(J * N * R^2).
//
// Vinv_list:   list of K  R x R  precision matrices
// pattern:     K x R integer (0/1, observed conditions per pattern)
// raw_sum:     K x J matrix
// n_k:         K-vector
// cm:          J x R per-condition column means
// csd:         J x R per-condition column SDs
// Vinvsuminv:  R x R inverse of sum_k n_k * Vinv_k
//
// Returns J x R x R array (Xbar).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::cube compute_Xbar_from_sums_rcpp(
    const Rcpp::List& Vinv_list,
    const arma::imat& pattern,
    const arma::mat& raw_sum,
    const arma::ivec& n_k,
    const arma::mat& cm,
    const arma::mat& csd,
    const arma::mat& Vinvsuminv) {

  unsigned int K = pattern.n_rows;
  unsigned int R = pattern.n_cols;
  unsigned int J = raw_sum.n_cols;

  // Cache Vinv matrices and per-pattern observed indices
  std::vector<arma::mat> Vinv(K);
  std::vector<std::vector<unsigned int>> obs_r_list(K);
  for (unsigned int k = 0; k < K; k++) {
    Vinv[k] = Rcpp::as<arma::mat>(Vinv_list[k]);
    for (unsigned int r = 0; r < R; r++) {
      if (pattern(k, r)) obs_r_list[k].push_back(r);
    }
  }

  arma::cube Xbar(J, R, R, fill::zeros);

  for (unsigned int j = 0; j < J; j++) {
    // A2 = sum_k sum_{r in obs_k} xsum_r * Vinv_k[,r]
    arma::mat A2(R, R, fill::zeros);
    for (unsigned int k = 0; k < K; k++) {
      for (unsigned int ri = 0; ri < obs_r_list[k].size(); ri++) {
        unsigned int r = obs_r_list[k][ri];
        double xsum_r = (raw_sum(k, j) - n_k(k) * cm(j, r)) / csd(j, r);
        A2.col(r) += Vinv[k].col(r) * xsum_r;
      }
    }

    // Xbar[j,,] = Vinvsuminv %*% A2
    arma::mat Xbar_j = Vinvsuminv * A2;
    for (unsigned int r1 = 0; r1 < R; r1++) {
      for (unsigned int r2 = 0; r2 < R; r2++) {
        Xbar(j, r1, r2) = Xbar_j(r1, r2);
      }
    }
  }

  return Xbar;
}


// ---------------------------------------------------------------------------
// 8. compute_bxxb_correction_3d_rcpp
//
// Computes the bxxb posterior correction for residual variance estimation
// in the 3D missing-data path.  This is the first term of the posterior
// correction: E_q[(X_3d b_l)' (X_3d b_l)] summed over all L effects.
//
// For each effect l and variable j with alpha_lj > threshold:
//   G_j[r,s] = sum_k cross_rs(k,j,r,s)    (per-outcome Gram matrix)
//   bxxb += alpha_lj * (G_j . mu_outer_lj) (Hadamard product)
//   bxxb[r,r] += alpha_lj * G_j[r,r] * post_var_lj[r]  (diagonal correction)
//
// where cross_rs uses the same formula as compute_svs_inv_3d_rcpp, and
// post_var_lj[r] = mu2_diag[l,j,r] - mu[l,j,r]^2.
//
// alpha:       L x J
// mu:          L x J x R  (posterior means)
// mu2_diag:    L x J x R  (diagonal posterior second moments)
// pattern:     K x R  (0/1 missingness patterns)
// raw_sq_sum:  K x J
// raw_sum:     K x J
// n_k:         K
// cm:          J x R  (per-outcome column means)
// csd:         J x R  (per-outcome column SDs)
//
// Returns R x R matrix: total bxxb summed over all L effects.
// The caller subtracts the b1_XtX_b1 terms in R.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat compute_bxxb_correction_3d_rcpp(
    const arma::mat& alpha,
    const arma::cube& mu,
    const arma::cube& mu2_diag,
    const arma::imat& pattern,
    const arma::mat& raw_sq_sum,
    const arma::mat& raw_sum,
    const arma::ivec& n_k,
    const arma::mat& cm,
    const arma::mat& csd) {

  unsigned int L = alpha.n_rows;
  unsigned int J = alpha.n_cols;
  unsigned int K = pattern.n_rows;
  unsigned int R = pattern.n_cols;

  // Cache per-pattern observed indices
  std::vector<std::vector<unsigned int>> obs_r_list(K);
  for (unsigned int k = 0; k < K; k++) {
    for (unsigned int r = 0; r < R; r++) {
      if (pattern(k, r)) obs_r_list[k].push_back(r);
    }
  }

  arma::mat total_bxxb(R, R, arma::fill::zeros);

  for (unsigned int l = 0; l < L; l++) {
    arma::mat bxxb_l(R, R, arma::fill::zeros);

    for (unsigned int j = 0; j < J; j++) {
      if (alpha(l, j) < 1e-10) continue;

      // Compute G_j: per-outcome Gram matrix from pattern summary stats
      arma::mat G_j(R, R, arma::fill::zeros);
      for (unsigned int k = 0; k < K; k++) {
        const auto& obs_r = obs_r_list[k];
        unsigned int n_obs = obs_r.size();
        if (n_obs == 0) continue;

        for (unsigned int ri = 0; ri < n_obs; ri++) {
          unsigned int r = obs_r[ri];
          for (unsigned int si = ri; si < n_obs; si++) {
            unsigned int s = obs_r[si];
            double cross_rs = (raw_sq_sum(k, j) -
              (cm(j, r) + cm(j, s)) * raw_sum(k, j) +
              n_k(k) * cm(j, r) * cm(j, s)) /
              (csd(j, r) * csd(j, s));
            G_j(r, s) += cross_rs;
            if (r != s) G_j(s, r) += cross_rs;
          }
        }
      }

      // Hadamard product: G_j . (mu_lj * mu_lj')
      // Plus diagonal correction for exact second moment
      double a_lj = alpha(l, j);
      for (unsigned int r = 0; r < R; r++) {
        for (unsigned int s = r; s < R; s++) {
          double mu_prod = mu(l, j, r) * mu(l, j, s);
          double val = a_lj * G_j(r, s) * mu_prod;
          bxxb_l(r, s) += val;
          if (r != s) bxxb_l(s, r) += val;
        }
        // Diagonal: replace mu^2 with exact mu2_diag
        // We already added a_lj * G_j(r,r) * mu(l,j,r)^2 above,
        // now add the posterior variance part:
        //   a_lj * G_j(r,r) * (mu2_diag(l,j,r) - mu(l,j,r)^2)
        double post_var_r = mu2_diag(l, j, r) - mu(l, j, r) * mu(l, j, r);
        bxxb_l(r, r) += a_lj * G_j(r, r) * post_var_r;
      }
    }

    total_bxxb += bxxb_l;
  }

  return total_bxxb;
}
