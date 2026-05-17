// Rcpp/Armadillo implementation of missing data imputation with
// pattern-based caching. Precomputes Cholesky decomposition once per
// unique missingness pattern, then applies vectorized imputation to
// all observations sharing that pattern.
//
// optimization: pattern-based grouping avoids redundant O(|M|^3)
// Cholesky per observation.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Log-determinant from upper Cholesky factor: log|A| = 2 sum(log(diag(chol(A))))
static double chol2ldet_cpp(const mat& R) {
  return 2.0 * accu(log(R.diag()));
}

// Pattern-based imputation with precomputed Cholesky per unique pattern.
//
// Y:           N x R matrix (will be modified in-place for imputed entries)
// mu:          N x R matrix of fitted values (model$Xr)
// Vinv:        R x R precision matrix (V^{-1})
// patterns:    list of logical vectors (length R each), one per unique pattern
// pattern_obs: list of integer vectors (1-indexed row indices), one per pattern
//
// Returns list with:
//   Y:                  N x R imputed matrix
//   Y_cov:              R x R accumulated imputation covariance
//   sum_neg_ent_Y_miss: scalar negative entropy sum
//
// [[Rcpp::export]]
List impute_missing_Y_rcpp(arma::mat Y, const arma::mat& mu,
                           const arma::mat& Vinv,
                           const List& patterns,
                           const List& pattern_obs) {
  unsigned int r = Y.n_cols;
  mat Y_cov(r, r, fill::zeros);
  double sum_neg_ent = 0.0;
  int n_patterns = patterns.size();

  for (int k = 0; k < n_patterns; k++) {
    LogicalVector mask_r = patterns[k];
    IntegerVector obs_idx_r = pattern_obs[k];

    // Build index vectors for missing and observed conditions
    uvec miss_idx, obs_idx_cols;
    for (int j = 0; j < mask_r.size(); j++) {
      if (mask_r[j]) {
        miss_idx.resize(miss_idx.n_elem + 1);
        miss_idx(miss_idx.n_elem - 1) = (unsigned)j;
      } else {
        obs_idx_cols.resize(obs_idx_cols.n_elem + 1);
        obs_idx_cols(obs_idx_cols.n_elem - 1) = (unsigned)j;
      }
    }

    // Convert R 1-indexed to C++ 0-indexed row indices
    uvec row_idx(obs_idx_r.size());
    for (int i = 0; i < obs_idx_r.size(); i++) {
      row_idx(i) = obs_idx_r[i] - 1;
    }

    unsigned int n_miss = miss_idx.n_elem;
    unsigned int n_obs_k = row_idx.n_elem;

    if (n_miss == 0 || n_obs_k == 0) continue;

    // Precompute for this pattern: O(|M|^3) Cholesky, done ONCE
    mat Vinv_mm = Vinv.submat(miss_idx, miss_idx);
    mat Vinv_mm_chol = chol(Vinv_mm);  // upper triangular
    mat cov_mm = solve(trimatu(Vinv_mm_chol),
                       solve(trimatl(trans(Vinv_mm_chol)),
                             eye(n_miss, n_miss)));

    // Accumulate Y_cov: n_obs_k copies of Lambda_{M,M}^{-1}
    Y_cov.submat(miss_idx, miss_idx) += n_obs_k * cov_mm;

    // Accumulate negative entropy
    sum_neg_ent += n_obs_k * 0.5 *
      (n_miss * log(1.0 / (2.0 * M_PI * exp(1.0))) +
       chol2ldet_cpp(Vinv_mm_chol));

    // Vectorized imputation for all observations with this pattern
    if (obs_idx_cols.n_elem > 0) {
      mat Vinv_mo = Vinv.submat(miss_idx, obs_idx_cols);
      // resid = Y[rows, obs_cols] - mu[rows, obs_cols]  (n_obs_k x |O|)
      mat resid = Y.submat(row_idx, obs_idx_cols) -
                  mu.submat(row_idx, obs_idx_cols);
      // correction = resid * Vinv_mo' * cov_mm  (n_obs_k x |M|)
      // Note: cov_mm is symmetric so cov_mm' = cov_mm
      mat correction = resid * trans(Vinv_mo) * cov_mm;
      // Y[rows, miss] = mu[rows, miss] - correction
      Y.submat(row_idx, miss_idx) =
        mu.submat(row_idx, miss_idx) - correction;
    } else {
      // All conditions missing: imputed value = mu (no conditioning)
      Y.submat(row_idx, miss_idx) = mu.submat(row_idx, miss_idx);
    }
  }

  return List::create(Named("Y") = Y,
                      Named("Y_cov") = Y_cov,
                      Named("sum_neg_ent_Y_miss") = sum_neg_ent);
}
