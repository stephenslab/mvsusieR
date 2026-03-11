// Inner EM loop for mixture weight updates.
//
// Solves:  max_{pi >= 0, sum(pi)=1}  sum_i w_i * log(sum_k pi_k * L_{ik})
// where L_{ik} = exp(llik[i,k]) and w_i are alpha weights from SuSiE.
//
// E-step: phi_{ik} = pi_k * L_{ik} / sum_{k'} pi_{k'} * L_{ik'}
// M-step: pi_k^{new} = sum_i w_i * phi_{ik} / sum_i w_i

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List inner_em_rcpp(const arma::mat& llik,      // (N x K) log-likelihood matrix
                  const arma::vec& weights,    // (N) observation weights (alpha)
                  arma::vec pi_init,           // (K) initial mixture weights
                  int max_iter,                // max inner iterations
                  double tol) {                // convergence tolerance

  int N = llik.n_rows;
  int K = llik.n_cols;

  arma::vec pi_cur = pi_init;
  arma::vec pi_new(K);
  arma::vec log_pi(K);
  arma::mat log_phi(N, K);
  arma::mat phi(N, K);
  arma::vec row_max(N);
  int n_iter = 0;

  // Precompute total weight (constant across iterations)
  double w_total = arma::accu(weights);

  for (int t = 0; t < max_iter; t++) {
    n_iter = t + 1;

    // --- E-step ---
    // log_phi[i,k] = llik[i,k] + log(pi_k)
    log_pi = arma::log(arma::clamp(pi_cur, 1e-300, arma::datum::inf));
    for (int k = 0; k < K; k++) {
      log_phi.col(k) = llik.col(k) + log_pi(k);
    }

    // Log-sum-exp normalization for numerical stability
    row_max = arma::max(log_phi, 1);
    for (int k = 0; k < K; k++) {
      log_phi.col(k) -= row_max;
    }
    phi = arma::exp(log_phi);

    // Normalize rows: phi[i,] /= sum(phi[i,])
    arma::vec row_sums = arma::sum(phi, 1);
    for (int k = 0; k < K; k++) {
      phi.col(k) /= row_sums;
    }

    // --- M-step ---
    // pi_new[k] = sum(weights * phi[,k]) / sum(weights)
    for (int k = 0; k < K; k++) {
      pi_new(k) = arma::dot(weights, phi.col(k)) / w_total;
    }

    // Check convergence: max absolute change in pi
    double max_diff = arma::max(arma::abs(pi_new - pi_cur));
    pi_cur = pi_new;
    if (max_diff < tol) break;
  }

  return List::create(
    Named("pi") = pi_cur,
    Named("n_iter") = n_iter,
    Named("converged") = (n_iter < max_iter)
  );
}
