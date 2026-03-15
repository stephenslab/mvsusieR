// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Helper: safe Cholesky decomposition with ridge fallback.
// Returns false only if all attempts fail.
// ---------------------------------------------------------------------------
static bool safe_chol_lower(arma::mat& L, const arma::mat& A) {
  if (arma::chol(L, A, "lower")) return true;
  unsigned int R = A.n_rows;
  if (arma::chol(L, A + 1e-10 * arma::eye(R, R), "lower")) return true;
  if (arma::chol(L, A + 1e-6  * arma::eye(R, R), "lower")) return true;
  return false;
}


// ---------------------------------------------------------------------------
// 1. precompute_eigen_cache_non_common_rcpp
//
// Batch eigendecomposition for the non-common-covariance path.
// For each variable j, computes chol(SVS_j) once and reuses it for
// all K prior structure matrices.
//
// Key optimization: the Cholesky of SVS_j is computed once per j
// (not K times), saving (K-1)*J Cholesky decompositions.
//
// svs_list:  list of J  R x R  matrices (per-variable sampling variance)
// U_list:    list of K  R x R  matrices (prior structure)
//
// Returns list with:
//   is_common_cov = FALSE
//   log_det_svs = J-vector
//   components = list of K items, each with:
//     Q = R x R x J cube
//     G = R x R x J cube
//     eigenvalues = R x J matrix
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List precompute_eigen_cache_non_common_rcpp(
    const Rcpp::List& svs_list,
    const Rcpp::List& U_list) {

  int J = svs_list.size();
  int K = U_list.size();

  // Cache U matrices
  std::vector<arma::mat> U(K);
  for (int k = 0; k < K; k++) {
    U[k] = Rcpp::as<arma::mat>(U_list[k]);
  }
  unsigned int R = U[0].n_rows;

  arma::vec log_det_svs(J);

  // Pre-allocate output arrays for each component
  std::vector<arma::cube> Q_all(K);
  std::vector<arma::cube> G_all(K);
  std::vector<arma::mat>  eig_all(K);
  for (int k = 0; k < K; k++) {
    Q_all[k]   = arma::cube(R, R, J);
    G_all[k]   = arma::cube(R, R, J);
    eig_all[k] = arma::mat(R, J);
  }

  arma::mat I_R = arma::eye(R, R);

  // Main loop: j outer (Cholesky cached), k inner
  for (int j = 0; j < J; j++) {
    arma::mat SVS = Rcpp::as<arma::mat>(svs_list[j]);

    // Cholesky: L L' = SVS (computed once per j)
    arma::mat L;
    if (!safe_chol_lower(L, SVS)) {
      // Last resort: identity (degenerate, but prevents crash)
      L = I_R;
    }
    log_det_svs(j) = 2.0 * arma::sum(arma::log(L.diag()));

    // L^{-1} via triangular solve
    arma::mat L_inv = arma::solve(arma::trimatl(L), I_R);
    arma::mat L_inv_t = L_inv.t();

    for (int k = 0; k < K; k++) {
      // M = L^{-1} U_k L^{-T}
      arma::mat M = L_inv * U[k] * L_inv_t;
      M = (M + M.t()) / 2.0;  // symmetrize

      // Eigendecompose (ascending order from eig_sym)
      arma::vec d;
      arma::mat P;
      arma::eig_sym(d, P, M);

      // Reverse to descending order (match R's eigen())
      d = arma::reverse(d);
      P = arma::fliplr(P);

      // Clamp negative eigenvalues to 0
      d = arma::clamp(d, 0.0, arma::datum::inf);

      // Q = L^{-T} P,  G = L P
      Q_all[k].slice(j) = L_inv_t * P;
      G_all[k].slice(j) = L * P;
      eig_all[k].col(j) = d;
    }
  }

  // Package output as list of K components
  Rcpp::List components(K);
  for (int k = 0; k < K; k++) {
    components[k] = Rcpp::List::create(
      Rcpp::Named("Q") = Q_all[k],
      Rcpp::Named("G") = G_all[k],
      Rcpp::Named("eigenvalues") = eig_all[k]
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("is_common_cov") = false,
    Rcpp::Named("log_det_svs") = log_det_svs,
    Rcpp::Named("components") = components
  );
}


// ---------------------------------------------------------------------------
// 2. loglik_non_common_rcpp
//
// Compute J x (K+1) log-likelihood matrix for the non-common-cov path.
// Column 1 = null component N(0, SVS).
// Columns 2..(K+1) = non-null components N(0, SVS + V*U_k).
//
// betahat:      J x R matrix
// V_scalar:     scalar prior variance multiplier
// log_det_svs:  J-vector of log|SVS_j|
// components:   list of K items (Q, G, eigenvalues)
//
// Returns J x (K+1) matrix.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat loglik_non_common_rcpp(
    const arma::mat& betahat,
    double V_scalar,
    const arma::vec& log_det_svs,
    const Rcpp::List& components) {

  unsigned int J = betahat.n_rows;
  unsigned int R = betahat.n_cols;
  int K = components.size();
  double const_term = -(double)R / 2.0 * std::log(2.0 * M_PI);

  arma::mat llik(J, K + 1, fill::zeros);

  // Cache component data
  std::vector<arma::cube> Q_all(K);
  std::vector<arma::mat>  eig_all(K);
  for (int k = 0; k < K; k++) {
    Rcpp::List comp_k = components[k];
    Q_all[k]   = Rcpp::as<arma::cube>(comp_k["Q"]);
    eig_all[k] = Rcpp::as<arma::mat>(comp_k["eigenvalues"]);
  }

  // Null component: SVS^{-1} = Q Q'  (using Q from first component)
  for (unsigned int j = 0; j < J; j++) {
    arma::vec b_rot = Q_all[0].slice(j).t() * betahat.row(j).t();
    llik(j, 0) = const_term - 0.5 * log_det_svs(j) - 0.5 * arma::dot(b_rot, b_rot);
  }

  // Non-null components
  for (int k = 0; k < K; k++) {
    for (unsigned int j = 0; j < J; j++) {
      arma::vec d_kj = eig_all[k].col(j);
      arma::vec b_rot = Q_all[k].slice(j).t() * betahat.row(j).t();

      double log_det = log_det_svs(j) + arma::sum(arma::log1p(V_scalar * d_kj));
      arma::vec inv_factors = 1.0 / (1.0 + V_scalar * d_kj);
      double mahal = arma::dot(b_rot % b_rot, inv_factors);

      llik(j, k + 1) = const_term - 0.5 * log_det - 0.5 * mahal;
    }
  }

  return llik;
}


// ---------------------------------------------------------------------------
// 3. posterior_non_common_rcpp
//
// Compute posterior moments for the non-common-cov path.
// For each (k, j): computes posterior mean, second moment, sign
// probabilities, and EM statistics.
//
// betahat:     J x R matrix
// V_scalar:    scalar prior variance multiplier
// components:  list of K items (Q, G, eigenvalues)
// pi_V_post:   J x (K+1) posterior mixture weights
// em_var_wt:   (K+1) x J matrix (or 0x0 if EM not needed)
//
// Returns list(mu, mu2, post_neg, post_zero,
//              prior_scale_em_update)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List posterior_non_common_rcpp(
    const arma::mat& betahat,
    double V_scalar,
    const Rcpp::List& components,
    const arma::mat& pi_V_post,
    const arma::mat& em_var_wt) {

  unsigned int J = betahat.n_rows;
  unsigned int R = betahat.n_cols;
  int K = components.size();
  bool do_em = (em_var_wt.n_rows > 0);

  arma::mat mu(J, R, fill::zeros);
  // Use R x R x J internally for efficient slice access
  arma::cube pm2_internal(R, R, J, fill::zeros);
  arma::mat post_neg(J, R, fill::zeros);
  arma::mat post_zero(J, R);
  arma::vec em_update(K + 1, fill::zeros);

  // Null component: beta = 0 (contributes to post_zero)
  for (unsigned int j = 0; j < J; j++) {
    for (unsigned int r = 0; r < R; r++) {
      post_zero(j, r) = pi_V_post(j, 0);
    }
  }

  // Cache component data
  std::vector<arma::cube> Q_all(K);
  std::vector<arma::cube> G_all(K);
  std::vector<arma::mat>  eig_all(K);
  for (int k = 0; k < K; k++) {
    Rcpp::List comp_k = components[k];
    Q_all[k]   = Rcpp::as<arma::cube>(comp_k["Q"]);
    G_all[k]   = Rcpp::as<arma::cube>(comp_k["G"]);
    eig_all[k] = Rcpp::as<arma::mat>(comp_k["eigenvalues"]);
  }

  double eps_thresh = std::sqrt(std::numeric_limits<double>::epsilon());

  for (int k = 0; k < K; k++) {
    for (unsigned int j = 0; j < J; j++) {
      arma::mat Q_j = Q_all[k].slice(j);
      arma::mat G_j = G_all[k].slice(j);
      arma::vec d_kj = eig_all[k].col(j);
      double w = pi_V_post(j, k + 1);

      // Shrinkage factors
      arma::vec shrink = V_scalar * d_kj / (1.0 + V_scalar * d_kj);
      arma::vec inv_factor = 1.0 / (1.0 + V_scalar * d_kj);

      // C_k = G diag(shrink) G'  (R x R)
      arma::mat G_shrunk = G_j * arma::diagmat(shrink);
      arma::mat C_k = G_shrunk * G_j.t();
      C_k = (C_k + C_k.t()) / 2.0;

      // Posterior mean: m_j = G (shrink .* Q' bhat_j)
      arma::vec b_rot = Q_j.t() * betahat.row(j).t();
      arma::vec m_j = G_j * (shrink % b_rot);

      // Accumulate posterior mean
      mu.row(j) += w * m_j.t();

      // Accumulate posterior second moment: C_k + m_j m_j'
      arma::mat update = w * (C_k + m_j * m_j.t());
      pm2_internal.slice(j) += update;

      // LFSR: sign probabilities
      arma::vec diag_Ck = C_k.diag();
      double max_diag = arma::max(diag_Ck);
      double thresh = eps_thresh * max_diag;
      for (unsigned int r = 0; r < R; r++) {
        if (diag_Ck(r) > thresh && diag_Ck(r) > 0) {
          double sd_r = std::sqrt(diag_Ck(r));
          post_neg(j, r) += w * R::pnorm(0.0, m_j(r), sd_r, 1, 0);
        } else {
          // Point mass at 0 for this component
          post_zero(j, r) += w;
        }
      }

      // EM update: sum inv_factor only over positive eigenvalues (non-null
      // subspace of the prior).  Null-space eigenvalues (d=0) give
      // inv_factor=1, which would inflate the trace by V per null dimension.
      arma::vec d_pos = arma::conv_to<arma::vec>::from(
          d_kj > std::max(d_kj.max() * 1e-8, 0.0));
      double tr_term = V_scalar * arma::dot(inv_factor, d_pos);
      double em_j = V_scalar * V_scalar *
        arma::dot(d_kj % inv_factor % inv_factor, b_rot % b_rot);
      double em_wt_j = do_em ? em_var_wt(k + 1, j) : w;
      em_update(k + 1) += em_wt_j * (tr_term + em_j);
    }
  }

  // Convert pm2 from R x R x J to J x R x R for R
  arma::cube mu2(J, R, R);
  for (unsigned int j = 0; j < J; j++) {
    for (unsigned int r1 = 0; r1 < R; r1++) {
      for (unsigned int r2 = 0; r2 < R; r2++) {
        mu2(j, r1, r2) = pm2_internal(r1, r2, j);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mu") = mu,
    Rcpp::Named("mu2") = mu2,
    Rcpp::Named("post_neg") = post_neg,
    Rcpp::Named("post_zero") = post_zero,
    Rcpp::Named("prior_scale_em_update") = em_update
  );
}


// ---------------------------------------------------------------------------
// 4. accumulate_mu2_common_rcpp
//
// Helper for the common-cov path: replaces the R-level J-loop for
// accumulating mu2 (posterior second moment).
//
// mu2:  J x R x R array (current values, modified in place)
// M_k:  J x R matrix of posterior means for component k
// C_k:  R x R posterior covariance (same for all j in common-cov)
// w_k:  J-vector of mixture weights for component k
//
// Returns updated J x R x R array.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::cube accumulate_mu2_common_rcpp(
    arma::cube mu2,
    const arma::mat& M_k,
    const arma::mat& C_k,
    const arma::vec& w_k) {

  unsigned int J = M_k.n_rows;
  unsigned int R = M_k.n_cols;

  for (unsigned int j = 0; j < J; j++) {
    arma::vec m = M_k.row(j).t();
    arma::mat update = w_k(j) * (C_k + m * m.t());
    for (unsigned int r1 = 0; r1 < R; r1++) {
      for (unsigned int r2 = 0; r2 < R; r2++) {
        mu2(j, r1, r2) += update(r1, r2);
      }
    }
  }

  return mu2;
}


// ---------------------------------------------------------------------------
// 5. loglik_common_rcpp
//
// C++ implementation of the common-cov loglik K-loop.
// All variables share the same SVS, so Q_k are R x R matrices (not cubes).
// Accepts optional BQ_cache to skip recomputing betahat %*% Q_k.
// Returns both the J x (K+1) log-likelihood matrix and BQ_cache.
//
// betahat:       J x R matrix
// V_scalar:      scalar prior variance multiplier
// log_det_svs:   scalar log|SVS|
// components:    list of K items, each with Q (R x R), eigenvalues (R-vector)
// BQ_cache_in:   list of K J x R matrices (or empty list to compute fresh)
//
// Returns list(llik = J x (K+1), BQ_cache = list of K J x R matrices)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List loglik_common_rcpp(
    const arma::mat& betahat,
    double V_scalar,
    double log_det_svs,
    const Rcpp::List& components,
    const Rcpp::List& BQ_cache_in) {

  unsigned int J = betahat.n_rows;
  unsigned int R = betahat.n_cols;
  int K = components.size();
  bool has_cache = (BQ_cache_in.size() == K);
  double const_term = -(double)R / 2.0 * std::log(2.0 * M_PI);

  arma::mat llik(J, K + 1, arma::fill::zeros);

  // Load or compute BQ products
  std::vector<arma::vec> d_all(K);
  std::vector<arma::mat> BQ_all(K);
  for (int k = 0; k < K; k++) {
    Rcpp::List comp_k = components[k];
    d_all[k] = Rcpp::as<arma::vec>(comp_k["eigenvalues"]);
    if (has_cache) {
      BQ_all[k] = Rcpp::as<arma::mat>(BQ_cache_in[k]);
    } else {
      arma::mat Q_k = Rcpp::as<arma::mat>(comp_k["Q"]);
      BQ_all[k] = betahat * Q_k;
    }
  }

  // Null component: uses BQ from component 1
  arma::vec mahal_null = arma::sum(BQ_all[0] % BQ_all[0], 1);
  llik.col(0) = const_term - 0.5 * log_det_svs - 0.5 * mahal_null;

  // Non-null components
  for (int k = 0; k < K; k++) {
    double log_det = log_det_svs + arma::sum(arma::log1p(V_scalar * d_all[k]));
    arma::vec inv_factors = 1.0 / (1.0 + V_scalar * d_all[k]);
    arma::vec mahal = (BQ_all[k] % BQ_all[k]) * inv_factors;
    llik.col(k + 1) = const_term - 0.5 * log_det - 0.5 * mahal;
  }

  // Return BQ_cache for reuse
  Rcpp::List BQ_cache_out(K);
  if (!has_cache) {
    for (int k = 0; k < K; k++)
      BQ_cache_out[k] = BQ_all[k];
  } else {
    BQ_cache_out = BQ_cache_in;  // pass through
  }

  return Rcpp::List::create(
    Rcpp::Named("llik") = llik,
    Rcpp::Named("BQ_cache") = BQ_cache_out
  );
}


// ---------------------------------------------------------------------------
// 6. posterior_common_rcpp
//
// C++ implementation of the common-cov posterior K-loop.
// Computes posterior means, sign probabilities (LFSR), EM statistics,
// and optionally reduced second-moment statistics (bxxb, alpha_mu2_sum,
// mu2_diag) to avoid allocating a J x R x R array.
//
// betahat:      J x R matrix
// V_scalar:     scalar prior variance multiplier
// components:   list of K items (Q, G, eigenvalues -- all R x R / R-vector)
// pi_V_post:    J x (K+1) posterior mixture weights
// em_var_wt:    (K+1) x J matrix (or 0x0 if EM not needed)
// BQ_cache:     list of K J x R matrices (or empty list to compute fresh)
// do_reduce:    if true, compute bxxb/alpha_mu2_sum/mu2_diag instead of mu2
// alpha:        J-vector (only used if do_reduce)
// d_var:        J-vector of d values (only used if do_reduce)
// v_inv:        R x R matrix (only used if do_reduce, for vbxxb)
//
// Returns list(mu, post_neg, post_zero, prior_scale_em_update,
//              and either mu2 or reduced stats)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List posterior_common_rcpp(
    const arma::mat& betahat,
    double V_scalar,
    const Rcpp::List& components,
    const arma::mat& pi_V_post,
    const arma::mat& em_var_wt,
    const Rcpp::List& BQ_cache_list,
    bool do_reduce,
    const arma::vec& alpha,
    const arma::vec& d_var,
    const arma::mat& v_inv) {

  unsigned int J = betahat.n_rows;
  unsigned int R = betahat.n_cols;
  int K = components.size();
  bool do_em = (em_var_wt.n_rows > 0);
  bool has_BQ_cache = (BQ_cache_list.size() == K);

  arma::mat mu(J, R, arma::fill::zeros);
  arma::mat post_neg(J, R, arma::fill::zeros);
  arma::mat post_zero(J, R);
  arma::vec em_update(K + 1, arma::fill::zeros);

  // Null component: beta = 0
  for (unsigned int j = 0; j < J; j++)
    for (unsigned int r = 0; r < R; r++)
      post_zero(j, r) = pi_V_post(j, 0);

  // Reduce accumulators
  arma::mat bxxb(R, R, arma::fill::zeros);
  arma::mat alpha_mu2_sum(R, R, arma::fill::zeros);
  arma::mat mu2_diag(J, R, arma::fill::zeros);

  // Full mu2 (only if !do_reduce)
  arma::cube pm2_internal;
  if (!do_reduce)
    pm2_internal = arma::cube(R, R, J, arma::fill::zeros);

  // Cache component data
  std::vector<arma::mat> Q_all(K), G_all(K);
  std::vector<arma::vec> d_all(K);
  for (int k = 0; k < K; k++) {
    Rcpp::List comp_k = components[k];
    Q_all[k] = Rcpp::as<arma::mat>(comp_k["Q"]);
    G_all[k] = Rcpp::as<arma::mat>(comp_k["G"]);
    d_all[k] = Rcpp::as<arma::vec>(comp_k["eigenvalues"]);
  }

  // Load BQ cache if available
  std::vector<arma::mat> BQ_all(K);
  if (has_BQ_cache) {
    for (int k = 0; k < K; k++)
      BQ_all[k] = Rcpp::as<arma::mat>(BQ_cache_list[k]);
  }

  double eps_thresh = std::sqrt(std::numeric_limits<double>::epsilon());

  for (int k = 0; k < K; k++) {
    arma::vec d_k = d_all[k];

    // Shrinkage factors
    arma::vec shrink = V_scalar * d_k / (1.0 + V_scalar * d_k);
    arma::vec inv_factor = 1.0 / (1.0 + V_scalar * d_k);

    // Posterior covariance (same for all j in common-cov)
    arma::mat G_scaled = G_all[k] * arma::diagmat(shrink);
    arma::mat C_k = G_scaled * G_all[k].t();
    C_k = (C_k + C_k.t()) / 2.0;

    // BQ: betahat %*% Q_k (J x R)
    arma::mat BQ = has_BQ_cache ? BQ_all[k] : betahat * Q_all[k];

    // Posterior means: M_k = BQ * diag(shrink) * G_k' (J x R)
    arma::mat BQ_shrunk = BQ * arma::diagmat(shrink);
    arma::mat M_k = BQ_shrunk * G_all[k].t();

    // Accumulate weighted posterior mean
    for (unsigned int j = 0; j < J; j++) {
      double w = pi_V_post(j, k + 1);
      mu.row(j) += w * M_k.row(j);
    }

    if (do_reduce) {
      // Accumulate bxxb, alpha_mu2_sum, mu2_diag
      arma::vec aw(J), daw(J);
      for (unsigned int j = 0; j < J; j++) {
        aw(j) = alpha(j) * pi_V_post(j, k + 1);
        daw(j) = d_var(j) * aw(j);
      }

      double sum_daw = arma::sum(daw);
      double sum_aw = arma::sum(aw);
      bxxb += sum_daw * C_k;
      alpha_mu2_sum += sum_aw * C_k;

      // Outer product contributions: sum_j wt_j * M_kj M_kj'
      arma::mat sqrt_daw_M = M_k;
      arma::mat sqrt_aw_M = M_k;
      for (unsigned int j = 0; j < J; j++) {
        sqrt_daw_M.row(j) *= std::sqrt(daw(j));
        sqrt_aw_M.row(j) *= std::sqrt(aw(j));
      }
      bxxb += sqrt_daw_M.t() * sqrt_daw_M;
      alpha_mu2_sum += sqrt_aw_M.t() * sqrt_aw_M;

      // mu2_diag
      arma::vec diag_Ck = C_k.diag();
      for (unsigned int j = 0; j < J; j++) {
        double w = pi_V_post(j, k + 1);
        for (unsigned int r = 0; r < R; r++) {
          mu2_diag(j, r) += w * (diag_Ck(r) + M_k(j, r) * M_k(j, r));
        }
      }
    } else {
      // Full mu2
      for (unsigned int j = 0; j < J; j++) {
        double w = pi_V_post(j, k + 1);
        arma::vec m = M_k.row(j).t();
        pm2_internal.slice(j) += w * (C_k + m * m.t());
      }
    }

    // Sign probabilities for LFSR
    arma::vec diag_Ck = C_k.diag();
    double max_diag = arma::max(diag_Ck);
    double thresh = eps_thresh * max_diag;
    for (unsigned int r = 0; r < R; r++) {
      if (diag_Ck(r) > thresh && diag_Ck(r) > 0) {
        double sd_r = std::sqrt(diag_Ck(r));
        for (unsigned int j = 0; j < J; j++) {
          double w = pi_V_post(j, k + 1);
          post_neg(j, r) += w * R::pnorm(0.0, M_k(j, r), sd_r, 1, 0);
        }
      } else {
        for (unsigned int j = 0; j < J; j++) {
          post_zero(j, r) += pi_V_post(j, k + 1);
        }
      }
    }

    // EM statistic: sum inv_factor only over positive eigenvalues (non-null
    // subspace of the prior).  Null-space eigenvalues (d=0) give inv_factor=1,
    // which would inflate the trace by V per null dimension.
    arma::vec d_pos = arma::conv_to<arma::vec>::from(
        d_k > std::max(d_k.max() * 1e-8, 0.0));
    double tr_term = V_scalar * arma::dot(inv_factor, d_pos);
    arma::vec em_per_var = V_scalar * V_scalar *
      ((BQ % BQ) * (d_k % inv_factor % inv_factor));
    for (unsigned int j = 0; j < J; j++) {
      double em_wt_j = do_em ? em_var_wt(k + 1, j) : pi_V_post(j, k + 1);
      em_update(k + 1) += em_wt_j * (tr_term + em_per_var(j));
    }
  }

  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("mu") = mu,
    Rcpp::Named("post_neg") = post_neg,
    Rcpp::Named("post_zero") = post_zero,
    Rcpp::Named("prior_scale_em_update") = em_update
  );

  if (do_reduce) {
    double vbxxb = arma::accu(v_inv % bxxb);
    result["bxxb"] = bxxb;
    result["vbxxb"] = vbxxb;
    result["alpha_mu2_sum"] = alpha_mu2_sum;
    result["mu2_diag"] = mu2_diag;
  } else {
    // Convert pm2 from R x R x J to J x R x R for R
    arma::cube mu2(J, R, R);
    for (unsigned int j = 0; j < J; j++)
      for (unsigned int r1 = 0; r1 < R; r1++)
        for (unsigned int r2 = 0; r2 < R; r2++)
          mu2(j, r1, r2) = pm2_internal(r1, r2, j);
    result["mu2"] = mu2;
  }

  return result;
}
