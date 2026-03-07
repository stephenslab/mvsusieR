context("Test eigendecomposition precomputation optimization")

# =============================================================================
# Unit test: eigendecompose_one_pair properties
# =============================================================================

test_that("eigendecompose_one_pair satisfies algebraic identities", {
  set.seed(1)
  R <- 4

  # Generate random SPD matrix SVS and PSD matrix U
  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.1
  B <- matrix(rnorm(R * R), R, R)
  U <- crossprod(B)

  decomp <- mvsusieR:::eigendecompose_one_pair(SVS, U)
  Q <- decomp$Q
  G <- decomp$G
  d <- decomp$eigenvalues

  # Q'G = I
  expect_equal(crossprod(Q, G), diag(R), tolerance = 1e-12)

  # G diag(d) G' = U
  expect_equal(G %*% diag(d) %*% t(G), U, tolerance = 1e-12)

  # Q Q' = SVS^{-1}  (note: Q Q', not Q'Q)
  expect_equal(tcrossprod(Q), solve(SVS), tolerance = 1e-12)

  # log_det correct
  expect_equal(decomp$log_det_svs, log(det(SVS)), tolerance = 1e-12)

  # Eigenvalues are non-negative
  expect_true(all(d >= 0))
})

test_that("eigendecompose_one_pair works with singular U", {
  set.seed(2)
  R <- 3
  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.1

  # Rank-1 U
  v <- rnorm(R)
  U <- tcrossprod(v)

  decomp <- mvsusieR:::eigendecompose_one_pair(SVS, U)
  Q <- decomp$Q
  G <- decomp$G
  d <- decomp$eigenvalues

  # Q'G = I still holds
  expect_equal(crossprod(Q, G), diag(R), tolerance = 1e-12)

  # G diag(d) G' = U
  expect_equal(G %*% diag(d) %*% t(G), U, tolerance = 1e-10)

  # Only 1 non-zero eigenvalue
  expect_equal(sum(d > 1e-10), 1)
})

# =============================================================================
# Unit test: loglik_precomputed vs direct computation
# =============================================================================

test_that("loglik_precomputed matches direct computation (common_cov, K=1)", {
  set.seed(3)
  R <- 3
  J <- 20

  # SPD covariance
  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.5

  # Prior structure
  B <- matrix(rnorm(R * R), R, R)
  U <- crossprod(B)
  V_structure <- list(U)

  # Simulated betahats
  betahat <- matrix(rnorm(J * R), J, R)

  # Build cache
  svs_list <- list(SVS)
  cache <- mvsusieR:::precompute_eigen_cache(svs_list, V_structure,
                                              is_common_cov = TRUE)

  for (V_scalar in c(0, 0.5, 1, 2.5, 10)) {
    llik_fast <- mvsusieR:::loglik_precomputed(betahat, V_scalar, cache)
    expect_equal(ncol(llik_fast), 2)  # null + 1 component

    # Direct computation
    llik_direct <- matrix(0, J, 2)
    for (j in seq_len(J)) {
      # Null: N(0, SVS)
      llik_direct[j, 1] <- mvtnorm::dmvnorm(betahat[j, ], sigma = SVS,
                                              log = TRUE)
      # Non-null: N(0, SVS + V*U)
      llik_direct[j, 2] <- mvtnorm::dmvnorm(betahat[j, ],
                                              sigma = SVS + V_scalar * U,
                                              log = TRUE)
    }
    expect_equal(llik_fast, llik_direct, tolerance = 1e-10,
                 info = paste("V_scalar =", V_scalar))
  }
})

test_that("loglik_precomputed matches direct computation (common_cov, K=3)", {
  set.seed(4)
  R <- 2
  J <- 15

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.3

  # K=3 prior components
  V_structure <- lapply(1:3, function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)

  svs_list <- list(SVS)
  cache <- mvsusieR:::precompute_eigen_cache(svs_list, V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5
  llik_fast <- mvsusieR:::loglik_precomputed(betahat, V_scalar, cache)
  expect_equal(ncol(llik_fast), 4)  # null + 3 components

  llik_direct <- matrix(0, J, 4)
  for (j in seq_len(J)) {
    llik_direct[j, 1] <- mvtnorm::dmvnorm(betahat[j, ], sigma = SVS,
                                            log = TRUE)
    for (k in 1:3) {
      llik_direct[j, k + 1] <- mvtnorm::dmvnorm(betahat[j, ],
                                                  sigma = SVS + V_scalar * V_structure[[k]],
                                                  log = TRUE)
    }
  }
  expect_equal(llik_fast, llik_direct, tolerance = 1e-10)
})

test_that("loglik_precomputed matches direct (non-common_cov)", {
  set.seed(5)
  R <- 2
  J <- 10

  # Each variable has different SVS
  svs_list <- lapply(1:J, function(j) {
    A <- matrix(rnorm(R * R), R, R)
    crossprod(A) + diag(R) * (0.1 + j * 0.05)
  })

  B <- matrix(rnorm(R * R), R, R)
  U <- crossprod(B)
  V_structure <- list(U)

  betahat <- matrix(rnorm(J * R), J, R)

  cache <- mvsusieR:::precompute_eigen_cache(svs_list, V_structure,
                                              is_common_cov = FALSE)

  V_scalar <- 2.0
  llik_fast <- mvsusieR:::loglik_precomputed(betahat, V_scalar, cache)

  llik_direct <- matrix(0, J, 2)
  for (j in seq_len(J)) {
    llik_direct[j, 1] <- mvtnorm::dmvnorm(betahat[j, ], sigma = svs_list[[j]],
                                            log = TRUE)
    llik_direct[j, 2] <- mvtnorm::dmvnorm(betahat[j, ],
                                            sigma = svs_list[[j]] + V_scalar * U,
                                            log = TRUE)
  }
  expect_equal(llik_fast, llik_direct, tolerance = 1e-10)
})

# =============================================================================
# Unit test: posterior_precomputed vs direct computation
# =============================================================================

test_that("posterior_precomputed matches direct computation (K=1)", {
  set.seed(6)
  R <- 2
  J <- 10

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.3

  B <- matrix(rnorm(R * R), R, R)
  U <- crossprod(B)
  V_structure <- list(U)

  betahat <- matrix(rnorm(J * R), J, R)

  svs_list <- list(SVS)
  cache <- mvsusieR:::precompute_eigen_cache(svs_list, V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5

  # Posterior weights: simulate mostly component 1 (non-null)
  pi_V_post <- cbind(rep(0.1, J), rep(0.9, J))

  post <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache, pi_V_post)

  # Direct computation for each variable
  VU <- V_scalar * U
  Sigma_total <- SVS + VU
  Sigma_total_inv <- solve(Sigma_total)
  # Posterior cov = VU - VU Sigma_total^{-1} VU
  post_cov <- VU - VU %*% Sigma_total_inv %*% VU
  post_cov <- (post_cov + t(post_cov)) / 2

  for (j in seq_len(J)) {
    # Component 1 posterior mean
    m_k <- drop(VU %*% Sigma_total_inv %*% betahat[j, ])
    # Mixture-weighted mean
    expected_mean <- 0.1 * rep(0, R) + 0.9 * m_k
    expect_equal(post$post_mean[j, ], expected_mean, tolerance = 1e-10,
                 info = paste("j =", j))

    # Mixture-weighted second moment
    expected_m2 <- 0.1 * matrix(0, R, R) + 0.9 * (post_cov + tcrossprod(m_k))
    expect_equal(post$post_mean2[j, , ], expected_m2, tolerance = 1e-10,
                 info = paste("j =", j))
  }
})

# =============================================================================
# E2E: precompute=TRUE vs precompute=FALSE (individual data, matrix prior)
# =============================================================================

# For fixed V, both paths execute the same IBSS loop with identical V.
# The only source of error is different numerical decompositions
# (eigendecomposition vs Cholesky in C++), which gives ~1e-10 differences.
# Tolerance 1e-8 is appropriate.

test_that("E2E: precompute=TRUE matches FALSE for individual data, fixed V", {
  sim <- simulate_multivariate(r = 2)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(fit_yes$b1, fit_no$b1, tolerance = 1e-8)
  expect_equal(fit_yes$b2, fit_no$b2, tolerance = 1e-8)
  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 1e-8)
  expect_equal(fit_yes$lbf, fit_no$lbf, tolerance = 1e-8)
  expect_equal(fit_yes$KL, fit_no$KL, tolerance = 1e-8)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 1e-8)
  expect_equal(fit_yes$lfsr, fit_no$lfsr, tolerance = 1e-8)
})

# For estimated V (optim/EM/uniroot), the optimizer evaluates the loglik at
# many V values. The eigendecomposition and Cholesky C++ paths differ by ~1e-10
# per evaluation, but the optimizer can amplify these tiny differences, leading
# to slightly different convergence points. This is expected and inherent to
# using different numerical paths for the same math. Tolerance 0.05 is
# appropriate to verify the results are substantively equivalent.

test_that("E2E: precompute=TRUE matches FALSE for individual data, optim V", {
  sim <- simulate_multivariate(r = 2)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 0.05)
})

test_that("E2E: precompute=TRUE matches FALSE for individual data, EM V", {
  sim <- simulate_multivariate(r = 2)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 0.05)
})

test_that("E2E: precompute=TRUE matches FALSE with estimated residual variance", {
  sim <- simulate_multivariate(r = 2)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     estimate_residual_variance = TRUE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     estimate_residual_variance = TRUE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
  expect_equal(fit_yes$sigma2, fit_no$sigma2, tolerance = 0.05)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 0.05)
})

# =============================================================================
# E2E: precompute with mash (mixture) prior
# =============================================================================

test_that("E2E: precompute=TRUE matches FALSE for mash prior, fixed V", {
  sim <- simulate_multivariate(r = 2)
  mash_init <- create_mash_prior(Ulist = list(sim$V), grid = c(0.5, 1, 2))

  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = mash_init,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = mash_init,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(fit_yes$b1, fit_no$b1, tolerance = 1e-8)
  expect_equal(fit_yes$b2, fit_no$b2, tolerance = 1e-8)
  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 1e-8)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 1e-8)
  expect_equal(fit_yes$lfsr, fit_no$lfsr, tolerance = 1e-8)
})

test_that("E2E: precompute=TRUE matches FALSE for mash prior with EM", {
  sim <- simulate_multivariate(r = 2)
  mash_init <- create_mash_prior(Ulist = list(sim$V), grid = 1)

  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = mash_init,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = mash_init,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 0.05)
})

# =============================================================================
# E2E: sufficient statistics path
# =============================================================================

test_that("E2E: precompute=TRUE matches FALSE for sufficient statistics, fixed V", {
  sim <- simulate_multivariate(r = 2)
  X_colmeans <- colMeans(sim$X)
  Y_colmeans <- colMeans(sim$y)
  X_centered <- scale(sim$X, center = TRUE, scale = FALSE)
  Y_centered <- scale(sim$y, center = TRUE, scale = FALSE)
  XtX <- crossprod(X_centered)
  XtY <- crossprod(X_centered, Y_centered)
  YtY <- crossprod(Y_centered)
  n <- nrow(sim$X)

  fit_no  <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                X_colmeans = X_colmeans,
                                Y_colmeans = Y_colmeans,
                                prior_variance = sim$V,
                                estimate_prior_variance = FALSE,
                                estimate_residual_variance = FALSE,
                                precompute_covariances = FALSE,
                                max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                X_colmeans = X_colmeans,
                                Y_colmeans = Y_colmeans,
                                prior_variance = sim$V,
                                estimate_prior_variance = FALSE,
                                estimate_residual_variance = FALSE,
                                precompute_covariances = TRUE,
                                max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(fit_yes$b1, fit_no$b1, tolerance = 1e-8)
  expect_equal(fit_yes$b2, fit_no$b2, tolerance = 1e-8)
  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 1e-8)
})

test_that("E2E: precompute=TRUE matches FALSE for sufficient statistics, optim V", {
  sim <- simulate_multivariate(r = 2)
  X_colmeans <- colMeans(sim$X)
  Y_colmeans <- colMeans(sim$y)
  X_centered <- scale(sim$X, center = TRUE, scale = FALSE)
  Y_centered <- scale(sim$y, center = TRUE, scale = FALSE)
  XtX <- crossprod(X_centered)
  XtY <- crossprod(X_centered, Y_centered)
  YtY <- crossprod(Y_centered)
  n <- nrow(sim$X)

  fit_no  <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                X_colmeans = X_colmeans,
                                Y_colmeans = Y_colmeans,
                                prior_variance = sim$V,
                                estimate_prior_variance = TRUE,
                                estimate_prior_method = "optim",
                                estimate_residual_variance = FALSE,
                                precompute_covariances = FALSE,
                                max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                X_colmeans = X_colmeans,
                                Y_colmeans = Y_colmeans,
                                prior_variance = sim$V,
                                estimate_prior_variance = TRUE,
                                estimate_prior_method = "optim",
                                estimate_residual_variance = FALSE,
                                precompute_covariances = TRUE,
                                max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
})

# =============================================================================
# E2E: R=3 higher dimension
# =============================================================================

test_that("E2E: precompute=TRUE matches FALSE for R=3, fixed V", {
  sim <- simulate_multivariate(r = 3)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(fit_yes$b1, fit_no$b1, tolerance = 1e-8)
  expect_equal(fit_yes$b2, fit_no$b2, tolerance = 1e-8)
  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 1e-8)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 1e-8)
  expect_equal(fit_yes$lfsr, fit_no$lfsr, tolerance = 1e-8)
})

# =============================================================================
# E2E: R=1 univariate
# =============================================================================

test_that("E2E: precompute=TRUE matches FALSE for R=1, fixed V", {
  sim <- simulate_univariate()
  # For R=1, use a matrix prior to go through the MV path
  V_mat <- matrix(sim$V, 1, 1)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = V_mat,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = V_mat,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(fit_yes$b1, fit_no$b1, tolerance = 1e-8)
  expect_equal(fit_yes$b2, fit_no$b2, tolerance = 1e-8)
  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 1e-8)
})

# =============================================================================
# E2E: uniroot estimation method
# =============================================================================

test_that("E2E: precompute=TRUE matches FALSE with uniroot", {
  sim <- simulate_multivariate(r = 2)
  fit_no  <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "uniroot",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = FALSE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "uniroot",
                     estimate_residual_variance = FALSE,
                     precompute_covariances = TRUE,
                     max_iter = 50, tol = 1e-3, verbosity = 0)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 0.05)
})
