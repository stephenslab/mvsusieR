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
    expect_equal(post$mu[j, ], expected_mean, tolerance = 1e-10,
                 info = paste("j =", j))

    # Mixture-weighted second moment
    expected_m2 <- 0.1 * matrix(0, R, R) + 0.9 * (post_cov + tcrossprod(m_k))
    expect_equal(post$mu2[j, , ], expected_m2, tolerance = 1e-10,
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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(get_b1(fit_yes), get_b1(fit_no), tolerance = 1e-8)
  expect_equal(get_b2(fit_yes), get_b2(fit_no), tolerance = 1e-8)
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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     estimate_residual_variance = TRUE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = mash_init,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(get_b1(fit_yes), get_b1(fit_no), tolerance = 1e-8)
  expect_equal(get_b2(fit_yes), get_b2(fit_no), tolerance = 1e-8)
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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = mash_init,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "EM",
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

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
                                precompute_cache = FALSE,
                                max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                X_colmeans = X_colmeans,
                                Y_colmeans = Y_colmeans,
                                prior_variance = sim$V,
                                estimate_prior_variance = FALSE,
                                estimate_residual_variance = FALSE,
                                precompute_cache = TRUE,
                                max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(get_b1(fit_yes), get_b1(fit_no), tolerance = 1e-8)
  expect_equal(get_b2(fit_yes), get_b2(fit_no), tolerance = 1e-8)
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
                                precompute_cache = FALSE,
                                max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                X_colmeans = X_colmeans,
                                Y_colmeans = Y_colmeans,
                                prior_variance = sim$V,
                                estimate_prior_variance = TRUE,
                                estimate_prior_method = "optim",
                                estimate_residual_variance = FALSE,
                                precompute_cache = TRUE,
                                max_iter = 5, tol = 1e-3, verbose = FALSE)

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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(get_b1(fit_yes), get_b1(fit_no), tolerance = 1e-8)
  expect_equal(get_b2(fit_yes), get_b2(fit_no), tolerance = 1e-8)
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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = V_mat,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_equal(fit_yes$alpha, fit_no$alpha, tolerance = 1e-8)
  expect_equal(get_b1(fit_yes), get_b1(fit_no), tolerance = 1e-8)
  expect_equal(get_b2(fit_yes), get_b2(fit_no), tolerance = 1e-8)
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
                     precompute_cache = FALSE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)
  fit_yes <- mvsusie(sim$X, sim$y, L = 5,
                     prior_variance = sim$V,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "uniroot",
                     estimate_residual_variance = FALSE,
                     precompute_cache = TRUE,
                     max_iter = 5, tol = 1e-3, verbose = FALSE)

  expect_equal(fit_yes$pip, fit_no$pip, tolerance = 0.05)
  expect_equal(fit_yes$V, fit_no$V, tolerance = 0.05)
  expect_equal(fit_yes$fitted, fit_no$fitted, tolerance = 0.05)
})

# =============================================================================
# Optimization #2: Cholesky caching in precompute_eigen_cache
# =============================================================================

test_that("#2: precompute_eigen_cache common-cov matches per-K eigendecompose", {
  set.seed(20)
  R <- 5
  K <- 4

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.2

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)
  for (k in seq_len(K)) {
    ref <- mvsusieR:::eigendecompose_one_pair(SVS, V_structure[[k]])

    expect_equal(cache$components[[k]]$eigenvalues, ref$eigenvalues,
                 tolerance = 1e-12, info = paste("eigenvalues, k =", k))

    expect_equal(crossprod(cache$components[[k]]$Q, cache$components[[k]]$G),
                 diag(R), tolerance = 1e-12, info = paste("Q'G=I, k =", k))

    d <- cache$components[[k]]$eigenvalues
    G <- cache$components[[k]]$G
    expect_equal(G %*% diag(d) %*% t(G), V_structure[[k]],
                 tolerance = 1e-12, info = paste("G*d*G'=U, k =", k))

    expect_equal(cache$components[[k]]$log_det_svs, ref$log_det_svs,
                 tolerance = 1e-12, info = paste("log_det_svs, k =", k))
  }
})

test_that("#2: precompute_eigen_cache_R matches cached version", {
  set.seed(21)
  R <- 4
  K <- 3

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.3

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  cache_main <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                                    is_common_cov = TRUE)
  cache_R <- mvsusieR:::precompute_eigen_cache_R(list(SVS), V_structure,
                                                   is_common_cov = TRUE)

  expect_equal(cache_main$log_det_svs, cache_R$log_det_svs, tolerance = 1e-12)
  for (k in seq_len(K)) {
    expect_equal(cache_main$components[[k]]$eigenvalues,
                 cache_R$components[[k]]$eigenvalues, tolerance = 1e-12)
  }
})

# =============================================================================
# Optimization #1: BQ_cache in loglik_precomputed
# =============================================================================

test_that("#1: loglik_precomputed with BQ_cache matches without cache", {
  set.seed(10)
  R <- 4
  J <- 25
  K <- 5

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.3

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  BQ_cache <- lapply(seq_len(K), function(k)
    betahat %*% cache$components[[k]]$Q)

  for (V_scalar in c(0.1, 0.5, 1, 3, 10)) {
    llik_no_cache <- mvsusieR:::loglik_precomputed(betahat, V_scalar, cache,
                                                    BQ_cache = NULL)
    llik_cached   <- mvsusieR:::loglik_precomputed(betahat, V_scalar, cache,
                                                    BQ_cache = BQ_cache)
    expect_equal(llik_cached, llik_no_cache, tolerance = 1e-14,
                 info = paste("V_scalar =", V_scalar))
  }
})

# =============================================================================
# Optimization #5: BQ_cache in posterior_precomputed
# =============================================================================

test_that("#5: posterior with BQ_cache matches without (no reduce)", {
  set.seed(50)
  R <- 3
  J <- 15
  K <- 3

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.4

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5
  pi_V_post <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_V_post / rowSums(pi_V_post)

  BQ_cache <- lapply(seq_len(K), function(k)
    betahat %*% cache$components[[k]]$Q)

  post_no  <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                pi_V_post, BQ_cache = NULL)
  post_yes <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                pi_V_post, BQ_cache = BQ_cache)

  expect_equal(post_yes$mu, post_no$mu, tolerance = 1e-10)
  expect_equal(post_yes$mu2, post_no$mu2, tolerance = 1e-10)
  expect_equal(post_yes$post_neg, post_no$post_neg, tolerance = 1e-10)
  expect_equal(post_yes$post_zero, post_no$post_zero, tolerance = 1e-10)
  expect_equal(post_yes$prior_scale_em_update,
               post_no$prior_scale_em_update, tolerance = 1e-10)
})

test_that("#5: posterior with BQ_cache matches without (reduce)", {
  set.seed(51)
  R <- 3
  J <- 15
  K <- 3

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.4

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5
  pi_V_post <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_V_post / rowSums(pi_V_post)

  alpha <- runif(J); alpha <- alpha / sum(alpha)
  d_var <- runif(J, 0.5, 2)
  v_inv <- solve(SVS)
  reduce_params <- list(alpha = alpha, d = d_var, v_inv = v_inv)

  BQ_cache <- lapply(seq_len(K), function(k)
    betahat %*% cache$components[[k]]$Q)

  post_no  <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                pi_V_post,
                                                reduce_params = reduce_params,
                                                BQ_cache = NULL)
  post_yes <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                pi_V_post,
                                                reduce_params = reduce_params,
                                                BQ_cache = BQ_cache)

  expect_equal(post_yes$mu, post_no$mu, tolerance = 1e-10)
  expect_equal(post_yes$bxxb, post_no$bxxb, tolerance = 1e-10)
  expect_equal(post_yes$vbxxb, post_no$vbxxb, tolerance = 1e-10)
  expect_equal(post_yes$alpha_mu2_sum, post_no$alpha_mu2_sum, tolerance = 1e-10)
  expect_equal(post_yes$mu2_diag, post_no$mu2_diag, tolerance = 1e-10)
})

# =============================================================================
# Optimization #3: loglik_common_rcpp matches R reference
# =============================================================================

test_that("#3: loglik_common_rcpp matches loglik_precomputed_R", {
  set.seed(30)
  R <- 4
  J <- 20
  K <- 5

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.3

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  for (V_scalar in c(0, 0.1, 1, 5)) {
    llik_cpp <- mvsusieR:::loglik_precomputed(betahat, V_scalar, cache)
    llik_R   <- mvsusieR:::loglik_precomputed_R(betahat, V_scalar, cache)
    expect_equal(llik_cpp, llik_R, tolerance = 1e-12,
                 info = paste("V_scalar =", V_scalar))
  }
})

# =============================================================================
# Optimization #4: posterior_common_rcpp matches R reference
# =============================================================================

test_that("#4: posterior_common_rcpp matches R reference (no reduce, no EM)", {
  set.seed(40)
  R <- 3
  J <- 15
  K <- 3

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.4

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5
  pi_V_post <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_V_post / rowSums(pi_V_post)

  post_cpp <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                pi_V_post)
  post_R   <- mvsusieR:::posterior_precomputed_R(betahat, V_scalar, cache,
                                                  pi_V_post)

  expect_equal(post_cpp$mu, post_R$mu, tolerance = 1e-10)
  expect_equal(post_cpp$mu2, post_R$mu2, tolerance = 1e-10)
  expect_equal(post_cpp$post_neg, post_R$post_neg, tolerance = 1e-10)
  expect_equal(post_cpp$post_zero, post_R$post_zero, tolerance = 1e-10)
  expect_equal(as.vector(post_cpp$prior_scale_em_update),
               as.vector(post_R$prior_scale_em_update), tolerance = 1e-10)
})

test_that("#4: posterior_common_rcpp matches R reference (with EM weights)", {
  set.seed(41)
  R <- 3
  J <- 15
  K <- 3

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.4

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5
  pi_V_post <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_V_post / rowSums(pi_V_post)

  em_var_wt <- matrix(runif((K + 1) * J), K + 1, J)
  em_var_wt <- em_var_wt / rowSums(em_var_wt)

  post_cpp <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                pi_V_post,
                                                em_var_wt = em_var_wt)
  post_R   <- mvsusieR:::posterior_precomputed_R(betahat, V_scalar, cache,
                                                  pi_V_post,
                                                  em_var_wt = em_var_wt)

  expect_equal(post_cpp$mu, post_R$mu, tolerance = 1e-10)
  expect_equal(post_cpp$mu2, post_R$mu2, tolerance = 1e-10)
  expect_equal(as.vector(post_cpp$prior_scale_em_update),
               as.vector(post_R$prior_scale_em_update), tolerance = 1e-10)
})

test_that("#4: posterior_common_rcpp reduce matches R full + reduce", {
  set.seed(42)
  R <- 3
  J <- 15
  K <- 3

  A <- matrix(rnorm(R * R), R, R)
  SVS <- crossprod(A) + diag(R) * 0.4

  V_structure <- lapply(seq_len(K), function(k) {
    B <- matrix(rnorm(R * R), R, R)
    crossprod(B)
  })

  betahat <- matrix(rnorm(J * R), J, R)
  cache <- mvsusieR:::precompute_eigen_cache(list(SVS), V_structure,
                                              is_common_cov = TRUE)

  V_scalar <- 1.5
  pi_V_post <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_V_post / rowSums(pi_V_post)

  alpha <- runif(J); alpha <- alpha / sum(alpha)
  d_var <- runif(J, 0.5, 2)
  v_inv <- solve(SVS)
  reduce_params <- list(alpha = alpha, d = d_var, v_inv = v_inv)

  post_reduce <- mvsusieR:::posterior_precomputed(betahat, V_scalar, cache,
                                                   pi_V_post,
                                                   reduce_params = reduce_params)

  post_R <- mvsusieR:::posterior_precomputed_R(betahat, V_scalar, cache,
                                                pi_V_post)

  bxxb <- matrix(0, R, R)
  alpha_mu2_sum <- matrix(0, R, R)
  mu2_diag <- matrix(0, J, R)
  for (j in seq_len(J)) {
    mu2_j <- post_R$mu2[j, , ]
    if (!is.matrix(mu2_j)) dim(mu2_j) <- c(R, R)
    bxxb <- bxxb + d_var[j] * alpha[j] * mu2_j
    alpha_mu2_sum <- alpha_mu2_sum + alpha[j] * mu2_j
    mu2_diag[j, ] <- diag(mu2_j)
  }
  vbxxb <- sum(v_inv * bxxb)

  expect_equal(post_reduce$mu, post_R$mu, tolerance = 1e-10)
  expect_equal(post_reduce$bxxb, bxxb, tolerance = 1e-10)
  expect_equal(post_reduce$vbxxb, vbxxb, tolerance = 1e-10)
  expect_equal(post_reduce$alpha_mu2_sum, alpha_mu2_sum, tolerance = 1e-10)
  expect_equal(post_reduce$mu2_diag, mu2_diag, tolerance = 1e-10)
})
