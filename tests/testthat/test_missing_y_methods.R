context("Test approximate and exact missing data methods")

# =============================================================================
# Approximate method: basic functionality
# =============================================================================

test_that("approximate method runs and converges with 20% MCAR", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 2, y_missing = 0.2)
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = sim$V,
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 30, tol = 1e-3, verbosity = 0)
  # Model should produce valid output
  expect_true(!is.null(fit$alpha))
  expect_true(!is.null(fit$pip))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  # At least some PIPs should be non-trivial (signal present)
  expect_true(max(fit$pip) > 0.1)
})

# =============================================================================
# Exact method: basic functionality
# =============================================================================

test_that("exact method runs and converges with 20% MCAR", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 2, y_missing = 0.2)
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = sim$V,
                 missing_y_method = "exact",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 30, tol = 1e-3, verbosity = 0)
  expect_true(!is.null(fit$alpha))
  expect_true(!is.null(fit$pip))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_true(max(fit$pip) > 0.1)
})

# =============================================================================
# No-missing identity: all methods on complete data give same results
# =============================================================================

test_that("approximate method equals standard when no data is missing", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 2)
  V <- cov(sim$y)
  fit_standard <- mvsusie(sim$X, sim$y, L = 5,
                          prior_variance = sim$V,
                          missing_y_method = "impute",
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          residual_variance = V,
                          max_iter = 10, verbosity = 0)
  fit_approx <- mvsusie(sim$X, sim$y, L = 5,
                         prior_variance = sim$V,
                         missing_y_method = "approximate",
                         estimate_prior_variance = FALSE,
                         estimate_residual_variance = FALSE,
                         residual_variance = V,
                         max_iter = 10, verbosity = 0)
  # When no data is missing, approximate should fall back to impute
  expect_equal(fit_standard$alpha, fit_approx$alpha)
  expect_equal(fit_standard$pip, fit_approx$pip)
})

# =============================================================================
# Approximate vs impute: PIPs should be correlated
# =============================================================================

test_that("approximate and impute PIPs are correlated with 20% MCAR", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 2, y_missing = 0.2)
  V <- cov(sim$y)
  fit_impute <- mvsusie(sim$X, sim$y_missing, L = 5,
                         prior_variance = sim$V,
                         missing_y_method = "impute",
                         estimate_prior_variance = TRUE,
                         estimate_residual_variance = FALSE,
                         residual_variance = V,
                         max_iter = 30, tol = 1e-3, verbosity = 0)
  fit_approx <- mvsusie(sim$X, sim$y_missing, L = 5,
                          prior_variance = sim$V,
                          missing_y_method = "approximate",
                          estimate_prior_variance = TRUE,
                          estimate_residual_variance = FALSE,
                          residual_variance = V,
                          max_iter = 30, tol = 1e-3, verbosity = 0)
  # PIPs from both methods should be positively correlated
  pip_cor <- cor(fit_impute$pip, fit_approx$pip)
  expect_true(pip_cor > 0.5)
})

# =============================================================================
# Exact ELBO monotonicity
# =============================================================================

test_that("exact method has monotonically non-decreasing ELBO", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 2, y_missing = 0.2)
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = sim$V,
                 missing_y_method = "exact",
                 estimate_prior_variance = FALSE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 20, tol = 1e-3, verbosity = 0)
  elbo <- fit$elbo[!is.na(fit$elbo)]
  if (length(elbo) > 1) {
    diffs <- diff(elbo)
    # Allow tiny decreases from numerical precision
    expect_true(all(diffs > -1e-3))
  }
})

# =============================================================================
# Null effects: pure noise should give near-zero prior variance
# =============================================================================

test_that("approximate correctly identifies null effects in pure noise", {
  set.seed(42)
  n <- 100; p <- 50; r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  # Introduce 20% MCAR
  Y_miss <- Y
  for (i in 1:n) Y_miss[i, runif(r) <= 0.2] <- NA
  fit <- mvsusie(X, Y_miss, L = 3,
                 prior_variance = diag(r),
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = diag(r),
                 max_iter = 20, tol = 1e-3, verbosity = 0)
  # All V should be zero or near-zero
  expect_true(all(fit$V < 1.0))
})

# =============================================================================
# R=1 fallback: missing_y_method is ignored for univariate
# =============================================================================

test_that("missing_y_method is ignored for R=1", {
  sim <- simulate_univariate()
  fit1 <- mvsusie(sim$X, sim$y, L = 5,
                  missing_y_method = "approximate",
                  estimate_prior_variance = FALSE,
                  max_iter = 5, verbosity = 0)
  fit2 <- mvsusie(sim$X, sim$y, L = 5,
                  missing_y_method = "impute",
                  estimate_prior_variance = FALSE,
                  max_iter = 5, verbosity = 0)
  expect_equal(fit1$alpha, fit2$alpha)
})

# =============================================================================
# estimate_residual_variance = TRUE silently overridden for approximate/exact
# =============================================================================

test_that("approximate/exact silently force estimate_residual_variance = FALSE", {
  sim <- simulate_multivariate(n = 50, p = 30, r = 2, y_missing = 0.2)
  # Should not error; estimate_residual_variance is silently set to FALSE
  expect_message(
    fit <- mvsusie(sim$X, sim$y_missing, L = 3,
                   prior_variance = sim$V,
                   missing_y_method = "approximate",
                   estimate_residual_variance = TRUE,
                   max_iter = 5, verbosity = 1),
    "estimate_residual_variance"
  )
  expect_true(!is.null(fit$alpha))
})

# =============================================================================
# Precomputation correctness: optimized svs_inv matches naive per-sample loop
# =============================================================================

test_that("precomputed svs_inv matches naive per-sample computation", {
  set.seed(123)
  n <- 30; p <- 10; r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  Y_missing <- matrix(FALSE, n, r)
  Y_missing[1:5, 1] <- TRUE
  Y_missing[6:10, 2] <- TRUE
  Y[Y_missing] <- NA
  V <- diag(r)

  # Create data with approximate method
  data <- mvsusieR:::create_mvsusie_data(X, Y, center = TRUE, scale = TRUE,
                                          missing_y_method = "approximate")
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)
  my <- data$miss3d

  # Naive per-sample svs_inv computation
  J <- data$p
  for (j in 1:J) {
    naive_svs_inv_j <- matrix(0, r, r)
    for (i in 1:n) {
      xij <- my$X_3d[i, j, ]
      Vinv_i <- my$Vinv[[my$pattern_assign[i]]]
      naive_svs_inv_j <- naive_svs_inv_j +
        t(Vinv_i * xij) * xij
    }
    expect_equal(data$svs_inv[[j]], naive_svs_inv_j, tolerance = 1e-10,
                 info = paste("svs_inv mismatch at j =", j))
  }
})

# =============================================================================
# XtR correctness: optimized matches naive
# =============================================================================

test_that("compute_XtR_3d matches naive per-variable computation", {
  set.seed(123)
  n <- 30; p <- 10; r <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * r), n, r)
  Y_missing <- matrix(FALSE, n, r)
  Y_missing[1:5, 1] <- TRUE
  Y_missing[6:10, 2] <- TRUE
  Y[Y_missing] <- NA
  V <- diag(r)

  data <- mvsusieR:::create_mvsusie_data(X, Y, center = TRUE, scale = TRUE,
                                          missing_y_method = "approximate")
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)
  my <- data$miss3d

  # Create a test residual
  R_mat <- matrix(rnorm(n * r), n, r)
  R_mat[!my$Y_non_missing] <- 0

  # Optimized
  XtR_opt <- mvsusieR:::compute_XtR_3d(data, R_mat)

  # Naive: per-variable computation
  Vinv <- my$Vinv
  VinvR <- matrix(0, n, r)
  for (i in 1:n) {
    VinvR[i, ] <- as.numeric(Vinv[[my$pattern_assign[i]]] %*% R_mat[i, ])
  }
  XtR_naive <- t(sapply(1:p, function(j) {
    colSums(my$X_3d[, j, ] * VinvR)
  }))

  expect_equal(XtR_opt, XtR_naive, tolerance = 1e-10)
})

# =============================================================================
# Cross-population: non-overlapping sample groups
# =============================================================================

test_that("approximate method handles non-overlapping sample groups", {
  set.seed(42)
  n <- 100; p <- 50; r <- 4
  X <- matrix(rnorm(n * p), n, p)
  beta <- matrix(0, p, r)
  beta[1:3, ] <- matrix(rnorm(12), 3, 4)
  Y <- X %*% beta + matrix(rnorm(n * r), n, r)

  # Samples 1-50 observed on conditions 1-2 only
  # Samples 51-100 observed on conditions 3-4 only
  Y_miss <- Y
  Y_miss[1:50, 3:4] <- NA
  Y_miss[51:100, 1:2] <- NA

  V <- diag(r)
  fit <- mvsusie(X, Y_miss, L = 5,
                 prior_variance = V,
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = V,
                 max_iter = 30, tol = 1e-3, verbosity = 0)
  # Model should converge and produce valid output
  expect_true(!is.null(fit$pip))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

# =============================================================================
# Diagonal V: approximate == exact (numerical identity)
# =============================================================================

test_that("approximate equals exact when residual variance is diagonal", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 2, y_missing = 0.2)
  V_diag <- diag(diag(cov(sim$y)))  # diagonal residual variance
  fit_approx <- mvsusie(sim$X, sim$y_missing, L = 5,
                         prior_variance = sim$V,
                         missing_y_method = "approximate",
                         estimate_prior_variance = FALSE,
                         estimate_residual_variance = FALSE,
                         residual_variance = V_diag,
                         max_iter = 10, verbosity = 0)
  fit_exact <- mvsusie(sim$X, sim$y_missing, L = 5,
                        prior_variance = sim$V,
                        missing_y_method = "exact",
                        estimate_prior_variance = FALSE,
                        estimate_residual_variance = FALSE,
                        residual_variance = V_diag,
                        max_iter = 10, verbosity = 0)
  # With diagonal V, approximate and exact should give very similar results
  expect_equal(fit_approx$alpha, fit_exact$alpha, tolerance = 0.1)
  expect_equal(fit_approx$pip, fit_exact$pip, tolerance = 0.1)
})

# =============================================================================
# R=3 with mixture prior: convergence check
# =============================================================================

test_that("approximate converges for R=3 with mixture prior", {
  sim <- simulate_multivariate(n = 100, p = 50, r = 3, y_missing = 0.15)
  mash_init <- create_mash_prior(Ulist = list(sim$V), grid = c(0.5, 1, 2))
  fit <- mvsusie(sim$X, sim$y_missing, L = 5,
                 prior_variance = mash_init,
                 missing_y_method = "approximate",
                 estimate_prior_variance = TRUE,
                 estimate_residual_variance = FALSE,
                 residual_variance = cov(sim$y),
                 max_iter = 30, tol = 1e-3, verbosity = 0)
  expect_true(!is.null(fit$pip))
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
})

# =============================================================================
# C++ vs R implementation comparison tests
# =============================================================================

test_that("C++ compute_VinvR_3d matches R implementation", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "approximate")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  mat <- matrix(rnorm(data$n * data$R), data$n, data$R)
  res_cpp <- mvsusieR:::compute_VinvR_3d(data, mat)
  res_r   <- mvsusieR:::compute_VinvR_3d_R(data, mat)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("C++ compute_XtR_3d matches R for approximate method", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "approximate")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  R_mat <- matrix(rnorm(data$n * data$R), data$n, data$R)
  res_cpp <- mvsusieR:::compute_XtR_3d(data, R_mat)
  res_r   <- mvsusieR:::compute_XtR_3d_R(data, R_mat)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("C++ compute_XtR_3d matches R for exact method", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "exact")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  R_mat <- matrix(rnorm(data$n * data$R), data$n, data$R)
  res_cpp <- mvsusieR:::compute_XtR_3d(data, R_mat)
  res_r   <- mvsusieR:::compute_XtR_3d_R(data, R_mat)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("C++ compute_Xb_3d matches R for approximate method", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "approximate")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  b <- matrix(rnorm(data$p * data$R), data$p, data$R)
  res_cpp <- mvsusieR:::compute_Xb_3d(data, b)
  res_r   <- mvsusieR:::compute_Xb_3d_R(data, b)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("C++ compute_Xb_3d matches R for exact method", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "exact")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  b <- matrix(rnorm(data$p * data$R), data$p, data$R)
  res_cpp <- mvsusieR:::compute_Xb_3d(data, b)
  res_r   <- mvsusieR:::compute_Xb_3d_R(data, b)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

test_that("C++ compute_betahat_3d matches R implementation", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "approximate")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  XtR <- matrix(rnorm(data$p * data$R), data$p, data$R)
  res_cpp <- mvsusieR:::compute_betahat_3d(data, XtR)
  res_r   <- mvsusieR:::compute_betahat_3d_R(data, XtR)
  expect_equal(res_cpp, res_r, tolerance = 1e-10)
})

# =============================================================================
# C++ vs R comparison tests: Xbar from precomputed sums
# =============================================================================

test_that("C++ compute_Xbar_from_sums matches naive N-loop", {
  sim <- simulate_multivariate(n = 80, p = 30, r = 3, y_missing = 0.2)
  data <- mvsusieR:::create_mvsusie_data(sim$X, sim$y_missing,
            missing_y_method = "exact")
  V <- cov(sim$y)
  data <- mvsusieR:::set_mvsusie_residual_variance(data, V)

  # The Xbar was computed during standardize_3d using C++.
  # Verify against the naive N-loop computation.
  my <- data$miss3d
  J <- data$p
  N <- data$n
  R_dim <- data$R

  # Naive N-loop Xbar
  Xbar_naive <- array(0, dim = c(J, R_dim, R_dim))
  Vinv <- my$Vinv
  X_3d <- my$X_3d
  pattern_assign <- my$pattern_assign
  for (j in seq_len(J)) {
    inner <- matrix(0, R_dim, R_dim)
    for (i in seq_len(N)) {
      inner <- inner + t(t(Vinv[[pattern_assign[i]]]) * X_3d[i, j, ])
    }
    Xbar_naive[j, , ] <- my$Vinvsuminv %*% inner
  }

  expect_equal(my$Xbar, Xbar_naive, tolerance = 1e-10)
})

# =============================================================================
# C++ vs R comparison tests: eigendecomposition precomputation
# =============================================================================

test_that("C++ precompute_eigen_cache matches R for non-common-cov", {
  set.seed(42)
  J <- 20
  R <- 3
  K <- 4

  # Generate J different SVS matrices (SPD)
  svs <- lapply(seq_len(J), function(j) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) + diag(R) * 0.5
  })

  # Generate K prior structure matrices (PSD)
  V_structure <- lapply(seq_len(K), function(k) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) / R
  })

  cache_cpp <- mvsusieR:::precompute_eigen_cache(svs, V_structure, FALSE)
  cache_r   <- mvsusieR:::precompute_eigen_cache_R(svs, V_structure, FALSE)

  expect_equal(as.numeric(cache_cpp$log_det_svs),
               as.numeric(cache_r$log_det_svs), tolerance = 1e-10)

  # Eigenvalues should match (up to ordering differences)
  for (k in seq_len(K)) {
    for (j in seq_len(J)) {
      eig_cpp <- sort(cache_cpp$components[[k]]$eigenvalues[, j], decreasing = TRUE)
      eig_r   <- sort(cache_r$components[[k]]$eigenvalues[, j], decreasing = TRUE)
      expect_equal(eig_cpp, eig_r, tolerance = 1e-10)
    }
  }
})

test_that("C++ loglik_non_common matches R implementation", {
  set.seed(42)
  J <- 20
  R <- 3
  K <- 4

  svs <- lapply(seq_len(J), function(j) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) + diag(R) * 0.5
  })
  V_structure <- lapply(seq_len(K), function(k) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) / R
  })

  betahat <- matrix(rnorm(J * R), J, R)
  V_scalar <- 1.5

  # Use C++ cache (same format) but compare loglik outputs
  cache <- mvsusieR:::precompute_eigen_cache(svs, V_structure, FALSE)

  llik_cpp <- mvsusieR:::loglik_non_common_cpp(
    betahat, V_scalar, cache$log_det_svs, cache$components)
  llik_r   <- mvsusieR:::loglik_precomputed_R(betahat, V_scalar, cache)

  expect_equal(llik_cpp, llik_r, tolerance = 1e-10)
})

test_that("C++ posterior_non_common matches R implementation", {
  set.seed(42)
  J <- 15
  R <- 3
  K <- 3

  svs <- lapply(seq_len(J), function(j) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) + diag(R) * 0.5
  })
  V_structure <- lapply(seq_len(K), function(k) {
    M <- matrix(rnorm(R * R), R, R)
    crossprod(M) / R
  })

  betahat <- matrix(rnorm(J * R), J, R)
  V_scalar <- 1.5

  cache <- mvsusieR:::precompute_eigen_cache(svs, V_structure, FALSE)

  # Generate random mixture weights (J x (K+1))
  pi_raw <- matrix(runif(J * (K + 1)), J, K + 1)
  pi_V_post <- pi_raw / rowSums(pi_raw)

  post_cpp <- mvsusieR:::posterior_non_common_cpp(
    betahat, V_scalar, cache$components, pi_V_post, matrix(0, 0, 0))
  post_r   <- mvsusieR:::posterior_precomputed_R(
    betahat, V_scalar, cache, pi_V_post)

  expect_equal(post_cpp$post_mean, post_r$post_mean, tolerance = 1e-8)
  expect_equal(post_cpp$post_mean2, post_r$post_mean2, tolerance = 1e-8)
  expect_equal(post_cpp$post_neg, post_r$post_neg, tolerance = 1e-6)
  expect_equal(post_cpp$post_zero, post_r$post_zero, tolerance = 1e-10)
  expect_equal(as.numeric(post_cpp$prior_scale_em_update),
               as.numeric(post_r$prior_scale_em_update), tolerance = 1e-8)
})

test_that("C++ accumulate_post_mean2_common matches R loop", {
  set.seed(42)
  J <- 30
  R <- 4

  post_mean2 <- array(rnorm(J * R * R), c(J, R, R))
  M_k <- matrix(rnorm(J * R), J, R)
  C_k_raw <- matrix(rnorm(R * R), R, R)
  C_k <- crossprod(C_k_raw) / R  # SPD
  w_k <- runif(J)

  # R loop
  pm2_r <- post_mean2
  for (j in seq_len(J)) {
    pm2_r[j, , ] <- pm2_r[j, , ] + w_k[j] * (C_k + tcrossprod(M_k[j, ]))
  }

  # C++ version
  pm2_cpp <- mvsusieR:::accumulate_post_mean2_common_cpp(post_mean2, M_k, C_k, w_k)

  expect_equal(pm2_cpp, pm2_r, tolerance = 1e-10)
})
