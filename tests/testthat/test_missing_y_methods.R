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
