context("Test uniroot prior variance estimation")

# =============================================================================
# Univariate (R=1): uniroot vs optim should give similar results
# =============================================================================

test_that("uniroot gives similar results to optim for R=1 individual data", {
  sim <- simulate_univariate()
  fit_optim <- mvsusie(sim$X, sim$y, L = 5,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       estimate_residual_variance = FALSE,
                       max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_uniroot <- mvsusie(sim$X, sim$y, L = 5,
                          estimate_prior_variance = TRUE,
                          estimate_prior_method = "uniroot",
                          estimate_residual_variance = FALSE,
                          max_iter = 50, tol = 1e-3, verbosity = 0)
  # Effects with signal should have similar V (null effects may differ
  # between optim and uniroot due to null-check threshold sensitivity)
  active <- which(fit_optim$V > 1)
  expect_equal(fit_uniroot$V[active], fit_optim$V[active], tolerance = 0.5)
  # PIPs should be similar
  expect_equal(fit_uniroot$pip, fit_optim$pip, tolerance = 0.05)
})

# =============================================================================
# Multivariate (R>1) with matrix prior: uniroot vs optim
# =============================================================================

test_that("uniroot gives similar results to optim for R>1 matrix prior", {
  sim <- simulate_multivariate(r = 2)
  fit_optim <- mvsusie(sim$X, sim$y, L = 5,
                       prior_variance = sim$V,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       estimate_residual_variance = FALSE,
                       max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_uniroot <- mvsusie(sim$X, sim$y, L = 5,
                          prior_variance = sim$V,
                          estimate_prior_variance = TRUE,
                          estimate_prior_method = "uniroot",
                          estimate_residual_variance = FALSE,
                          max_iter = 50, tol = 1e-3, verbosity = 0)
  # Effects with signal should have similar V (null effects may differ
  # between optim and uniroot due to null-check threshold sensitivity)
  active <- which(fit_optim$V > 1)
  expect_equal(fit_uniroot$V[active], fit_optim$V[active], tolerance = 0.5)
  # PIPs should be similar
  expect_equal(fit_uniroot$pip, fit_optim$pip, tolerance = 0.1)
})

# =============================================================================
# Multivariate with mixture (mash) prior: uniroot vs optim
# =============================================================================

test_that("uniroot gives similar results to optim for mixture prior", {
  sim <- simulate_multivariate(r = 2)
  mash_init <- create_mash_prior(Ulist = list(sim$V), grid = 1)
  fit_optim <- mvsusie(sim$X, sim$y, L = 5,
                       prior_variance = mash_init,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       estimate_residual_variance = FALSE,
                       max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_uniroot <- mvsusie(sim$X, sim$y, L = 5,
                          prior_variance = mash_init,
                          estimate_prior_variance = TRUE,
                          estimate_prior_method = "uniroot",
                          estimate_residual_variance = FALSE,
                          max_iter = 50, tol = 1e-3, verbosity = 0)
  # Effects with signal should have similar V
  active <- which(fit_optim$V > 1)
  expect_equal(fit_uniroot$V[active], fit_optim$V[active], tolerance = 0.5)
  # PIPs should be similar
  expect_equal(fit_uniroot$pip, fit_optim$pip, tolerance = 0.1)
})

# =============================================================================
# Sufficient statistics: uniroot vs optim
# =============================================================================

test_that("uniroot gives similar results to optim for sufficient statistics", {
  sim <- simulate_multivariate(r = 2)
  X_colmeans <- colMeans(sim$X)
  Y_colmeans <- colMeans(sim$y)
  X_centered <- scale(sim$X, center = TRUE, scale = FALSE)
  Y_centered <- scale(sim$y, center = TRUE, scale = FALSE)
  XtX <- crossprod(X_centered)
  XtY <- crossprod(X_centered, Y_centered)
  YtY <- crossprod(Y_centered)
  n <- nrow(sim$X)

  fit_optim <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                  X_colmeans = X_colmeans,
                                  Y_colmeans = Y_colmeans,
                                  prior_variance = sim$V,
                                  estimate_prior_variance = TRUE,
                                  estimate_prior_method = "optim",
                                  estimate_residual_variance = FALSE,
                                  max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_uniroot <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                                    X_colmeans = X_colmeans,
                                    Y_colmeans = Y_colmeans,
                                    prior_variance = sim$V,
                                    estimate_prior_variance = TRUE,
                                    estimate_prior_method = "uniroot",
                                    estimate_residual_variance = FALSE,
                                    max_iter = 50, tol = 1e-3, verbosity = 0)
  # V estimates should be close
  expect_equal(fit_uniroot$V, fit_optim$V, tolerance = 0.2)
  # Alpha should be similar
  expect_equal(fit_uniroot$alpha, fit_optim$alpha, tolerance = 0.1)
})

# =============================================================================
# RSS: uniroot vs optim
# =============================================================================

test_that("uniroot gives similar results to optim for RSS", with(simulate_multivariate(r = 2), {
  z <- sapply(1:ncol(y), function(j) {
    ss <- susieR:::univariate_regression(X, y[, j])
    ss$betahat / ss$sebetahat
  })
  R_ld <- cor(X)
  n <- nrow(X)

  fit_optim <- mvsusie_rss(z, R_ld, N = n, L = 5,
                           prior_variance = V,
                           estimate_prior_variance = TRUE,
                           estimate_prior_method = "optim",
                           estimate_residual_variance = FALSE,
                           max_iter = 50, tol = 1e-3, verbosity = 0)
  fit_uniroot <- mvsusie_rss(z, R_ld, N = n, L = 5,
                              prior_variance = V,
                              estimate_prior_variance = TRUE,
                              estimate_prior_method = "uniroot",
                              estimate_residual_variance = FALSE,
                              max_iter = 50, tol = 1e-3, verbosity = 0)
  # V estimates should be close
  expect_equal(fit_uniroot$V, fit_optim$V, tolerance = 0.2)
  # Alpha should be similar
  expect_equal(fit_uniroot$alpha, fit_optim$alpha, tolerance = 0.1)
}))

# =============================================================================
# Edge case: uniroot with no signal should set V to 0
# =============================================================================

test_that("uniroot correctly identifies null effects", {
  set.seed(42)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)  # pure noise, no signal
  fit <- mvsusie(X, y, L = 3,
                 estimate_prior_variance = TRUE,
                 estimate_prior_method = "uniroot",
                 estimate_residual_variance = FALSE,
                 max_iter = 50, tol = 1e-3, verbosity = 0)
  # All V should be zero or near-zero (no true effects)
  expect_true(all(fit$V < 0.5))
})

# =============================================================================
# Multivariate R=3 with mixture prior: uniroot convergence
# =============================================================================

test_that("uniroot converges for R=3 with mixture prior", {
  sim <- simulate_multivariate(r = 3)
  mash_init <- create_mash_prior(Ulist = list(sim$V), grid = c(0.5, 1, 2))
  fit <- mvsusie(sim$X, sim$y, L = 5,
                 prior_variance = mash_init,
                 estimate_prior_variance = TRUE,
                 estimate_prior_method = "uniroot",
                 estimate_residual_variance = FALSE,
                 max_iter = 100, tol = 1e-3, verbosity = 0)
  # Model should converge
  expect_true(fit$convergence$converged)
  # ELBO should be monotonically non-decreasing (or very close)
  elbo <- fit$elbo[!is.na(fit$elbo)]
  if (length(elbo) > 1) {
    diffs <- diff(elbo)
    # Allow tiny decreases due to numerical precision
    expect_true(all(diffs > -1e-3))
  }
})
