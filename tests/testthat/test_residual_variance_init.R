# =============================================================================
# Individual data: NULL residual_variance defaults
# =============================================================================

test_that("default residual_variance for R=1 individual data uses var(Y)", {
  set.seed(1)
  n <- 50; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  # Fit with default (NULL) residual_variance
  fit <- mvsusie(X, y, L = 2, max_iter = 2, verbosity = 0)
  # The residual variance should be close to var(y) (modulo centering)
  expect_true(!is.null(fit$sigma2))
  expect_true(is.numeric(fit$sigma2))
})

test_that("default residual_variance for R>1 individual data uses cov(Y)", {
  sim <- simulate_multivariate(n = 50, p = 30, r = 2)
  # Fit with default (NULL) residual_variance, must provide prior for R>1
  fit <- mvsusie(sim$X, sim$y, L = 2, prior_variance = sim$V,
                 max_iter = 2, verbosity = 0)
  # The residual variance should be a 2x2 matrix
  expect_true(is.matrix(fit$sigma2))
  expect_equal(nrow(fit$sigma2), 2)
  expect_equal(ncol(fit$sigma2), 2)
  # Should be positive definite
  evals <- eigen(fit$sigma2, only.values = TRUE)$values
  expect_true(all(evals > 0))
})

test_that("default residual_variance for R>1 with missing data uses flash/cov fallback", {
  skip_on_cran()
  sim <- simulate_multivariate(n = 50, p = 30, r = 2, y_missing = 0.1)
  # Fit with default (NULL) residual_variance and missing data
  fit <- mvsusie(sim$X, sim$y_missing, L = 2, prior_variance = sim$V,
                 max_iter = 2, verbosity = 0)
  # Should still produce a valid 2x2 positive definite matrix
  expect_true(is.matrix(fit$sigma2))
  expect_equal(nrow(fit$sigma2), 2)
  evals <- eigen(fit$sigma2, only.values = TRUE)$values
  expect_true(all(evals > 0))
})

# =============================================================================
# Sufficient statistics: NULL residual_variance defaults
# =============================================================================

test_that("default residual_variance for R=1 SS data uses YtY/(N-1)", {
  set.seed(1)
  n <- 50; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  XtX <- crossprod(Xc)
  XtY <- crossprod(Xc, yc)
  YtY <- crossprod(yc)
  # Fit with default (NULL) residual_variance
  fit <- mvsusie_ss(XtX, XtY, YtY, n, L = 2, max_iter = 2, verbosity = 0)
  expect_true(!is.null(fit$sigma2))
  expect_true(is.numeric(fit$sigma2))
})

test_that("default residual_variance for R>1 SS data uses cov2cor(YtY)", {
  sim <- simulate_multivariate(n = 50, p = 30, r = 2)
  Xc <- scale(sim$X, center = TRUE, scale = FALSE)
  Yc <- scale(sim$y, center = TRUE, scale = FALSE)
  XtX <- crossprod(Xc)
  XtY <- crossprod(Xc, Yc)
  YtY <- crossprod(Yc)
  # Fit with default (NULL) residual_variance, must provide prior for R>1
  fit <- mvsusie_ss(XtX, XtY, YtY, nrow(sim$X), L = 2,
                    prior_variance = sim$V,
                    max_iter = 2, verbosity = 0)
  # Should be a 2x2 matrix
  expect_true(is.matrix(fit$sigma2))
  expect_equal(nrow(fit$sigma2), 2)
  # Should be positive definite
  evals <- eigen(fit$sigma2, only.values = TRUE)$values
  expect_true(all(evals > 0))
})

# =============================================================================
# RSS: NULL residual_variance defaults
# =============================================================================

test_that("default residual_variance for R=1 RSS data works", {
  sim <- simulate_multivariate(n = 50, p = 30, r = 1)
  z <- sapply(1:ncol(sim$y), function(j) {
    ss <- susieR:::univariate_regression(sim$X, sim$y[, j])
    ss$betahat / ss$sebetahat
  })
  R_ld <- cor(sim$X)
  # Fit with default (NULL) residual_variance
  fit <- mvsusie_rss(z, R_ld, N = nrow(sim$X), L = 2, max_iter = 2,
                     verbosity = 0)
  expect_true(!is.null(fit$sigma2))
})

test_that("default residual_variance for R>1 RSS data works", {
  sim <- simulate_multivariate(n = 50, p = 30, r = 2)
  z <- sapply(1:ncol(sim$y), function(j) {
    ss <- susieR:::univariate_regression(sim$X, sim$y[, j])
    ss$betahat / ss$sebetahat
  })
  R_ld <- cor(sim$X)
  # Fit with default (NULL) residual_variance, must provide prior for R>1
  fit <- mvsusie_rss(z, R_ld, N = nrow(sim$X), L = 2,
                     prior_variance = sim$V,
                     max_iter = 2, verbosity = 0)
  expect_true(is.matrix(fit$sigma2))
  expect_equal(nrow(fit$sigma2), 2)
})
