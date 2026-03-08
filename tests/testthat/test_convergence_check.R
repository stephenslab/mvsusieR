context("Test ELBO convergence")

test_that("mvsusie converges with monotonic ELBO (dense)", {
  set.seed(1)
  n <- 100; p <- 80; L <- 5
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  V <- 0.2 * var(y)

  fit <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    tol = 1e-6, verbosity = 0)
  expect_true(all(diff(fit$elbo) >= -1e-6))

  fit2 <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    tol = 1e-6, verbosity = 0)
  expect_true(all(diff(fit2$elbo) >= -1e-6))
})

test_that("mvsusie converges with monotonic ELBO (suff stat)", {
  set.seed(1)
  n <- 100; p <- 80; L <- 5
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  XtX <- crossprod(Xc)
  Xty <- crossprod(Xc, yc)
  yty <- crossprod(yc)
  V <- 0.2 * as.numeric(yty / (n - 1))

  fit <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    tol = 1e-6, verbosity = 0)
  expect_true(all(diff(fit$elbo) >= -1e-6))

  fit2 <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    tol = 1e-6, verbosity = 0)
  expect_true(all(diff(fit2$elbo) >= -1e-6))
})
