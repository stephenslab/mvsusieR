context("Test check for convergence using ELBO or not")

# NOTE: Tests using R6 SingleEffectModel/BayesianSimpleRegression/SuSiE/DenseData
# have been replaced with S3 path equivalents using mvsusie() / mvsusie_ss().

test_that("mvsusie gets same result checking ELBO or not (dense)", {
  set.seed(1)
  n <- 100; p <- 200; L <- 10
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  V <- 0.2 * var(y)

  # Fixed prior
  A <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    compute_objective = TRUE,
    tol = 1e-6, verbosity = 0)
  B <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    compute_objective = FALSE,
    tol = 1e-6, verbosity = 0)
  expect_susie_equal(A, B, FALSE, FALSE, tol = 1e-3)

  # Estimated prior (optim)
  A <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    compute_objective = TRUE,
    tol = 1e-6, verbosity = 0)
  B <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    compute_objective = FALSE,
    tol = 1e-6, verbosity = 0)
  expect_susie_equal(A, B, TRUE, FALSE, tol = 5e-4)
})

test_that("mvsusie gets same result checking ELBO or not (suff stat)", {
  set.seed(1)
  n <- 100; p <- 200; L <- 10
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  XtX <- crossprod(Xc)
  Xty <- crossprod(Xc, yc)
  yty <- crossprod(yc)
  V <- 0.2 * as.numeric(yty / (n - 1))

  # Fixed prior
  A <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    compute_objective = TRUE,
    tol = 1e-6, verbosity = 0)
  B <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    compute_objective = FALSE,
    tol = 1e-6, verbosity = 0)
  expect_susie_equal(A, B, FALSE, FALSE, tol = 1e-3)

  # Estimated prior (optim)
  A <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    compute_objective = TRUE,
    tol = 1e-6, verbosity = 0)
  B <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    compute_objective = FALSE,
    tol = 1e-6, verbosity = 0)
  expect_susie_equal(A, B, TRUE, FALSE, tol = 5e-4)
})
