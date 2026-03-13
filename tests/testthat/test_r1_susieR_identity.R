# =============================================================================
# R=1 identity tests: mvsusie S3 dispatch must match susieR::susie
#
# When R=1, the multivariate S3 path (mv_individual dispatch through
# susieR workhorse) must produce numerically identical results to
# susieR::susie (individual dispatch). This validates that every
# S3 method in mv_methods.R correctly reduces to the univariate case.
#
# Parameter mapping:
#   susieR::susie(X, y, scaled_prior_variance = V/var(y), residual_variance = 1)
#   mvsusie_core(X, Y, prior_variance = V, residual_variance = matrix(1,1,1))
#
# The small numerical differences (~1e-5) arise from:
#   - dmvnorm with 1x1 matrices vs dnorm-based BF in susieR
#   - Matrix algebra with 1x1 matrices vs scalar arithmetic
#
# NOTE: Residual variance estimation (estimate_residual_variance = TRUE) is
# not yet properly unified between mvsusie and susieR, so all tests here
# use fixed residual variance. Residual variance estimation tests will be
# added once that feature is properly implemented for mvsusie.
# =============================================================================


# Helper: run mvsusie_core with R=1 and compare against susieR::susie
compare_r1_fits <- function(fit_susie, fit_mv, tol = 1e-4,
                            check_prior = FALSE,
                            label = "") {
  info <- function(field) paste(label, "-", field)

  # alpha (L x J)
  expect_equal(fit_susie$alpha, unname(fit_mv$alpha),
               tolerance = tol, info = info("alpha"))

  # mu (L x J)
  expect_equal(fit_susie$mu, unname(fit_mv$mu),
               tolerance = tol, info = info("mu"))

  # b1 = alpha * mu (L x J)
  expect_equal(fit_susie$alpha * fit_susie$mu, unname(fit_mv$b1),
               tolerance = tol, info = info("b1"))

  # b2 = alpha * mu2 (L x J)
  expect_equal(fit_susie$alpha * fit_susie$mu2, unname(fit_mv$b2),
               tolerance = tol, info = info("b2"))

  # KL divergence (length L)
  if (!any(is.na(fit_susie$KL)) && !any(is.na(fit_mv$KL)))
    expect_equal(unname(fit_susie$KL), unname(fit_mv$KL),
                 tolerance = tol, info = info("KL"))

  # lbf (length L)
  if (!any(is.na(fit_susie$lbf)) && !any(is.na(fit_mv$lbf)))
    expect_equal(unname(fit_susie$lbf), unname(fit_mv$lbf),
                 tolerance = tol, info = info("lbf"))

  # lbf_variable (L x J)
  if (!any(is.na(fit_susie$lbf_variable)) && !any(is.na(fit_mv$lbf_variable)))
    expect_equal(fit_susie$lbf_variable, unname(fit_mv$lbf_variable),
                 tolerance = tol, info = info("lbf_variable"))

  # ELBO
  if (!any(is.na(fit_susie$elbo)) && !any(is.na(fit_mv$elbo)))
    expect_equal(fit_susie$elbo, fit_mv$elbo,
                 tolerance = tol, info = info("elbo"))

  # Number of iterations
  expect_equal(fit_susie$niter, fit_mv$niter,
               info = info("niter"))

  # pip
  expect_equal(unname(fit_susie$pip), unname(fit_mv$pip),
               tolerance = tol, info = info("pip"))

  # fitted values
  expect_equal(as.vector(fit_susie$fitted), as.vector(fit_mv$fitted),
               tolerance = tol, info = info("fitted"))

  # coef (intercept + J coefficients)
  expect_equal(as.vector(coef(fit_susie)), as.vector(fit_mv$coef),
               tolerance = tol, info = info("coef"))

  # Prior variance (V)
  # susieR V is a length-L vector; mvsusie V is now also a numeric vector
  if (check_prior) {
    expect_equal(fit_susie$V, fit_mv$V, tolerance = tol, info = info("V"))
  }
}


test_that("R=1: fixed prior, fixed residual", {
  set.seed(1)
  n <- 100; p <- 80; L <- 5
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  y <- y - mean(y)
  V <- 0.2 * var(y)

  fit_susie <- susieR::susie(X, y, L = L,
    scaled_prior_variance = V / var(y),
    residual_variance = 1,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE)

  fit_mv <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    verbose = FALSE)

  compare_r1_fits(fit_susie, fit_mv,
                  label = "fixed prior fixed residual")
})


test_that("R=1: estimated prior (optim), fixed residual", {
  set.seed(1)
  n <- 100; p <- 80; L <- 5
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  y <- y - mean(y)
  V <- 0.2 * var(y)

  fit_susie <- susieR::susie(X, y, L = L,
    scaled_prior_variance = V / var(y),
    residual_variance = 1,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim")

  fit_mv <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    verbose = FALSE)

  compare_r1_fits(fit_susie, fit_mv,
                  check_prior = TRUE,
                  label = "estimated prior (optim) fixed residual")
})


test_that("R=1: estimated prior (EM), fixed residual", {
  set.seed(1)
  n <- 100; p <- 80; L <- 5
  beta <- rep(0, p); beta[1:4] <- 1
  X <- matrix(rnorm(n * p, 3, 4), n, p)
  y <- c(X %*% beta + rnorm(n))
  y <- y - mean(y)
  V <- 0.2 * var(y)

  fit_susie <- susieR::susie(X, y, L = L,
    scaled_prior_variance = V / var(y),
    residual_variance = 1,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "EM")

  fit_mv <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "EM",
    verbose = FALSE)

  compare_r1_fits(fit_susie, fit_mv,
                  check_prior = TRUE,
                  label = "estimated prior (EM) fixed residual")
})


test_that("R=1: sufficient statistics path matches susieR::susie_ss", {
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

  # susieR suff stat
  fit_ss <- susieR::susie_ss(XtX = XtX, Xty = as.vector(Xty),
    yty = as.numeric(yty), n = n, L = L,
    scaled_prior_variance = V / (as.numeric(yty) / (n - 1)),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE)

  # mvsusie suff stat S3
  fit_mv <- mvsusie_ss_core(
    XtX = XtX, XtY = Xty,
    YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    verbose = FALSE)

  # Compare key fields
  expect_equal(fit_ss$alpha, unname(fit_mv$alpha),
               tolerance = 1e-4, info = "ss alpha")
  expect_equal(fit_ss$alpha * fit_ss$mu, unname(fit_mv$b1),
               tolerance = 1e-4, info = "ss b1")
  if (!any(is.na(fit_ss$elbo)) && !any(is.na(fit_mv$elbo)))
    expect_equal(fit_ss$elbo, fit_mv$elbo,
                 tolerance = 1e-4, info = "ss elbo")
  expect_equal(fit_ss$niter, fit_mv$niter, info = "ss niter")
  expect_equal(unname(fit_ss$pip), unname(fit_mv$pip),
               tolerance = 1e-4, info = "ss pip")
  if (!any(is.na(fit_ss$KL)) && !any(is.na(fit_mv$KL)))
    expect_equal(fit_ss$KL, fit_mv$KL,
                 tolerance = 1e-4, info = "ss KL")
})


test_that("R=1: suff stat with estimated prior (EM), fixed residual", {
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

  fit_ss <- susieR::susie_ss(XtX = XtX, Xty = as.vector(Xty),
    yty = as.numeric(yty), n = n, L = L,
    scaled_prior_variance = V / (as.numeric(yty) / (n - 1)),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "EM")

  fit_mv <- mvsusie_ss_core(
    XtX = XtX, XtY = Xty,
    YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "EM",
    verbose = FALSE)

  expect_equal(fit_ss$alpha, unname(fit_mv$alpha),
               tolerance = 1e-4, info = "ss EM alpha")
  expect_equal(fit_ss$alpha * fit_ss$mu, unname(fit_mv$b1),
               tolerance = 1e-4, info = "ss EM b1")
  if (!any(is.na(fit_ss$elbo)) && !any(is.na(fit_mv$elbo)))
    expect_equal(fit_ss$elbo, fit_mv$elbo,
                 tolerance = 1e-4, info = "ss EM elbo")
  expect_equal(fit_ss$niter, fit_mv$niter, info = "ss EM niter")
  expect_equal(unname(fit_ss$pip), unname(fit_mv$pip),
               tolerance = 1e-4, info = "ss EM pip")
})
