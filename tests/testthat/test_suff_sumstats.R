# NOTE: Tests using R6 SSData/DenseData/BayesianSimpleRegression/
# BayesianMultivariateRegression have been removed because R6 classes
# have been deleted. The remaining tests compare mvsusie() vs
# mvsusie_ss() through the S3 path.

# =============================================================================
# Default residual_variance: individual vs SS paths must agree
# =============================================================================

test_that("R=1 individual vs SS agree with default residual_variance", with(simulate_multivariate(r = 1, center_scale = FALSE), {
  X_colmeans <- colMeans(X)
  Y_colmeans <- colMeans(y)
  Xc <- t(t(X) - X_colmeans)
  yc <- t(t(y) - Y_colmeans)
  XtX <- crossprod(Xc)
  XtY <- crossprod(Xc, yc)
  YtY <- crossprod(yc)
  prior_var <- V[1, 1]
  # Both use default residual_variance = NULL
  fit_ind <- mvsusie(X, y, L = L, prior_variance = prior_var,
                     intercept = TRUE, standardize = TRUE,
                     estimate_residual_variance = FALSE,
                     estimate_prior_variance = FALSE, verbose = FALSE)
  fit_ss <- mvsusie_ss(XtX, XtY, YtY, n, L = L,
                       X_colmeans = X_colmeans, Y_colmeans = Y_colmeans,
                       prior_variance = prior_var, standardize = TRUE,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE, verbose = FALSE)
  expect_equal(fit_ind$alpha, fit_ss$alpha)
  expect_equal(fit_ind$elbo, fit_ss$elbo)
  expect_equal(get_b1(fit_ind), get_b1(fit_ss))
  expect_equal(coef(fit_ind), coef(fit_ss))
  expect_equal(fit_ind$V, fit_ss$V)
}))

test_that("R>1 individual vs SS agree with default residual_variance", with(simulate_multivariate(r = 3, center_scale = FALSE), {
  X_colmeans <- colMeans(X)
  Y_colmeans <- colMeans(y)
  Xc <- t(t(X) - X_colmeans)
  yc <- t(t(y) - Y_colmeans)
  XtX <- crossprod(Xc)
  XtY <- crossprod(Xc, yc)
  YtY <- crossprod(yc)
  # Both use default residual_variance = NULL
  fit_ind <- mvsusie(X, y, L = L, prior_variance = V,
                     intercept = TRUE, standardize = TRUE,
                     estimate_residual_variance = FALSE,
                     estimate_prior_variance = FALSE, verbose = FALSE)
  fit_ss <- mvsusie_ss(XtX, XtY, YtY, n, L = L,
                       X_colmeans = X_colmeans, Y_colmeans = Y_colmeans,
                       prior_variance = V, standardize = TRUE,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = FALSE, verbose = FALSE)
  expect_equal(fit_ind$alpha, fit_ss$alpha)
  expect_equal(fit_ind$elbo, fit_ss$elbo)
  expect_equal(get_b1(fit_ind), get_b1(fit_ss))
  expect_equal(coef(fit_ind), coef(fit_ss))
  expect_equal(fit_ind$V, fit_ss$V)
}))

test_that("R>1 individual vs SS agree with optim and default residual_variance", with(simulate_multivariate(r = 2, center_scale = FALSE), {
  X_colmeans <- colMeans(X)
  Y_colmeans <- colMeans(y)
  Xc <- t(t(X) - X_colmeans)
  yc <- t(t(y) - Y_colmeans)
  XtX <- crossprod(Xc)
  XtY <- crossprod(Xc, yc)
  YtY <- crossprod(yc)
  # Both use default residual_variance = NULL, with prior estimation
  fit_ind <- mvsusie(X, y, L = 5, prior_variance = V,
                     intercept = TRUE, standardize = TRUE,
                     estimate_residual_variance = FALSE,
                     estimate_prior_variance = TRUE,
                     estimate_prior_method = "optim",
                     max_iter = 20, verbose = FALSE)
  fit_ss <- mvsusie_ss(XtX, XtY, YtY, n, L = 5,
                       X_colmeans = X_colmeans, Y_colmeans = Y_colmeans,
                       prior_variance = V, standardize = TRUE,
                       estimate_residual_variance = FALSE,
                       estimate_prior_variance = TRUE,
                       estimate_prior_method = "optim",
                       max_iter = 20, verbose = FALSE)
  expect_equal(fit_ind$alpha, fit_ss$alpha)
  expect_equal(fit_ind$elbo, fit_ss$elbo)
  expect_equal(get_b1(fit_ind), get_b1(fit_ss))
  expect_equal(coef(fit_ind), coef(fit_ss))
  # V is estimated by iterative optimization (optim); individual data uses
  # crossprod(X, residual) while SS uses XtX %*% b, which accumulate different
  # floating point errors over 20 iterations. Use slightly relaxed tolerance.
  expect_equal(fit_ind$V, fit_ss$V, tolerance = 1e-7)
}))

# =============================================================================
# Existing tests with explicit residual_variance
# =============================================================================

test_that("With full observations, matrix vs mash prior agree (suff stat)", with(simulate_multivariate(r=3, center_scale = F),{
  prior_var = V
  residual_var = cov(y)
  X_colmeans = colMeans(X)
  Y_colmeans = colMeans(y)
  X.c = t(t(X) - X_colmeans)
  y.c = t(t(y) - Y_colmeans)
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)

  # Mash regression: compare matrix prior vs mash prior through mvsusie_ss
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit3_matrix = mvsusie_ss(XtX, XtY, YtY, n, L=1, X_colmeans, Y_colmeans,
                          prior_variance=prior_var, residual_variance = residual_var,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=FALSE)
  fit3_mash = mvsusie_ss(XtX, XtY, YtY, n, L=1, X_colmeans, Y_colmeans,
                          prior_variance=mash_init, residual_variance = residual_var,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=FALSE,
                          precompute_cache = T)
  expect_susie_equal(fit3_matrix, fit3_mash, F, F)
}))

test_that("When R = 1, estimated prior variance with ss data agrees with full data", with(simulate_multivariate(r=1, center_scale = F), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = mvsusie(X,y, L = L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  X_colmeans = colMeans(X)
  Y_colmeans = colMeans(y)
  X.c = t(t(X) - X_colmeans)
  y.c = t(t(y) - Y_colmeans)
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = mvsusie_ss(XtX, XtY, YtY, n, L=L, X_colmeans, Y_colmeans,
                prior_variance=prior_var, residual_variance = residual_var,
                standardize = TRUE,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(get_b1(fit1), get_b1(fit2))
  expect_equal(get_b2(fit1), get_b2(fit2))
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(get_b1(fit1),get_b1(fit2))
  expect_equal(fit1$V, fit2$V)
}))

test_that("With full observation, the estimated prior variance are same for SSData and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = mvsusie(X, y, L = L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  X_colmeans = colMeans(X)
  Y_colmeans = colMeans(y)
  X.c = t(t(X) - X_colmeans)
  y.c = t(t(y) - Y_colmeans)
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = mvsusie_ss(XtX, XtY, YtY, n, L=L, X_colmeans, Y_colmeans,
                          prior_variance=prior_var, residual_variance = residual_var,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(get_b1(fit1), get_b1(fit2))
  expect_equal(get_b2(fit1), get_b2(fit2))
  expect_equal(coef(fit1), coef(fit2))
  expect_equal(get_b1(fit1), get_b1(fit2))
  expect_equal(fit1$V, fit2$V)

  # Mash regression
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit3 = mvsusie_ss(XtX, XtY, YtY, n, L=L, X_colmeans, Y_colmeans,
                          prior_variance=mash_init, residual_variance = residual_var,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit3$alpha)
  expect_equal(fit1$lbf, fit3$lbf)
  expect_equal(get_b1(fit1), get_b1(fit3))
  expect_equal(get_b2(fit1), get_b2(fit3))
  expect_equal(coef(fit1), coef(fit3))
  expect_equal(get_b1(fit1), get_b1(fit3))
  expect_equal(fit1$V, fit3$V)
}))

test_that("When R = 1, the elbo using sufficient data agrees with full data", with(simulate_multivariate(r=1, center_scale = F), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = mvsusie(X,y, L = L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  X_colmeans = colMeans(X)
  Y_colmeans = colMeans(y)
  X.c = t(t(X) - X_colmeans)
  y.c = t(t(y) - Y_colmeans)
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = mvsusie_ss(XtX, XtY, YtY, n, L=L,X_colmeans, Y_colmeans,
                prior_variance=prior_var, residual_variance = residual_var,
                standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  expect_equal(fit1$elbo, fit2$elbo)
}))

test_that("With full observation, the elbo are same for SSData and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = mvsusie(X, y, L = L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  X_colmeans = colMeans(X)
  Y_colmeans = colMeans(y)
  X.c = t(t(X) - X_colmeans)
  y.c = t(t(y) - Y_colmeans)
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = mvsusie_ss(XtX, XtY, YtY, n, L=L,X_colmeans, Y_colmeans,
                          prior_variance=prior_var, residual_variance = residual_var,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=F)
  # Mash regression
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit3 = mvsusie_ss(XtX, XtY, YtY, n, L=L,X_colmeans, Y_colmeans,
                          prior_variance=mash_init, residual_variance = residual_var,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=F,
                          precompute_cache = T)

  expect_equal(fit1$elbo, fit2$elbo)
  expect_equal(fit1$elbo, fit3$elbo)
}))

# ---- Tests merged from test_convergence_check.R ----

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
    tol = 1e-6, verbose = FALSE)
  expect_true(all(diff(fit$elbo) >= -1e-6))

  fit2 <- mvsusie_core(X, matrix(y, ncol = 1), L = L,
    prior_variance = V,
    residual_variance = matrix(1, 1, 1),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    tol = 1e-6, verbose = FALSE)
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
    tol = 1e-6, verbose = FALSE)
  expect_true(all(diff(fit$elbo) >= -1e-6))

  fit2 <- mvsusie_ss_core(XtX = XtX, XtY = Xty, YtY = yty, N = n, L = L,
    prior_variance = V,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    estimate_prior_method = "optim",
    tol = 1e-6, verbose = FALSE)
  expect_true(all(diff(fit2$elbo) >= -1e-6))
})

# ---- Tests merged from test_missing_sumstats.R ----

test_that("With full observations, matrix vs mash prior agree (with missing data path)", with(simulate_multivariate(r=3, center_scale = F),{
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)

  # Mash regression: compare matrix prior vs mash prior through mvsusie with L=1
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit4_matrix = mvsusie(X, y, L=1, prior_variance=prior_var,
                        residual_variance=residual_var,
                        estimate_residual_variance=FALSE,
                        estimate_prior_variance=FALSE,
                        )
  fit4_mash = mvsusie(X, y, L=1, prior_variance=mash_init,
                      residual_variance=residual_var,
                      estimate_residual_variance=FALSE,
                      estimate_prior_variance=FALSE,
                      precompute_cache=TRUE)
  expect_susie_equal(fit4_matrix, fit4_mash, F, F)
}))

test_that("When R = 1, estimated prior variance with missing data agrees with full data", with(simulate_multivariate(r=1, center_scale = F, y_missing = 0.5), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = mvsusie(X[!is.na(y_missing),],y_missing[!is.na(y_missing),,drop=F], L = L,
              prior_variance=prior_var, residual_variance = residual_var,
              intercept=T, standardize = T,
              estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  fit2 = mvsusie(X, y_missing, L=L,
              prior_variance=prior_var, residual_variance = residual_var,
              intercept=T, standardize = T,
              estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit2$alpha, tolerance = 1E-8)
  expect_equal(fit1$lbf, fit2$lbf, tolerance = 1E-8)
  expect_equal(get_b1(fit1), get_b1(fit2), tolerance = 1E-8)
  expect_equal(get_b2(fit1), get_b2(fit2), tolerance = 1E-8)
  expect_equal(coef(fit1), coef(fit2), tolerance = 1E-8)
  expect_equal(fit1$V, fit2$V, tolerance = 1E-8)
}))

test_that("With full observation, the estimated prior variance: matrix vs mash prior agree", with(simulate_multivariate(r=3, center_scale = F), {
  # Mash regression: compare S3 matrix prior vs S3 mash prior through mvsusie
  prior_var = V
  residual_var = cov(y)
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit4_matrix = mvsusie(X, y, L=L, prior_variance=V,
                        residual_variance=residual_var,
                        intercept=T, standardize=T,
                        estimate_residual_variance=F, estimate_prior_variance=T,
                        estimate_prior_method='EM')
  fit4_mash = mvsusie(X, y, L=L, prior_variance=mash_init,
                      residual_variance=residual_var,
                      intercept=T, standardize=T,
                      estimate_residual_variance=F, estimate_prior_variance=T,
                      estimate_prior_method='EM')

  expect_equal(fit4_matrix$alpha,fit4_mash$alpha,tolerance = 1e-8,ignore_attr = TRUE)
  expect_equal(fit4_matrix$lbf,fit4_mash$lbf,tolerance = 1e-8,ignore_attr = TRUE)
  expect_equal(get_b1(fit4_matrix),get_b1(fit4_mash),tolerance = 1e-8,ignore_attr = TRUE)
  expect_equal(get_b2(fit4_matrix),get_b2(fit4_mash),tolerance = 1e-8,ignore_attr = TRUE)
  expect_equal(coef(fit4_matrix),coef(fit4_mash),tolerance = 1e-8,ignore_attr = TRUE)
  expect_equal(fit4_matrix$V,fit4_mash$V,tolerance = 1e-8,ignore_attr = TRUE)
}))

test_that("When R = 1, the elbo with missing data agrees with full data", with(simulate_multivariate(r=1, center_scale = F, y_missing = 0.5), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = mvsusie(X[!is.na(y_missing),],y_missing[!is.na(y_missing),,drop=F], L = L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  fit2 = mvsusie(X, y_missing, L=L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  # Compare converged ELBO values (final element); convergence speed
  # may differ between complete-data and missing-data paths.
  expect_equal(tail(fit1$elbo, 1), tail(fit2$elbo, 1), tolerance = 1e-4)
}))

test_that("With full observation, the elbo: matrix vs mash prior agree", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = mvsusie(X, y, L = L,
                prior_variance=prior_var, residual_variance = residual_var,
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  # Mash regression: compare matrix prior vs mash prior through mvsusie
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit4 = mvsusie(X, y, L=L, prior_variance=mash_init,
                 residual_variance=residual_var,
                 intercept=T, standardize=T,
                 estimate_residual_variance=F, estimate_prior_variance=F)

  expect_equal(fit1$elbo, fit4$elbo)
}))

# =============================================================================
# RSS with Bhat/Shat/varY: slope coefficients match individual data exactly
# =============================================================================

test_that("R>1 mvsusie_rss with Bhat/Shat/varY gives same slopes as mvsusie", {
  set.seed(42)
  n <- 200; p <- 50; R <- 2
  X <- matrix(rnorm(n * p), n, p)
  B <- matrix(0, p, R)
  B[1, ] <- c(1, 0.5)
  B[2, ] <- c(0.5, 1)
  Y <- X %*% B + matrix(rnorm(n * R), n, R)

  prior_var <- diag(R) * 0.2
  resid_var <- cov(Y)

  # Individual data fit
  fit_ind <- mvsusie(X, Y, L = 5, prior_variance = prior_var,
                     residual_variance = resid_var,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     max_iter = 100, verbose = FALSE)

  # Compute Bhat/Shat
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  Bhat <- matrix(0, p, R)
  Shat <- matrix(0, p, R)
  for (r in seq_len(R)) {
    reg <- susieR:::univariate_regression(Xc, Yc[, r])
    Bhat[, r] <- reg$betahat
    Shat[, r] <- reg$sebetahat
  }

  # RSS fit with Bhat/Shat/varY
  fit_rss <- mvsusie_rss(Bhat = Bhat, Shat = Shat, R = cor(X),
                          N = n, varY = cov(Y),
                          L = 5, prior_variance = prior_var,
                          residual_variance = resid_var,
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          max_iter = 100, verbose = FALSE)

  # Slope coefficients should match to machine precision
  coef_ind <- coef(fit_ind)
  coef_rss <- coef(fit_rss)
  expect_equal(coef_ind[-1, ], coef_rss[-1, ], tolerance = 1e-10)
  expect_equal(fit_ind$alpha, fit_rss$alpha, tolerance = 1e-10)
  expect_equal(fit_ind$pip, fit_rss$pip, tolerance = 1e-10)
})

test_that("R=1 mvsusie_rss with Bhat/Shat/varY gives same slopes as mvsusie", {
  set.seed(42)
  n <- 200; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  b <- rep(0, p); b[1:2] <- 1
  y <- X %*% b + rnorm(n)

  prior_var <- 0.2
  resid_var <- var(as.vector(y))

  fit_ind <- mvsusie(X, y, L = 5, prior_variance = prior_var,
                     residual_variance = resid_var,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     max_iter = 100, verbose = FALSE)

  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  reg <- susieR:::univariate_regression(Xc, yc)

  fit_rss <- mvsusie_rss(Bhat = reg$betahat, Shat = reg$sebetahat,
                          R = cor(X), N = n, varY = var(as.vector(y)),
                          L = 5, prior_variance = prior_var,
                          residual_variance = resid_var,
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          max_iter = 100, verbose = FALSE)

  coef_ind <- coef(fit_ind)
  coef_rss <- coef(fit_rss)
  expect_equal(coef_ind[-1, ], coef_rss[-1, ], tolerance = 1e-10)
  expect_equal(fit_ind$alpha, fit_rss$alpha, tolerance = 1e-10)
  expect_equal(fit_ind$pip, fit_rss$pip, tolerance = 1e-10)
})

test_that("R>1 mvsusie_rss with z-scores only gives same PIPs as mvsusie", {
  set.seed(42)
  n <- 200; p <- 50; R <- 2
  X <- matrix(rnorm(n * p), n, p)
  B <- matrix(0, p, R)
  B[1, ] <- c(1, 0.5)
  B[2, ] <- c(0.5, 1)
  Y <- X %*% B + matrix(rnorm(n * R), n, R)

  prior_var <- diag(R) * 0.2
  resid_var <- cov(Y)

  fit_ind <- mvsusie(X, Y, L = 5, prior_variance = prior_var,
                     residual_variance = resid_var,
                     estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE,
                     max_iter = 100, verbose = FALSE)

  Z <- calc_z(X, Y, center = TRUE, scale = TRUE)
  LD <- cor(X)

  fit_rss <- mvsusie_rss(Z, LD, N = n,
                          L = 5, prior_variance = prior_var,
                          residual_variance = resid_var,
                          estimate_prior_variance = FALSE,
                          estimate_residual_variance = FALSE,
                          max_iter = 100, verbose = FALSE)

  # When only z-scores and LD are provided (without varY), the RSS path
  # uses cov2cor(residual_variance) as YtY, which is an approximation.
  # PIPs should be highly correlated but not necessarily identical.
  expect_gt(cor(fit_ind$pip, fit_rss$pip), 0.95)
  # Slopes are on different scales: standardized vs original.
  # Check high correlation for non-trivial coefficients.
  slopes_ind <- as.vector(coef(fit_ind)[-1, ])
  slopes_rss <- as.vector(coef(fit_rss)[-1, ])
  expect_gt(cor(slopes_ind, slopes_rss), 0.95)
})
