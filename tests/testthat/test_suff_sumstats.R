context("Test sufficient data computation")

# NOTE: Tests using R6 SSData/DenseData/BayesianSimpleRegression/
# BayesianMultivariateRegression have been removed because R6 classes
# have been deleted. The remaining tests compare mvsusie() vs
# mvsusie_ss() through the S3 path.

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
                          compute_objective=F, standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=FALSE)
  fit3_mash = mvsusie_ss(XtX, XtY, YtY, n, L=1, X_colmeans, Y_colmeans,
                          prior_variance=mash_init, residual_variance = residual_var,
                          compute_objective=F, standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=FALSE,
                          precompute_covariances = T)
  expect_susie_equal(fit3_matrix, fit3_mash, F, F)
}))

test_that("When R = 1, estimated prior variance with ss data agrees with full data", with(simulate_multivariate(r=1, center_scale = F), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = mvsusie(X,y, L = L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F,
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
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F,
                standardize = TRUE,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(fit1$b1, fit2$b1)
  expect_equal(fit1$b2, fit2$b2)
  expect_equal(fit1$coef, fit2$coef)
  expect_equal(fit1$b1_rescaled,fit2$b1_rescaled)
  expect_equal(fit1$V, fit2$V)
}))

test_that("With full observation, the estimated prior variance are same for SSData and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = mvsusie(X, y, L = L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F,
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
                          prior_variance=prior_var, residual_variance = residual_var, compute_objective=F,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(fit1$b1, fit2$b1)
  expect_equal(fit1$b2, fit2$b2)
  expect_equal(fit1$coef, fit2$coef)
  expect_equal(fit1$b1_rescaled, fit2$b1_rescaled)
  expect_equal(fit1$V, fit2$V)

  # Mash regression
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit3 = mvsusie_ss(XtX, XtY, YtY, n, L=L, X_colmeans, Y_colmeans,
                          prior_variance=mash_init, residual_variance = residual_var, compute_objective=F,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')

  expect_equal(fit1$alpha, fit3$alpha)
  expect_equal(fit1$lbf, fit3$lbf)
  expect_equal(fit1$b1, fit3$b1)
  expect_equal(fit1$b2, fit3$b2)
  expect_equal(fit1$coef, fit3$coef)
  expect_equal(fit1$b1_rescaled, fit3$b1_rescaled)
  expect_equal(fit1$V, fit3$V)
}))

test_that("When R = 1, the elbo using sufficient data agrees with full data", with(simulate_multivariate(r=1, center_scale = F), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = mvsusie(X,y, L = L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=T,
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
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=T,
                standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)

  expect_equal(fit1$elbo, fit2$elbo)
}))

test_that("With full observation, the elbo are same for SSData and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = mvsusie(X, y, L = L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=T,
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
                          prior_variance=prior_var, residual_variance = residual_var, compute_objective=T,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=F)
  # Mash regression
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  fit3 = mvsusie_ss(XtX, XtY, YtY, n, L=L,X_colmeans, Y_colmeans,
                          prior_variance=mash_init, residual_variance = residual_var, compute_objective=T,
                          standardize = T,
                          estimate_residual_variance=F, estimate_prior_variance=F,
                          precompute_covariances = T)

  expect_equal(fit1$elbo, fit2$elbo)
  expect_equal(fit1$elbo, fit3$elbo)
}))
