context("Test missing data computation")

test_that("summary statistics are consistent in DenseData and SSData objects", with(simulate_multivariate(r=3), {
  residual_variance = cov(y)
  # mix-up X and don't scale it such that sbhat will be different
  X = matrix(runif(ncol(X) * nrow(X)), nrow(X), ncol(X))
  # missing value computations, for a test.
  # set `missing_code` to NULL to force using missing data computation routines
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  N = n
  d1 = SSData$new(XtX,XtY,YtY,N)
  d1$set_residual_variance(residual_variance, quantities = 'residual_variance')
  d1$standardize(FALSE)
  d1$set_residual_variance(quantities = 'effect_variance')
  d2 = DenseData$new(X,y)
  d2$set_residual_variance(residual_variance, quantities = 'residual_variance')
  d2$standardize(TRUE,FALSE)
  d2$set_residual_variance(quantities = 'effect_variance')

  expect_equal(d1$svs, d2$svs)
  expect_equal(d1$is_common_cov, F)
  expect_equal(d2$is_common_cov, F)
}))

test_that("When R = 1, result from ss data is same as the result using full data", with(simulate_multivariate(r=1, center_scale = F),{
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  data1 = DenseData$new(X,y)
  data1$set_residual_variance(residual_var, quantities = 'residual_variance')
  data1$standardize(TRUE,TRUE)
  data1$set_residual_variance(quantities = 'effect_variance')
  fit1 = BayesianSimpleRegression$new(ncol(X), prior_var)
  fit1$fit(data1, save_summary_stats = T)
  
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  data2 = SSData$new(XtX, XtY, YtY, n)
  data2$set_residual_variance(residual_var, quantities = 'residual_variance')
  data2$standardize(TRUE)
  data2$set_residual_variance(quantities = 'effect_variance')
  fit2 = BayesianSimpleRegression$new(ncol(X), prior_var)
  fit2$fit(data2, save_summary_stats = T)
  
  expect_equal(fit1$bhat, fit2$bhat)
  expect_equal(fit1$sbhat, fit2$sbhat)
  expect_equal(fit1$posterior_b1, fit2$posterior_b1)
  expect_equal(fit1$posterior_b2, fit2$posterior_b2)
  expect_equal(fit1$lbf, fit2$lbf)
}))

test_that("With full observations, the results are same for SSData and DenseData.", with(simulate_multivariate(r=3, center_scale = F),{
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  data1 = DenseData$new(X, y)
  data1$set_residual_variance(residual_var, quantities = 'residual_variance')
  data1$standardize(TRUE, TRUE)
  data1$set_residual_variance(quantities = 'effect_variance')
  fit1 = BayesianMultivariateRegression$new(ncol(X), prior_var)
  fit1$fit(data1, save_summary_stats = T)
  
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  data2 = SSData$new(XtX, XtY, YtY, n)
  data2$set_residual_variance(residual_var, quantities = 'residual_variance')
  data2$standardize(TRUE)
  data2$set_residual_variance(quantities = 'effect_variance')
  fit2 = BayesianMultivariateRegression$new(ncol(X), prior_var)
  fit2$fit(data2, save_summary_stats = T)
  
  expect_equal(fit1$bhat, fit2$bhat)
  expect_equal(fit1$sbhat, fit2$sbhat)
  expect_equal(fit1$posterior_b1, fit2$posterior_b1)
  expect_equal(fit1$posterior_b2, fit2$posterior_b2)
  expect_equal(fit1$lbf, fit2$lbf)
  
  # Mash regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1-null_weight, null_weight)
  mash_init$precompute_cov_matrices(data2, algorithm = 'cpp')
  fit3 = MashRegression$new(ncol(X), mash_init)
  fit3$fit(data2, save_summary_stats = T)
  expect_equal(fit1$bhat, fit3$bhat)
  expect_equal(fit1$posterior_b1, fit3$posterior_b1)
  expect_equal(fit1$posterior_b2, fit3$posterior_b2, check.attributes = FALSE)
  expect_equal(fit1$lbf, fit3$lbf)
}))

test_that("When R = 1, estimated prior variance with ss data agrees with full data", with(simulate_multivariate(r=1, center_scale = F), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = msusie(X,y, L = L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = msusie_suff_stat(XtX, XtY, YtY, n, L=L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
                standardize = T, 
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(fit1$b1, fit2$b1)
  expect_equal(fit1$b2, fit2$b2)
  expect_equal(fit1$coef[-1], fit2$coef[-1])
  expect_equal(fit1$V, fit2$V)
}))

test_that("With full observation, the estimated prior variance are same for SSData and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = msusie(X, y, L = L, 
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
                intercept=T, standardize = T, 
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = msusie_suff_stat(XtX, XtY, YtY, n, L=L,
                          prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
                          standardize = T, 
                          estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  expect_equal(fit1$alpha, fit2$alpha)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(fit1$b1, fit2$b1)
  expect_equal(fit1$b2, fit2$b2)
  expect_equal(fit1$coef[-1,], fit2$coef[-1,])
  expect_equal(fit1$V, fit2$V)
  
  # Mash regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1-null_weight, null_weight)
  fit3 = msusie_suff_stat(XtX, XtY, YtY, n, L=L,
                          prior_variance=mash_init, residual_variance = residual_var, compute_objective=F, 
                          standardize = T, 
                          estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM',
                          precompute_covariances=T)
  
  expect_equal(fit1$alpha, fit3$alpha)
  expect_equal(fit1$lbf, fit3$lbf)
  expect_equal(fit1$b1, fit3$b1)
  expect_equal(fit1$b2, fit3$b2)
  expect_equal(fit1$coef[-1,], fit3$coef[-1,])
  expect_equal(fit1$V, fit3$V)
}))

test_that("When R = 1, the elbo using sufficient data agrees with full data", with(simulate_multivariate(r=1, center_scale = F), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = msusie(X,y, L = L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=T, 
                intercept=T, standardize = T,
                estimate_residual_variance=F, estimate_prior_variance=F)
  
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = msusie_suff_stat(XtX, XtY, YtY, n, L=L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=T, 
                standardize = T, 
                estimate_residual_variance=F, estimate_prior_variance=F)
  
  expect_equal(fit1$elbo, fit2$elbo)
}))

test_that("With full observation, the elbo are same for SSData and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = msusie(X, y, L = L, 
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=T, 
                intercept=T, standardize = T, 
                estimate_residual_variance=F, estimate_prior_variance=F)
  
  X.c = t(t(X) - colMeans(X))
  y.c = t(t(y) - colMeans(y))
  XtY = crossprod(X.c, y.c)
  XtX = crossprod(X.c)
  YtY = crossprod(y.c)
  fit2 = msusie_suff_stat(XtX, XtY, YtY, n, L=L,
                          prior_variance=prior_var, residual_variance = residual_var, compute_objective=T, 
                          standardize = T, 
                          estimate_residual_variance=F, estimate_prior_variance=F)
  # Mash regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1-null_weight, null_weight)
  fit3 = msusie_suff_stat(XtX, XtY, YtY, n, L=L,
                          prior_variance=mash_init, residual_variance = residual_var, compute_objective=T, 
                          standardize = T, 
                          estimate_residual_variance=F, estimate_prior_variance=F,
                          precompute_covariances = T)
  
  expect_equal(fit1$elbo, fit2$elbo)
  expect_equal(fit1$elbo, fit3$elbo)
}))

