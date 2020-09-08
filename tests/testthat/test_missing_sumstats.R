context("Test missing data computation")

test_that("summary statistics are consistent in DenseData and DenseDataYMissing objects", with(simulate_multivariate(r=3), {
    residual_variance = cov(y)
    # mix-up X and don't scale it such that sbhat will be different
    X = matrix(runif(ncol(X) * nrow(X)), nrow(X), ncol(X))
    # missing value computations, for a test.
    # set `missing_code` to NULL to force using missing data computation routines
    d1 = expect_warning(DenseDataYMissing$new(X,y))
    d1$set_residual_variance(residual_variance, quantities = 'residual_variance')
    d1$standardize(TRUE,FALSE)
    d1$set_residual_variance(quantities = 'effect_variance')
    d2 = DenseData$new(X,y)
    d2$set_residual_variance(residual_variance, quantities = 'residual_variance')
    d2$standardize(TRUE,FALSE)
    d2$set_residual_variance(quantities = 'effect_variance')
    # for missing data
    expect_equal(d1$Y_has_missing, TRUE)
    # for complete data regular computation
    expect_equal(d2$Y_has_missing, FALSE)
    expect_equal(d1$svs, d2$svs)
    expect_equal(d1$is_common_cov, F)
    expect_equal(d2$is_common_cov, F)
    
    d3 = expect_warning(DenseDataYMissing$new(X,y,approximate=FALSE))
    d3$set_residual_variance(diag(diag(residual_variance)), quantities = 'residual_variance')
    d3$standardize(TRUE,FALSE)
    d3$set_residual_variance(quantities = 'effect_variance')
    d4 = expect_warning(DenseDataYMissing$new(X,y,approximate=TRUE))
    d4$set_residual_variance(diag(diag(residual_variance)), quantities = 'residual_variance')
    d4$standardize(TRUE,FALSE)
    d4$set_residual_variance(quantities = 'effect_variance')
    d5 = DenseData$new(X,y)
    d5$set_residual_variance(diag(diag(residual_variance)), quantities = 'residual_variance')
    d5$standardize(TRUE,FALSE)
    d5$set_residual_variance(quantities = 'effect_variance')
    # for missing data
    expect_equal(d3$Y_has_missing, TRUE)
    expect_equal(d4$Y_has_missing, TRUE)
    # for complete data regular computation
    expect_equal(d5$Y_has_missing, FALSE)
    expect_equal(d3$svs, d4$svs)
    expect_equal(d4$svs, d5$svs)
    expect_equal(d3$XtY, d4$XtY)
    expect_equal(d3$is_common_cov, F)
    expect_equal(d4$is_common_cov, F)
    expect_equal(d5$is_common_cov, F)
}))

test_that("When R = 1, result from missing data is same as the result using observed data", with(simulate_multivariate(r=1, center_scale = F, y_missing = 0.5),{
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  data1 = DenseData$new(X[!is.na(y_missing),],y_missing[!is.na(y_missing),,drop=F])
  data1$set_residual_variance(residual_var, quantities = 'residual_variance')
  data1$standardize(TRUE,FALSE)
  data1$set_residual_variance(quantities = 'effect_variance')
  fit1 = BayesianSimpleRegression$new(ncol(X), prior_var)
  fit1$fit(data1, save_summary_stats = T)
  
  data2 = DenseDataYMissing$new(X,y_missing)
  data2$set_residual_variance(residual_var, quantities = 'residual_variance')
  data2$standardize(TRUE,FALSE)
  data2$set_residual_variance(quantities = 'effect_variance')
  fit2 = BayesianSimpleRegression$new(ncol(X), prior_var)
  fit2$fit(data2, save_summary_stats = T)
  
  data3 = DenseDataYMissing$new(X,y_missing,approximate=TRUE)
  data3$set_residual_variance(residual_var, quantities = 'residual_variance')
  data3$standardize(TRUE,FALSE)
  data3$set_residual_variance(quantities = 'effect_variance')
  fit3 = BayesianSimpleRegression$new(ncol(X), prior_var)
  fit3$fit(data3, save_summary_stats = T)
  
  expect_equal(fit1$bhat, fit2$bhat)
  expect_equal(fit1$bhat, fit3$bhat)
  expect_equal(fit1$sbhat, fit2$sbhat)
  expect_equal(fit1$sbhat, fit3$sbhat)
  expect_equal(fit1$posterior_b1, fit2$posterior_b1)
  expect_equal(fit1$posterior_b1, fit3$posterior_b1)
  expect_equal(fit1$posterior_b2, fit2$posterior_b2)
  expect_equal(fit1$posterior_b2, fit3$posterior_b2)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(fit1$lbf, fit3$lbf)
}))

test_that("With full observations, the results are same for DenseDataYMissing and DenseData.", with(simulate_multivariate(r=3, center_scale = F),{
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  data1 = DenseData$new(X, y)
  data1$set_residual_variance(residual_var, quantities = 'residual_variance')
  data1$standardize(TRUE, FALSE)
  data1$set_residual_variance(quantities = 'effect_variance')
  fit1 = BayesianMultivariateRegression$new(ncol(X), prior_var)
  fit1$fit(data1, save_summary_stats = T)

  data2 = expect_warning(DenseDataYMissing$new(X,y))
  data2$set_residual_variance(residual_var, quantities = 'residual_variance')
  data2$standardize(TRUE,FALSE)
  data2$set_residual_variance(quantities = 'effect_variance')
  fit2 = BayesianMultivariateRegression$new(ncol(X), prior_var)
  fit2$fit(data2, save_summary_stats = T)
  
  data3 = expect_warning(DenseDataYMissing$new(X,y,approximate=TRUE))
  data3$set_residual_variance(residual_var, quantities = 'residual_variance')
  data3$standardize(TRUE,FALSE)
  data3$set_residual_variance(quantities = 'effect_variance')
  fit3 = BayesianMultivariateRegression$new(ncol(X), prior_var)
  fit3$fit(data3, save_summary_stats = T)
  
  expect_equal(fit1$bhat, fit2$bhat)
  expect_equal(fit1$bhat, fit3$bhat)
  expect_equal(fit1$sbhat, fit2$sbhat)
  expect_equal(fit1$sbhat, fit3$sbhat)
  expect_equal(fit1$posterior_b1, fit2$posterior_b1)
  expect_equal(fit1$posterior_b1, fit3$posterior_b1)
  expect_equal(fit1$posterior_b2, fit2$posterior_b2)
  expect_equal(fit1$posterior_b2, fit3$posterior_b2)
  expect_equal(fit1$lbf, fit2$lbf)
  expect_equal(fit1$lbf, fit3$lbf)

  # Mash regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1-null_weight, null_weight)
  mash_init$precompute_cov_matrices(data2, algorithm = 'cpp')
  fit4 = MashRegression$new(ncol(X), mash_init)
  fit4$fit(data2, save_summary_stats = T)
  expect_equal(fit1$bhat, fit4$bhat)
  expect_equal(fit1$posterior_b1, fit4$posterior_b1)
  expect_equal(fit1$posterior_b2, fit4$posterior_b2, check.attributes = FALSE)
  expect_equal(fit1$lbf, fit4$lbf)
  
  mash_init$precompute_cov_matrices(data3, algorithm = 'cpp')
  fit5 = MashRegression$new(ncol(X), mash_init)
  fit5$fit(data3, save_summary_stats = T)
  expect_equal(fit1$bhat, fit5$bhat)
  expect_equal(fit1$posterior_b1, fit5$posterior_b1)
  expect_equal(fit1$posterior_b2, fit5$posterior_b2, check.attributes = FALSE)
  expect_equal(fit1$lbf, fit5$lbf)
}))

test_that("When R = 1, estimated prior variance with missing data agrees with full data", with(simulate_multivariate(r=1, center_scale = F, y_missing = 0.5), {
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  fit1 = msusie(X[!is.na(y_missing),],y_missing[!is.na(y_missing),,drop=F], L = L,
              prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
              intercept=T, standardize = F,
              estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  fit2 = msusie(X, y_missing, L=L,
              prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
              intercept=T, standardize = F, 
              estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  fit3 = msusie(X, y_missing, L=L,
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
                intercept=T, standardize = F, 
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM',
                approximate=TRUE)
  
  expect_equal(fit1$alpha, fit2$alpha, tolerance = 1E-8)
  expect_equal(fit2$alpha, fit3$alpha)
  expect_equal(fit1$lbf, fit2$lbf, tolerance = 1E-8)
  expect_equal(fit2$lbf, fit3$lbf)
  expect_equal(fit1$b1, fit2$b1, tolerance = 1E-8)
  expect_equal(fit2$b1, fit3$b1)
  expect_equal(fit1$b2, fit2$b2, tolerance = 1E-8)
  expect_equal(fit2$b2, fit3$b2)
  expect_equal(fit1$coef, fit2$coef, tolerance = 1E-8)
  expect_equal(fit2$coef, fit3$coef)
  expect_equal(fit1$V, fit2$V, tolerance = 1E-8)
  expect_equal(fit2$V, fit3$V)
}))

test_that("With full observation, the estimated prior variance are same for DenseDataYMissing and DenseData", with(simulate_multivariate(r=3, center_scale = F), {
  # Multivariate regression
  prior_var = V
  residual_var = cov(y)
  fit1 = msusie(X, y, L = L, 
                prior_variance=prior_var, residual_variance = residual_var, compute_objective=F, 
                intercept=T, standardize = F, 
                estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
  
  data2 = expect_warning(DenseDataYMissing$new(X,y))
  data2$set_residual_variance(residual_var, quantities = 'residual_variance')
  data2$standardize(TRUE,FALSE)
  data2$set_residual_variance(quantities = 'effect_variance')
  fit2 = mmbr_core(data2, s_init=NULL, L=L, prior_variance=prior_var, prior_weights=c(rep(1/ncol(X), ncol(X))),
                   estimate_residual_variance=F, estimate_prior_variance=T, estimate_prior_method='EM', check_null_threshold=0,
                   precompute_covariances=F, compute_objective=F, max_iter=100, tol=1e-3, track_fit=F, verbose=T, n_thread=1)
  
  data3 = expect_warning(DenseDataYMissing$new(X,y,approximate=TRUE))
  data3$set_residual_variance(residual_var, quantities = 'residual_variance')
  data3$standardize(TRUE,FALSE)
  data3$set_residual_variance(quantities = 'effect_variance')
  fit3 = mmbr_core(data3, s_init=NULL, L=L, prior_variance=prior_var, prior_weights=c(rep(1/ncol(X), ncol(X))),
                   estimate_residual_variance=F, estimate_prior_variance=T, estimate_prior_method='EM', check_null_threshold=0,
                   precompute_covariances=F, compute_objective=F, max_iter=100, tol=1e-3, track_fit=F, verbose=T, n_thread=1)
  
  expect_equal(fit1$alpha, fit2$alpha, tolerance = 1E-8)
  expect_equal(fit2$alpha, fit3$alpha)
  expect_equal(fit1$lbf, fit2$lbf, tolerance = 1E-8)
  expect_equal(fit2$lbf, fit3$lbf)
  expect_equal(fit1$b1, fit2$b1, tolerance = 1E-8)
  expect_equal(fit2$b1, fit3$b1)
  expect_equal(fit1$b2, fit2$b2, tolerance = 1E-8)
  expect_equal(fit2$b2, fit3$b2)
  expect_equal(fit1$coef, fit2$coef, tolerance = 1E-8, check.attributes = FALSE)
  expect_equal(fit2$coef, fit3$coef, check.attributes = FALSE)
  expect_equal(fit1$V, fit2$V, tolerance = 1E-8)
  expect_equal(fit2$V, fit3$V)
  
  # Mash regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1-null_weight, null_weight)
  fit4 = expect_warning(mmbr_core(data2, s_init=NULL, L=L, prior_variance=mash_init, prior_weights=c(rep(1/ncol(X), ncol(X))),
                   estimate_residual_variance=F, estimate_prior_variance=T, estimate_prior_method='EM', check_null_threshold=0,
                   precompute_covariances=F, compute_objective=F, max_iter=100, tol=1e-3, track_fit=F, verbose=T, n_thread=1))
  
  fit5 = expect_warning(mmbr_core(data3, s_init=NULL, L=L, prior_variance=mash_init, prior_weights=c(rep(1/ncol(X), ncol(X))),
                                  estimate_residual_variance=F, estimate_prior_variance=T, estimate_prior_method='EM', check_null_threshold=0,
                                  precompute_covariances=F, compute_objective=F, max_iter=100, tol=1e-3, track_fit=F, verbose=T, n_thread=1))
  
  expect_equal(fit1$alpha, fit4$alpha, tolerance = 1E-8)
  expect_equal(fit1$lbf, fit4$lbf, tolerance = 1E-8)
  expect_equal(fit1$b1, fit4$b1, tolerance = 1E-8)
  expect_equal(fit1$b2, fit4$b2, tolerance = 1E-8)
  expect_equal(fit1$coef, fit4$coef, tolerance = 1E-8, check.attributes = FALSE)
  expect_equal(fit1$V, fit4$V, tolerance = 1E-8)
  
  expect_equal(fit1$alpha, fit5$alpha, tolerance = 1E-8)
  expect_equal(fit1$lbf, fit5$lbf, tolerance = 1E-8)
  expect_equal(fit1$b1, fit5$b1, tolerance = 1E-8)
  expect_equal(fit1$b2, fit5$b2, tolerance = 1E-8)
  expect_equal(fit1$coef, fit5$coef, tolerance = 1E-8, check.attributes = FALSE)
  expect_equal(fit1$V, fit5$V, tolerance = 1E-8)
}))

