context("Test missing data computation")

# NOTE: Tests using R6 DenseData/DenseDataYMissing/BayesianSimpleRegression/
# BayesianMultivariateRegression/mvsusie_core have been removed because R6
# classes have been deleted. The remaining tests exercise the S3 path
# through mvsusie() which handles missing data internally.

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
  expect_equal(fit1$b1, fit2$b1, tolerance = 1E-8)
  expect_equal(fit1$b2, fit2$b2, tolerance = 1E-8)
  expect_equal(fit1$coef, fit2$coef, tolerance = 1E-8)
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

  expect_equal(fit4_matrix$alpha,fit4_mash$alpha,tolerance = 1e-8,check.attributes = FALSE)
  expect_equal(fit4_matrix$lbf,fit4_mash$lbf,tolerance = 1e-8,check.attributes = FALSE)
  expect_equal(fit4_matrix$b1,fit4_mash$b1,tolerance = 1e-8,check.attributes = FALSE)
  expect_equal(fit4_matrix$b2,fit4_mash$b2,tolerance = 1e-8,check.attributes = FALSE)
  expect_equal(fit4_matrix$coef,fit4_mash$coef,tolerance = 1e-8,check.attributes = FALSE)
  expect_equal(fit4_matrix$V,fit4_mash$V,tolerance = 1e-8,check.attributes = FALSE)
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
