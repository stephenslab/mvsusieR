context("Test Estimate Residual Variance")

test_that("estimated residual variance in multivariate case is identical to univariate case", with(simulate_multivariate(r = 1), {
  # fit susie
  A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = as.numeric(cov(y)), prior_weights = NULL, 
                    estimate_residual_variance = T, estimate_prior_variance = FALSE)
  d = DenseData$new(X,y)
  
  # BayesianSimpleRegression
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, as.numeric(cov(y)), as.numeric(V))
  B = SuSiE$new(SER, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  B$fit(d.copy)
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A, BA, F, F)
  expect_equal(A$sigma2, BA$sigma2)
  
  # BayesianMultivariateRegression
  SER = SingleEffectModel(BayesianMultivariateRegression)$new(d$n_effect, cov(y), V)
  C = SuSiE$new(SER, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  C$fit(d.copy)
  CA = report_susie_model(d.copy, C)
  expect_susieR_equal(A, CA, F, F)
  
  # Mash Regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
  residual_covar = cov(y)
  SER = SingleEffectModel(MashRegression)$new(ncol(X), residual_covar, mash_init)
  D = SuSiE$new(SER, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  D$fit(d.copy)
  # compare result
  expect_equal(B$posterior_b1, D$posterior_b1)
  expect_equal(B$posterior_b2, D$posterior_b2)
  expect_equal(B$lbf, D$lbf)
  expect_equal(B$residual_variance, D$residual_variance)
}))
