context("Test Estimate Residual Variance")

test_that("estimated residual variance in multivariate case is identical to univariate case", with(simulate_multivariate(r = 1), {
  # fit susie
  A = susieR::susie(X, y, L = L,
                    scaled_prior_variance = V/var(y),
                    residual_variance = as.numeric(cov(y)), prior_weights = NULL,
                    estimate_residual_variance = T, estimate_prior_variance = FALSE)
  d = DenseData$new(X,y)
  d$standardize(TRUE,TRUE)
  # BayesianSimpleRegression
  SER.simple = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, as.numeric(cov(y)), as.numeric(V))
  B = SuSiE$new(SER.simple, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  B$fit(d.copy)
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A, BA, F, F)
  expect_equal(A$sigma2, BA$sigma2)

  # BayesianMultivariateRegression
  SER.multi = SingleEffectModel(BayesianMultivariateRegression)$new(d$n_effect, cov(y), V)
  C = SuSiE$new(SER.multi, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  C$fit(d.copy)
  CA = report_susie_model(d.copy, C)
  expect_susie_equal(BA, CA, F, F)
  expect_equal(BA$sigma2, as.numeric(CA$sigma2))

  # Mash Regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
  residual_covar = cov(y)
  SER.mash = SingleEffectModel(MashRegression)$new(ncol(X), residual_covar, mash_init)
  D = SuSiE$new(SER.mash, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  D$fit(d.copy)
  DA = report_susie_model(d.copy, D)
  expect_susie_equal(CA, DA, F, T)
}))

test_that("estimated residual variance in multivariate case is identical to univariate case RSS", with(simulate_multivariate(r = 1), {
  # fit susie_rss
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  residual_var = as.numeric(var(y))

  # A = susieR::susie_rss(z, R, L = L, prior_variance = as.numeric(V),
  #                       residual_variance = 1, prior_weights = NULL,
  #                       estimate_residual_variance = T, estimate_prior_variance = FALSE, check_z = F)
  d = RSSData$new(z,R,1e-08)

  # BayesianSimpleRegression
  SER.simple = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, as.numeric(cov(y)), as.numeric(V))
  B = SuSiE$new(SER.simple, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  B$fit(d.copy)
  BA = report_susie_model(d.copy, B)
  # expect_susieR_equal(A, BA, F, F)
  # expect_equal(A$sigma2, BA$sigma2)

  # BayesianMultivariateRegression
  SER.multi = SingleEffectModel(BayesianMultivariateRegression)$new(d$n_effect, cov(y), V)
  C = SuSiE$new(SER.multi, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  C$fit(d.copy)
  CA = report_susie_model(d.copy, C)
  expect_susie_equal(BA, CA, F, F)
  expect_equal(BA$sigma2, as.numeric(CA$sigma2))

  # Mash Regression
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
  residual_covar = cov(y)
  SER.mash = SingleEffectModel(MashRegression)$new(ncol(X), residual_covar, mash_init)
  D = SuSiE$new(SER.mash, L, estimate_residual_variance = T)
  d.copy = d$clone(T)
  D$fit(d.copy)
  DA = report_susie_model(d.copy, D)
  expect_susie_equal(CA, DA, F, T)
}))