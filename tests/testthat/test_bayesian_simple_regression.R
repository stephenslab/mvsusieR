context("Test Bayesian multiple regression")

test_that("mmbr is identical to susieR", with(simulate_univariate(), {
    # Test fixed prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = NULL)
    B = BayesianSimpleRegression$new(d$n_effect, V)
    B$fit(d$clone(T), prior_weights = NULL)
    expect_equal(A$mu, as.vector(B$posterior_b1))
    expect_equal(A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf, as.vector(B$lbf))
    # Test estimated prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = "optim")
    B = BayesianSimpleRegression$new(d$n_effect, V)
    B$fit(d$clone(T), prior_weights = NULL, estimate_prior_variance_method='optim')
    expect_equal(A$V, B$prior_variance)
    expect_equal(A$mu, as.vector(B$posterior_b1))
    expect_equal(A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf, as.vector(B$lbf))
}))

test_that("mmbr_RSS is identical to susieR_RSS", with(simulate_univariate(summary=T), {
  # Test fixed prior
  attr(R, 'eigen') = eigen(R, symmetric = T)
  R = susieR:::set_R_attributes(R, 1e-08)
  attr(R, 'lambda') = 0
  Sigma = susieR:::update_Sigma(R, 1, z)
  A = susieR:::single_effect_regression_rss(z, Sigma, V, prior_weights = NULL, optimize_V = "none")
  B = BayesianSimpleRegression$new(d$n_effect, V)
  B$fit(d$clone(T), prior_weights = NULL)
  expect_equal(A$mu, as.vector(B$posterior_b1))
  expect_equal(A$mu2, as.vector(B$posterior_b2))
  expect_equal(A$lbf, as.vector(B$lbf))
  # Test estimated prior
  A = susieR:::single_effect_regression_rss(z, Sigma, V, prior_weights = NULL, optimize_V = "optim")
  B = BayesianSimpleRegression$new(d$n_effect, V)
  B$fit(d$clone(T), prior_weights = NULL, estimate_prior_variance_method='optim')
  expect_equal(A$V, B$prior_variance, tol = 1e-6)
  expect_equal(A$mu, as.vector(B$posterior_b1))
  expect_equal(A$mu2, as.vector(B$posterior_b2))
  expect_equal(A$lbf, as.vector(B$lbf), tol = 1e-6)
}))
