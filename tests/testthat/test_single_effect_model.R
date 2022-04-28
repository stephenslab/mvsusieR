context("Test Single Effect regression")

test_that("mvsusieR is identical to susieR", with(simulate_univariate(), {
    # Test fixed prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = NULL)
    kl = susieR:::SER_posterior_e_loglik(X,y,1,A$alpha*A$mu,A$alpha*A$mu2)- A$loglik
    B = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    d.copy = d$clone(T)
    B$fit(d.copy)
    B$compute_kl(d.copy)
    expect_equal(A$alpha * A$mu, as.vector(B$posterior_b1))
    expect_equal(A$alpha * A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf_model, B$lbf)
    expect_equal(kl, B$kl)
    # Test estimated prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = "optim")
    kl = susieR:::SER_posterior_e_loglik(X,y,1,A$alpha*A$mu,A$alpha*A$mu2)- A$loglik
    B = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    d.copy = d$clone(T)
    B$fit(d.copy, estimate_prior_variance_method='optim')
    B$compute_kl(d.copy)
    expect_equal(A$alpha * A$mu, as.vector(B$posterior_b1))
    expect_equal(A$alpha * A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf_model, B$lbf)
    expect_equal(kl, B$kl)
}))

test_that("mvsusieR is identical to susieR (RSS)", with(simulate_univariate(summary = T), {
  # Test fixed prior
  adj = (n-1)/(z^2 + n - 2)
  ztilde = sqrt(adj) * z
  A = susieR:::single_effect_regression_ss(ztilde * sqrt(n-1), (n-1)*diag(R), V, prior_weights = NULL, optimize_V = "none")
  kl = susieR:::SER_posterior_e_loglik_ss((n-1)*diag(R), ztilde * sqrt(n-1), 1, A$alpha*A$mu,A$alpha*A$mu2)- A$lbf_model
  B = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  d.copy = d$clone(T)
  B$fit(d.copy)
  B$compute_kl(d.copy)
  expect_equal(A$alpha * A$mu, as.vector(B$posterior_b1))
  expect_equal(A$alpha * A$mu2, as.vector(B$posterior_b2))
  expect_equal(A$lbf_model, B$lbf)
  expect_equal(kl, B$kl)
  # Test estimated prior
  A = susieR:::single_effect_regression_ss(ztilde * sqrt(n-1), (n-1)*diag(R), V, prior_weights = NULL, optimize_V = "optim")
  kl = susieR:::SER_posterior_e_loglik_ss((n-1)*diag(R),sqrt(n-1)*ztilde,1,A$alpha*A$mu,A$alpha*A$mu2)- A$lbf_model
  B = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  d.copy = d$clone(T)
  B$fit(d.copy, estimate_prior_variance_method='optim')
  B$compute_kl(d.copy)
  expect_equal(A$alpha * A$mu, as.vector(B$posterior_b1))
  expect_equal(A$alpha * A$mu2, as.vector(B$posterior_b2))
  expect_equal(A$lbf_model, B$lbf)
  expect_equal(kl, B$kl)
}))
