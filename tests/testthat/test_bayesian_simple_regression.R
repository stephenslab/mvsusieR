context("Test Bayesian multiple regression")

test_that("mmbr is identical to susieR", with(simulate_univariate(), {
    # Test fixed prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = NULL)
    B = BayesianSimpleRegression$new(d$n_effect, 1, V, estimate_prior_variance = FALSE)
    B$fit(d$clone(T), prior_weights = NULL)
    expect_equal(A$mu, as.vector(B$posterior_b1))
    expect_equal(A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf, as.vector(B$lbf))
    # Test estimated prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = "optim")
    B = BayesianSimpleRegression$new(d$n_effect, 1, V, estimate_prior_variance = TRUE)
    B$fit(d$clone(T), prior_weights = NULL)
    expect_equal(A$V, B$prior_variance)
    expect_equal(A$mu, as.vector(B$posterior_b1))
    expect_equal(A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf, as.vector(B$lbf))
}))