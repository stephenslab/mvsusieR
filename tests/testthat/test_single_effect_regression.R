context("Test Single Effect regression")

test_that("mmbr is identical to susieR", with(simulate_univariate(), {
    # Test fixed prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = NULL)
    kl = susieR:::SER_posterior_e_loglik(X,y,1,A$alpha*A$mu,A$alpha*A$mu2)- A$loglik
    B = SingleEffectRegression(BayesianMultipleRegression)$new(d$n_effect, 1, V, estimate_prior_variance = FALSE, prior_weights = NULL)
    d.copy = d$clone(T)
    B$fit(d.copy)
    B$compute_kl(d.copy)
    expect_equal(A$alpha * A$mu, as.vector(B$posterior_b1))
    expect_equal(A$alpha * A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf_model, B$lbf_single_effect)
    expect_equal(kl, B$kl)
    # Test estimated prior
    A = susieR:::single_effect_regression(y, X, V, residual_variance = 1, prior_weights = NULL, optimize_V = "optim")
    kl = susieR:::SER_posterior_e_loglik(X,y,1,A$alpha*A$mu,A$alpha*A$mu2)- A$loglik
    B = SingleEffectRegression(BayesianMultipleRegression)$new(d$n_effect, 1, V, estimate_prior_variance = TRUE, prior_weights = NULL)
    d.copy = d$clone(T)
    B$fit(d.copy)
    B$compute_kl(d.copy)
    expect_equal(A$alpha * A$mu, as.vector(B$posterior_b1))
    expect_equal(A$alpha * A$mu2, as.vector(B$posterior_b2))
    expect_equal(A$lbf_model, B$lbf_single_effect)
    expect_equal(kl, B$kl)
}))