context("Test MASH multiple regression")

test_that("Degenerated mash regression is identical to univariate BMR", with(simulate_multivariate(r=1), {
    # Run univariate BMR
    prior_var = V[1,1]
    residual_var = as.numeric(var(y))
    data = DenseData$new(X,y)
    A = BayesianMultipleRegression$new(ncol(X), residual_var, prior_var)
    A$fit(data, save_summary_stats = T)
    # Run MASH BMR
    # null_weight = 1 - 1 / ncol(X)
    null_weight = 0
    mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
    residual_covar = cov(y)
    B = MashMultipleRegression$new(ncol(X), residual_covar, mash_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$sbhat, B$sbhat[,1])
    expect_equal(A$posterior_b1, B$posterior_b1)
    expect_equal(A$posterior_b2[,1], as.vector(B$posterior_b2))
    expect_equal(A$lbf, B$lbf)
}))