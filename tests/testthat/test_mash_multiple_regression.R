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
    mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight, alpha = 0)
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

test_that("Mash regression + precomputed cov is identical to not precompute", with(simulate_multivariate(r=2), {
    # Run univariate BMR
    prior_var = V[1,1]
    data = DenseData$new(X,y)
    null_weight = 0
    residual_covar = cov(y)
    #
    A_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight, alpha = 1)
    A = MashMultipleRegression$new(ncol(X), residual_covar, A_init)
    A$fit(data, save_summary_stats = T)
    B_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight, alpha = 1)
    B_init$precompute_cov_matrices(data, residual_covar, algorithm = 'cpp')
    B = MashMultipleRegression$new(ncol(X), residual_covar, B_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$sbhat, B$sbhat)
    expect_equal(A$posterior_b1, B$posterior_b1)
    expect_equal(A$posterior_b2, B$posterior_b2)
    expect_equal(A$lbf, B$lbf)
    #
    A_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight, alpha = 0)
    A = MashMultipleRegression$new(ncol(X), residual_covar, A_init)
    A$fit(data, save_summary_stats = T)
    B_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight, alpha = 0)
    B_init$precompute_cov_matrices(data, residual_covar, algorithm = 'cpp')
    B = MashMultipleRegression$new(ncol(X), residual_covar, B_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$sbhat, B$sbhat)
    expect_equal(A$posterior_b1, B$posterior_b1)
    expect_equal(A$posterior_b2, B$posterior_b2)
    expect_equal(A$lbf, B$lbf)
}))