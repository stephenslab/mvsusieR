context("Test summary statistics computation with missing data")

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    residual_variance = cov(y)
    V = cov2cor(residual_variance)
    # V = diag(ncol(y))
    # mix-up X and don't scale it such that sbhat will be different
    X = matrix(runif(ncol(X) * nrow(X)), nrow(X), ncol(X))
    # set one y to be zero and missing code to zero to trick data object into thinking and using
    # missing value computations, for a test.
    d1 = DenseData$new(X,y,scale=F,missing_code='test')
    d2 = DenseData$new(X,y,scale=F)
    # use alpha = 0
    a = 0
    # code for missing data
    expect_equal(d1$Y_has_missing, TRUE)
    res1 = d1$get_sumstats(diag(residual_variance), V, a)
    # code for complete data regular computation
    expect_equal(d2$Y_has_missing, FALSE)
    res2 = d2$get_sumstats(diag(residual_variance), V, a)
    expect_equal(res1$sbhat0, res2$sbhat0)
    expect_equal(res1$svs, res2$svs)
    expect_equal(res1$is_common_sbhat, ifelse(a, T, F))
    expect_equal(res2$is_common_sbhat, ifelse(a, T, F))
    # use alpha = 1
    a = 1
    # code for missing data
    expect_equal(d1$Y_has_missing, TRUE)
    res1 = d1$get_sumstats(diag(residual_variance), V, a)
    # code for complete data regular computation
    expect_equal(d2$Y_has_missing, FALSE)
    res2 = d2$get_sumstats(diag(residual_variance), V, a)
    expect_equal(res1$sbhat0, res2$sbhat0)
    expect_equal(res1$svs, res2$svs)
    expect_equal(res1$is_common_sbhat, ifelse(a, T, F))
    expect_equal(res2$is_common_sbhat, ifelse(a, T, F))
}))