context("Test summary statistics computation with missing data")

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    residual_variance = cov(y)
    V = cov2cor(residual_variance)
    # V = diag(ncol(y))
    # mix-up X and don't scale it such that sbhat will be different
    X = matrix(runif(ncol(X) * nrow(X)), nrow(X), ncol(X))
    # missing value computations, for a test.
    # set `missing_code` to NULL to force using missing data computation routines
    d1 = DenseData$new(X,y,scale=F,missing_code=NULL)
    d2 = DenseData$new(X,y,scale=F)
    # for missing data
    expect_equal(d1$Y_has_missing, TRUE)
    res1 = d1$get_sumstats(diag(residual_variance), V)
    # for complete data regular computation
    expect_equal(d2$Y_has_missing, FALSE)
    res2 = d2$get_sumstats(diag(residual_variance), V)
    expect_equal(res1$sbhat0, res2$sbhat0)
    expect_equal(res1$svs, res2$svs)
    expect_equal(res1$is_common_sbhat, F)
    expect_equal(res2$is_common_sbhat, F)
    # now test it when X is scaled
    d1 = DenseData$new(X,y,missing_code=NULL)
    d2 = DenseData$new(X,y,)
    # for missing data
    expect_equal(d1$Y_has_missing, TRUE)
    res1 = d1$get_sumstats(diag(residual_variance), V)
    # for complete data regular computation
    expect_equal(d2$Y_has_missing, FALSE)
    res2 = d2$get_sumstats(diag(residual_variance), V)
    expect_equal(res1$sbhat0, res2$sbhat0)
    expect_equal(res1$svs, res2$svs)
    expect_equal(res1$is_common_sbhat, T)
    expect_equal(res2$is_common_sbhat, T)
}))

test_that("marginal statistics agree between mmbr and R's `.lm.fit` when there is missing data", with(simulate_multivariate(r=3,center_scale=F,y_missing=0.5), {
   resid_Y <- compute_cov_diag(y)
   resid_Y_miss <- compute_cov_diag(y_missing)
   # can only compare bhat not sbhat
   # full data, center and scale
   d = DenseData$new(X, y,scale=T)
   res = d$get_sumstats(resid_Y)
   b = res$bhat
   for (i in 1:ncol(b)) {
       expect_equal(b[,i], susieR:::univariate_regression(X, y[,i], scale=T)$betahat)
   }
   # missing data center and scale
   d = DenseData$new(X, y_missing, scale=T)
   res = d$get_sumstats(resid_Y_miss)
   b = res$bhat
   for (i in 1:ncol(b)) {
       expect_equal(b[,i], susieR:::univariate_regression(X, y_missing[,i], scale=T)$betahat)
   }
   # full data not scale
   d = DenseData$new(X, y, scale=F)
   res = d$get_sumstats(resid_Y)
   b = res$bhat
   for (i in 1:ncol(b)) {
       expect_equal(b[,i], susieR:::univariate_regression(X, y[,i], scale=F)$betahat)
   }
   # missing data not scale
   d = DenseData$new(X, y_missing, scale=F)
   res = d$get_sumstats(resid_Y_miss)
   b = res$bhat
   for (i in 1:ncol(b)) {
       expect_equal(b[,i], susieR:::univariate_regression(X, y_missing[,i], scale=F)$betahat)
   }
}))