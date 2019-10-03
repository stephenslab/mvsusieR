context("Test summary statistics computation with missing data")

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    d = DenseData$new(X,y)
    residual_variance = cov(y)
    V = cov2cor(residual_variance)
    # V = diag(ncol(y))
    # use alpha = 0
    a = 0
    # code for missing data
    res = get_sumstats_missing_data(d$X, d$Y, diag(residual_variance), V, a)
    svs1 = res$svs
    sbhat1 = res$sbhat0
    # code for regular computation
    sigma2 = diag(residual_variance)
    sbhat0 = sqrt(do.call(rbind, lapply(1:length(d$d), function(j) sigma2 / d$d[j])))
    sbhat = sbhat0 ^ (1-a)
    common_sbhat = is_mat_common(sbhat)
    if (common_sbhat) svs2 = sbhat[1,] * t(V * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
    else svs2 = lapply(1:nrow(sbhat), function(j) sbhat[j,] * t(V * sbhat[j,]))
    expect_equal(sbhat1, sbhat0)
    expect_equal(svs1, svs2)
    expect_equal(common_sbhat, ifelse(a, T, F))
    # use alpha = 1
    a = 1
    # code for missing data
    res = get_sumstats_missing_data(d$X, d$Y, diag(residual_variance), V, a)
    svs1 = res$svs
    sbhat1 = res$sbhat0
    # code for regular computation
    sigma2 = diag(residual_variance)
    sbhat0 = sqrt(do.call(rbind, lapply(1:length(d$d), function(j) sigma2 / d$d[j])))
    sbhat = sbhat0 ^ (1-a)
    common_sbhat = is_mat_common(sbhat)
    if (common_sbhat) svs2 = sbhat[1,] * t(V * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
    else svs2 = lapply(1:nrow(sbhat), function(j) sbhat[j,] * t(V * sbhat[j,]))
    expect_equal(sbhat1, sbhat0)
    expect_equal(svs1, svs2)
    expect_equal(common_sbhat, ifelse(a, T, F))
}))