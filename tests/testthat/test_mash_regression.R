context("Test MASH regression")

test_that("Degenerated mash regression is identical to univariate BR", with(simulate_multivariate(r=1), {
    # Compare matrix prior path vs mash prior path with L=1
    residual_var = cov(y)
    null_weight = 0
    mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
    # Matrix prior (S3 path)
    A = mvsusie(X, y, L=1, prior_variance=V,
                residual_variance=residual_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    # Mash prior path (S3 path)
    B = mvsusie(X, y, L=1, prior_variance=mash_init,
                residual_variance=residual_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    expect_susie_equal(A, B, F, F)
}))

test_that("Single component mash regression is identical to multivariate BR", with(simulate_multivariate(r=3), {
    # Run multivariate regression (matrix prior)
    residual_var = cov(y)
    # Run MASH regression EE model with single component
    mash_init = create_mash_prior(Ulist = list(V), grid = 1)
    A = mvsusie(X, y, L=1, prior_variance=V,
                residual_variance=residual_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    B = mvsusie(X, y, L=1, prior_variance=mash_init,
                residual_variance=residual_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    expect_susie_equal(A, B, F, F)
}))

test_that("Mash regression + precomputed cov is identical to not precompute", with(simulate_multivariate(r=2), {
    residual_covar = cov(y)
    null_weight = 0
    mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
    A = mvsusie(X, y, L=1, prior_variance=mash_init,
                residual_variance=residual_covar,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE,
                precompute_covariances=FALSE)
    B = mvsusie(X, y, L=1, prior_variance=mash_init,
                residual_variance=residual_covar,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE,
                precompute_covariances=TRUE)
    expect_susie_equal(A, B, F, F)
}))

test_that("Degenerated mash regression is identical to univariate BR RSS", with(simulate_multivariate(r=1), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R_ld = cor(X)
  n = nrow(X)
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V*n), grid = 1, null_weight = null_weight)
  # Matrix prior (S3 path)
  A = mvsusie_rss(z, R_ld, N=n, L=1, prior_variance=V*n,
                  estimate_prior_variance=FALSE,
                  compute_objective=FALSE)
  # Mash prior path (S3 path)
  B = mvsusie_rss(z, R_ld, N=n, L=1, prior_variance=mash_init,
                  estimate_prior_variance=FALSE,
                  compute_objective=FALSE)
  expect_susie_equal(A, B, F, F)
}))

test_that("Single component mash regression is identical to multivariate BR RSS", with(simulate_multivariate(r=3), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
    })
  R = cor(X)
  n = nrow(X)
  residual_var = cov(y)
  # Run MASH regression EE model with single component
  mash_init = create_mash_prior(Ulist = list(V), grid = 1)
  A = mvsusie_rss(z, R, N=n, L=1, prior_variance=V,
                  estimate_prior_method="EM",
                  compute_objective=FALSE)
  B = mvsusie_rss(z, R, N=n, L=1, prior_variance=mash_init,
                  estimate_prior_method="EM",
                  compute_objective=FALSE)
  expect_susie_equal(A, B, F, F)
}))

test_that("Mash regression + precomputed cov is identical to not precompute (RSS)", with(simulate_multivariate(r=2), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  n = nrow(X)
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  A = mvsusie_rss(z, R, N=n, L=1, prior_variance=mash_init,
                  estimate_prior_method="EM",
                  compute_objective=FALSE,
                  precompute_covariances=FALSE)
  B = mvsusie_rss(z, R, N=n, L=1, prior_variance=mash_init,
                  estimate_prior_method="EM",
                  compute_objective=FALSE,
                  precompute_covariances=TRUE)
  expect_susie_equal(A, B, F, F)
}))
