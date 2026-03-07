context("Test SuSiE regression")

# NOTE: Tests comparing R6 SuSiE/SingleEffectModel/BayesianSimpleRegression
# against susieR have been removed because R6 classes have been deleted.
# Those tests are superseded by test_r1_susieR_identity.R which compares
# the S3 path against susieR for R=1.

test_that("mash regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    # Compare matrix prior vs mash prior (both S3 path)
    residual_var = cov(y)
    m_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    A = mvsusie(X,y,L=L,prior_variance=V,
                residual_variance=residual_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    B = mvsusie(X,y,L=L,prior_variance=m_init,
                residual_variance=residual_var,
                compute_objective=FALSE, estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    expect_susie_equal(A,B,F,F)
}))

test_that("multivariate regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    # Compare matrix prior vs mash prior (both S3 path) with precomputed cov
    residual_var = cov(y)
    m_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    A = mvsusie(X,y,L=L,prior_variance=V,
                residual_variance=residual_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE, precompute_covariances=TRUE)
    B = mvsusie(X,y,L=L,prior_variance=m_init,
                residual_variance=residual_var,
                compute_objective=FALSE, estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    m_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
    A = mvsusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE,
                              estimate_residual_variance=FALSE, estimate_prior_variance=FALSE,
                              precompute_covariances=FALSE)
    B = mvsusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE,
               precompute_covariances=TRUE, estimate_prior_variance=FALSE, estimate_residual_variance=FALSE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with BMR using one component prior matrix", with(simulate_multivariate(r=3), {
    m_init = create_mixture_prior(mixture_prior = list(matrices=list(V)))
    # don't compare ELBO
    A = mvsusie(X,y,L=L,prior_variance=m_init, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    B = mvsusie(X,y,L=L,prior_variance=V, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE)
    expect_susie_equal(A,B,F,F)
    # compare ELBO
    A = mvsusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    B = mvsusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=FALSE)
    expect_susie_equal(A,B,F,F)
    # compare estimate prior variance "EM" method (mash_prior vs matrix prior)
    A = mvsusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM', precompute_covariances=TRUE)
    B = mvsusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
    expect_susie_equal(A,B,F,F)
}))

test_that("customized initialization interface", with(simulate_multivariate(r=3), {
    # not sure what to test here ...
    m_init = create_mixture_prior(mixture_prior = list(matrices = list(V), weights = 1), null_weight=0)
    A = mvsusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE, estimate_prior_variance=FALSE)
    B = mvsusie(X,y,L=L,prior_variance=m_init,s_init=A,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE, estimate_prior_variance=FALSE)
    # let's just test of null is null ...
    null_weight = 0.2
    m_init = create_mixture_prior(R = ncol(y),null_weight = null_weight, max_mixture_len=-1)
    expect_equal(m_init$null_weight, null_weight)
}))

test_that("mash regression in SuSiE is identical to univariate case (RSS)", with(simulate_multivariate(r=1), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R_ld = cor(X)
  n = nrow(X)
  m_init = create_mash_prior(Ulist = list(V*n), grid = 1, null_weight = 0)
  # Compare matrix prior vs mash prior (both S3 path)
  A = mvsusie_rss(z,R_ld,N=n,L=L,prior_variance=V*n,
             estimate_prior_variance=FALSE,
             compute_objective=FALSE)
  B = mvsusie_rss(z,R_ld,N=n,L=L,prior_variance=m_init,compute_objective=FALSE,estimate_prior_variance=FALSE,
                 precompute_covariances=TRUE)
  expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed (RSS)", with(simulate_multivariate(r=3), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  n = nrow(X)
  m_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = 0)
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE, precompute_covariances=FALSE, estimate_prior_method="EM")
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE, precompute_covariances=TRUE, estimate_prior_method="EM")
  expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with BMR using one component prior matrix (RSS)", with(simulate_multivariate(r=3), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  n = nrow(X)
  m_init = create_mixture_prior(mixture_prior = list(matrices=list(V)))
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init, estimate_prior_method="EM")
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=V, estimate_prior_method="EM")
  expect_susie_equal(A,B,F,F)
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init, compute_objective=T, estimate_prior_method="EM")
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=V, compute_objective=T, estimate_prior_method="EM")
  expect_susie_equal(A,B,F,F)
}))

test_that("customized initialization interface (RSS)", with(simulate_multivariate(r=3), {
  # not sure what to test here ...
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  n = nrow(X)
  m_init = create_mixture_prior(mixture_prior = list(matrices = list(V), weights = 1), null_weight=0)
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_prior_method="EM")
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,s_init=A,compute_objective=FALSE, estimate_prior_method="EM")
  # let's just test of null is null ...
  null_weight = 0.2
  m_init = create_mixture_prior(R = ncol(y),null_weight = null_weight, max_mixture_len=-1)
  expect_equal(m_init$null_weight, null_weight)
}))
