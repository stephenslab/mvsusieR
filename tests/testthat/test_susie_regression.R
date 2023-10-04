context("Test SuSiE regression")

test_that("mvsusieR is identical to susieR when prior is a scalar", with(simulate_univariate(), {
    # Test fixed prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = FALSE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(TRUE)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A,BA,F,F,tol = 1e-4)
    # Test estimated prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = TRUE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(TRUE)
    B$fit(d.copy,estimate_prior_variance_method='optim')
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A,BA,TRUE,FALSE,tol = 1e-4)
    # Test fixed prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = FALSE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A,BA,F,T,tol = 1e-4)
    # Test estimated prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = TRUE, estimate_prior_method='optim')
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy, estimate_prior_variance_method='optim')
    BA = report_susie_model(d.copy, B)
    # FIXME: have to use bigger tolerance level ...
    expect_susieR_equal(A,BA,TRUE,TRUE,tol = 1e-4)
    # Test estimated prior using EM algorithm and estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = TRUE, estimate_prior_method='EM')
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy, estimate_prior_variance_method='EM')
    BA = report_susie_model(d.copy, B)
    # FIXME: have to use bigger tolerance level ...
    skip('susieR estimate prior variance is inconsistent with mvsusieR')
    expect_susieR_equal(A, BA, T, T, 1E-6)
}))

test_that("mash regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    # Test fixed prior fixed residual
    residual_var = as.numeric(var(y))
    sigma = sd(y)
    n = length(y)
    sigma = sigma / sqrt(n)
    prior_var = as.numeric(scale_covariance(V[1,1]/residual_var, sigma))
    A = mvsusie(X,y,L=L,prior_variance=prior_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    m_init = MashInitializer$new(list(V), 1, 1, 0)
    B = mvsusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, 
               estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    expect_susie_equal(A,B,F,F)
}))

test_that("multivariate regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    # Test fixed prior fixed residual
    residual_var = as.numeric(var(y))
    sigma = sd(y)
    n = length(y)
    sigma = sigma / sqrt(n)
    prior_var = as.numeric(scale_covariance(V[1,1]/residual_var, sigma))
    A = mvsusie(X,y,L=L,prior_variance=prior_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE, precompute_covariances=TRUE)
    B = mvsusie(X,y,L=L,prior_variance=V,compute_objective=FALSE, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    m_init = MashInitializer$new(list(V), 1, 1, 0)
    A = expect_warning(mvsusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, 
                              estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, 
                              precompute_covariances=FALSE))
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
    # compare estimate prior variance "simple" method
    A = mvsusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'simple', precompute_covariances=TRUE)
    B = mvsusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'simple')
    expect_susie_equal(A,B,F,F)
    # compare estimate prior variance "EM" method
    A = mvsusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
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
    expect_equal(m_init$prior_variance$pi[1], null_weight)
}))

test_that("mvsusieR is identical to susieR (RSS)", with(simulate_univariate(summary = T), {
  # Test fixed prior fixed residual
  V = V/var(y)
  A = susieR::susie_rss(z, R, n = n, L = L, scaled_prior_variance = V, estimate_prior_variance = F)
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
  d.copy = d$clone(T)
  B$fit(d.copy)
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A,BA,F,F,rss = TRUE,tol = 1e-4)
  # Test estimated prior fixed residual
  A = susieR::susie_rss(z, R, n=n, L = L, scaled_prior_variance = V, estimate_prior_variance = TRUE)
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
  d.copy = d$clone(T)
  B$fit(d.copy,estimate_prior_variance_method='optim')
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A,BA,TRUE,FALSE,rss = TRUE,tol = 1e-4)
}))

test_that("mash regression in SuSiE is identical to univariate case (RSS)", with(simulate_multivariate(r=1), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  n=nrow(X)
  residual_var = as.numeric(var(y))
  prior_var = V[1,1]
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=prior_var,
             estimate_prior_variance=FALSE,
             compute_objective=FALSE)
  m_init = MashInitializer$new(list(V*n), 1, 1, 0)
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE,estimate_prior_variance=FALSE,
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
  m_init = MashInitializer$new(list(V), 1, 1, 0)
  A = expect_warning(mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE, precompute_covariances=FALSE, estimate_prior_method="simple"))
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE, precompute_covariances=TRUE, estimate_prior_method="simple")
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
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init)
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=V)
  expect_susie_equal(A,B,F,F)
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init, compute_objective=T)
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=V, compute_objective=T)
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
  A = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,compute_objective=FALSE)
  B = mvsusie_rss(z,R,N=n,L=L,prior_variance=m_init,s_init=A,compute_objective=FALSE)
  # let's just test of null is null ...
  null_weight = 0.2
  m_init = create_mixture_prior(R = ncol(y),null_weight = null_weight, max_mixture_len=-1)
  expect_equal(m_init$prior_variance$pi[1], null_weight)
}))
