context("Test SuSiE regression")

test_that("mmbr is identical to susieR when prior is a scalar", with(simulate_univariate(), {
    # Test fixed prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = FALSE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A, BA, F, F)
    # Test estimated prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = TRUE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(T)
    B$fit(d.copy,estimate_prior_variance_method='optim')
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A, BA, T, F, 1E-5)
    # Test fixed prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = FALSE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A, BA, F, T)
    # Test estimated prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = TRUE, estimate_prior_method='optim')
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy, estimate_prior_variance_method='optim')
    BA = report_susie_model(d.copy, B)
    # FIXME: have to use bigger tolerance level ...
    expect_susieR_equal(A, BA, T, T, 1E-6)
    # Test estimated prior using EM algorithm and estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = TRUE, estimate_prior_method='EM')
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy, estimate_prior_variance_method='EM')
    BA = report_susie_model(d.copy, B)
    # FIXME: have to use bigger tolerance level ...
    skip('susieR estimate prior variance is inconsistent with mmbr')
    expect_susieR_equal(A, BA, T, T, 1E-6)
}))

test_that("mash regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    # Test fixed prior fixed residual
    residual_var = as.numeric(var(y))
    scaled_prior_var = V[1,1] / residual_var
    A = msusie(X,y,L=L,prior_variance=scaled_prior_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE)
    m_init = MashInitializer$new(list(V), 1, 1, 0)
    B = msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    expect_susie_equal(A,B,F,F)
}))

test_that("multivariate regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    # Test fixed prior fixed residual
    residual_var = as.numeric(var(y))
    scaled_prior_var = V[1,1] / residual_var
    A = msusie(X,y,L=L,prior_variance=scaled_prior_var,
                estimate_residual_variance=FALSE,
                estimate_prior_variance=FALSE,
                compute_objective=FALSE, precompute_covariances=TRUE)
    B = msusie(X,y,L=L,prior_variance=V,compute_objective=FALSE, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    m_init = MashInitializer$new(list(V), 1, 1, 0)
    A = expect_warning(msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, precompute_covariances=FALSE))
    B = msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, precompute_covariances=TRUE, estimate_prior_variance=FALSE, estimate_residual_variance=FALSE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with BMR using one component prior matrix", with(simulate_multivariate(r=3), {
    m_init = create_mash_prior(mixture_prior = list(matrices=list(V)))
    # don't compare ELBO
    A = msusie(X,y,L=L,prior_variance=m_init, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    B = msusie(X,y,L=L,prior_variance=V, estimate_residual_variance=FALSE, estimate_prior_variance=FALSE)
    expect_susie_equal(A,B,F,F)
    # compare ELBO
    A = msusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=FALSE, precompute_covariances=TRUE)
    B = msusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=FALSE)
    expect_susie_equal(A,B,F,F)
    # compare estimate prior variance "simple" method
    A = msusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'simple', precompute_covariances=TRUE)
    B = msusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'simple')
    expect_susie_equal(A,B,F,F)
    # compare estimate prior variance "EM" method
    A = msusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM', precompute_covariances=TRUE)
    B = msusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, estimate_prior_variance=TRUE, estimate_prior_method = 'EM')
    expect_susie_equal(A,B,F,F)
}))

test_that("customized initialization interface", with(simulate_multivariate(r=3), {
    # not sure what to test here ...
    m_init = create_mash_prior(mixture_prior = list(matrices = list(V), weights = 1), null_weight=0)
    A = msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE, estimate_prior_variance=FALSE)
    B = msusie(X,y,L=L,prior_variance=m_init,s_init=A,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE, estimate_prior_variance=FALSE)
    # let's just test of null is null ...
    null_weight = 0.2
    m_init = create_mash_prior(sample_data = list(X=X,Y=y,center=T,scale=T,residual_variance=cov(y)),
                                null_weight = null_weight, max_mixture_len=-1)
    expect_equal(m_init$prior_variance$pi[1], null_weight)
}))

test_that("mmbr is identical to susieR (RSS)", with(simulate_univariate(summary = T), {
  # Test fixed prior fixed residual
  A = susieR::susie_rss(z, R, L = L, prior_variance = V, residual_variance = 1, prior_weights = NULL,
                        estimate_residual_variance = FALSE, estimate_prior_variance = F, check_z = F)
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
  d.copy = d$clone(T)
  B$fit(d.copy)
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A, BA, F, F)
  # Test estimated prior fixed residual
  A = susieR::susie_rss(z, R, L = L, prior_variance = V, residual_variance = 1, prior_weights = NULL,
                        estimate_residual_variance = FALSE, estimate_prior_variance = TRUE, check_z = F)
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
  d.copy = d$clone(T)
  B$fit(d.copy,estimate_prior_variance_method='optim')
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A, BA, T, F, 1E-5)
  # Test fixed prior estimated residual
  A = susieR::susie_rss(z, R, L = L, prior_variance = V, residual_variance = 1, prior_weights = NULL,
                        estimate_residual_variance = TRUE, estimate_prior_variance = FALSE, check_z = F)
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
  d.copy = d$clone(T)
  B$fit(d.copy)
  BA = report_susie_model(d.copy, B)
  expect_susieR_equal(A, BA, F, T)
  # Test estimated prior estimated residual
  A = susieR::susie_rss(z, R, L = L, prior_variance = V, residual_variance = 1, prior_weights = NULL,
                        estimate_residual_variance = TRUE, estimate_prior_variance = TRUE, check_z = F)
  SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, V)
  B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
  d.copy = d$clone(T)
  B$fit(d.copy, estimate_prior_variance_method='optim')
  BA = report_susie_model(d.copy, B)
  # FIXME: have to use bigger tolerance level ...
  expect_susieR_equal(A, BA, T, T, 1E-6)
}))

test_that("mash regression in SuSiE is identical to univariate case (RSS)", with(simulate_multivariate(r=1), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  residual_var = as.numeric(var(y))
  prior_var = V[1,1]
  A = msusie_rss(z,R,L=L,prior_variance=prior_var,
             estimate_residual_variance=FALSE,
             estimate_prior_variance=FALSE,
             compute_objective=FALSE)
  m_init = MashInitializer$new(list(V), 1, 1, 0)
  B = msusie_rss(z,R,L=L,prior_variance=m_init,compute_objective=FALSE,estimate_prior_variance=FALSE,
                 estimate_residual_variance=FALSE, precompute_covariances=TRUE)
  expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed (RSS)", with(simulate_multivariate(r=3), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  m_init = MashInitializer$new(list(V), 1, 1, 0)
  A = expect_warning(msusie_rss(z,R,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=FALSE))
  B = msusie_rss(z,R,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE)
  expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with BMR using one component prior matrix (RSS)", with(simulate_multivariate(r=3), {
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  m_init = create_mash_prior(mixture_prior = list(matrices=list(V)))
  A = msusie_rss(z,R,L=L,prior_variance=m_init, estimate_residual_variance=FALSE, precompute_covariances=TRUE)
  B = msusie_rss(z,R,L=L,prior_variance=V, estimate_residual_variance=FALSE, precompute_covariances=TRUE)
  expect_susie_equal(A,B,F,F)
  A = msusie_rss(z,R,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F, precompute_covariances=TRUE)
  B = msusie_rss(z,R,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F, precompute_covariances=TRUE)
  expect_susie_equal(A,B,F,F)
}))

test_that("customized initialization interface (RSS)", with(simulate_multivariate(r=3), {
  # not sure what to test here ...
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  m_init = create_mash_prior(mixture_prior = list(matrices = list(V), weights = 1), null_weight=0)
  A = msusie_rss(z,R,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE)
  B = msusie_rss(z,R,L=L,prior_variance=m_init,s_init=A,compute_objective=FALSE, estimate_residual_variance=FALSE, precompute_covariances=TRUE)
  # let's just test of null is null ...
  null_weight = 0.2
  m_init = create_mash_prior(sample_data = list(X=X,Y=y,center=T,scale=T,residual_variance=cov(y)),
                             null_weight = null_weight, max_mixture_len=-1)
  expect_equal(m_init$prior_variance$pi[1], null_weight)
}))
