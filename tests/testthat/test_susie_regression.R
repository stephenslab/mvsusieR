context("Test SuSiE regression")

test_that("mmbr is identical to susieR", with(simulate_univariate(), {
    # Test fixed prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = FALSE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, 1, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A, BA, F, F)
    # Test estimated prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = TRUE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, 1, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(T)
    B$fit(d.copy,estimate_prior_variance_method='optim')
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A, BA, T, F, 1E-5)
    # Test fixed prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = FALSE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, 1, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B)
    expect_susieR_equal(A, BA, F, T)
    # Test estimated prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = TRUE)
    SER = SingleEffectModel(BayesianSimpleRegression)$new(d$n_effect, 1, V)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy, estimate_prior_variance_method='optim')
    BA = report_susie_model(d.copy, B)
    # FIXME: have to use bigger tolerance level ...
    expect_susieR_equal(A, BA, T, T, 1E-6)
}))

test_that("mash regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
    residual_var = as.numeric(var(y))
    scaled_prior_var = V[1,1] / residual_var
    A = msusie(X,y,L=L,prior_variance=scaled_prior_var,
                estimate_residual_variance=FALSE,
                compute_objective=FALSE)
    m_init = MashInitializer$new(list(V), 1, 1, 0)
    B = msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with when various covariance quantities are precomputed", with(simulate_multivariate(r=3), {
    m_init = MashInitializer$new(list(V), 1, 1, 0)
    A = msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, estimate_residual_variance=FALSE)
    B = msusie(X,y,L=L,prior_variance=m_init,compute_objective=FALSE, precompute_covariances=TRUE, estimate_residual_variance=FALSE)
    expect_susie_equal(A,B,F,F)
}))

test_that("mash regression in SuSiE agrees with BMR using one component prior matrix", with(simulate_multivariate(r=3), {
    m_init = create_mash_prior(mixture_prior = list(matrices=list(V)))
    A = msusie(X,y,L=L,prior_variance=m_init, estimate_residual_variance=FALSE)
    B = msusie(X,y,L=L,prior_variance=V, estimate_residual_variance=FALSE)
    expect_susie_equal(A,B,F,F)
    A = msusie(X,y,L=L,prior_variance=m_init, compute_objective=T, estimate_residual_variance=F)
    B = msusie(X,y,L=L,prior_variance=V, compute_objective=T, estimate_residual_variance=F)
    expect_susie_equal(A,B,F,F)
}))

test_that("customized initialization interface", with(simulate_multivariate(r=3), {
    # not sure what to test here ...
    m_init = create_mash_prior(mixture_prior = list(matrices = list(V), weights = 1), null_weight=0)
    A = msusie(X,y,L=L,prior_variance=m_init,s_init=list(coef_index=c(2,3,4),coef_value=matrix(1,3,3)),compute_objective=FALSE, estimate_residual_variance=FALSE)
    # let's just test of null is null ...
    null_weight = 0.2
    m_init = create_mash_prior(sample_data = list(X=X,Y=y,center=T,scale=T,residual_variance=cov(y)),
                                null_weight = null_weight, max_mixture_len=-1)
    expect_equal(m_init$prior_variance$pi[1], null_weight)
}))