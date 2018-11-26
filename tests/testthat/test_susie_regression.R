context("Test SuSiE regression")

test_that("mmbr is identical to susieR", with(simulate_univariate(), {
    # Test fixed prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = FALSE)
    SER = SingleEffectRegression(BayesianMultipleRegression)$new(d$n_effect, 1, V, estimate_prior_variance = FALSE, prior_weights = NULL)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B) 
    expect_susieR_equal(A, BA, F, F)
    # Test estimated prior fixed residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = FALSE, estimate_prior_variance = TRUE)
    SER = SingleEffectRegression(BayesianMultipleRegression)$new(d$n_effect, 1, V, estimate_prior_variance = TRUE, prior_weights = NULL)
    B = SuSiE$new(SER, L, estimate_residual_variance = FALSE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B) 
    expect_susieR_equal(A, BA, T, F)
    # Test fixed prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = FALSE)
    SER = SingleEffectRegression(BayesianMultipleRegression)$new(d$n_effect, 1, V, estimate_prior_variance = FALSE, prior_weights = NULL)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B) 
    expect_susieR_equal(A, BA, F, T)
    # Test estimated prior estimated residual
    A = susieR::susie(X, y, L = L, scaled_prior_variance = V/var(y), residual_variance = 1, prior_weights = NULL, estimate_residual_variance = TRUE, estimate_prior_variance = TRUE)
    SER = SingleEffectRegression(BayesianMultipleRegression)$new(d$n_effect, 1, V, estimate_prior_variance = TRUE, prior_weights = NULL)
    B = SuSiE$new(SER, L, estimate_residual_variance = TRUE)
    d.copy = d$clone(T)
    B$fit(d.copy)
    BA = report_susie_model(d.copy, B) 
    expect_susieR_equal(A, BA, T, T)
}))

# test_that("mash regression in SuSiE is identical to univariate case", with(simulate_multivariate(r=1), {
#     prior_var = V[1,1]
#     residual_var = as.numeric(var(y))
#     A = susie(X,y,L=L,V=prior_var,compute_objective=FALSE)
#     residual_cov = cov(y)
#     m_init = MashInitializer$new(list(V), 1, 1, 0)
#     B = susie(X,y,L=L,V=m_init,compute_objective=FALSE)
#     expect_susie_equal(A,B,F,F)
# }))