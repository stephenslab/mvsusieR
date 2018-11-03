# 
BR_model = BayesBaysianRegression()
SER_model = SingleEffectRegression(BR_model, data$J, prior_weights)
SuSiE_model = SuSiE(SER_model, L, scaled_prior_variance, residual_variance,estimate_prior_variance, 
                    estimate_residual_variance, max_iter,tol,track_pip,track_lbf)