# SUm of SIngle Effect (SuSiE) regression
SuSiE <- R6Class("SuSiE",
  public = list(
    sigma2 = NULL, # residual variance
    initialize = function(m, L,
        scaled_prior_variance, residual_variance,
        estimate_prior_variance, estimate_residual_variance, 
        max_iter=100,tol=1e-3,track_pip=FALSE,track_lbf=FALSE) 
    {
        private$elbo = vector()
        private$niter = max_iter
        private$tol = tol
        if (track_pip) private$pip_history = list()
        if (track_lbf) private$lbf_history = list()
        # initialize single effect regression models
        private$L = L
        private$kl = rep(NA,L)
        private$SER_models = lapply(1:private$L, function(j) m$clone(deep=T))
        self$set_prior_b(scaled_prior_variance * residual_variance)
    },
    fit = function(d) {
        for(i in 1:private$niter) {
            private$save_history()
            for (l in 1:L) {
                fitted_l = private$SER_models[l].predict(d)
                d.remove_from_fitted(fitted_l)
                d.comp_residual()
                private$SER_models[l]$fit(d)
                private$SER_models[l]$comp_posterior_expected_loglik(d, self$sigma2)
                private$kl[l] = private$SER_models[l]$eloglik - private$SER_models[l]$lbf
                fitted_l = private$SER_models[l].predict(d)
                d.add_back_fitted(fitted_l)
            }
            private$elbo[i] = private$compute_objective(d)
            private$estimate_residual_variance(d.get_n_sample())
            if (private$is_converged()) {
                private$save_history()
                private$niter = i
                break
            }
        }
    },
    predict = function(x) {

    },
    coef = function(x) {

    },
    # some get private numbers functions
    get_prior_b = function() { 
        # get prior effect size, because it might be updated during iterations
        private$exit()
    },
    get_objective = function() {},
    get_pip = function() {}, # posterior inclusion probability, L by J matrix
    get_pip_history = function() return(private$pip_history)
    get_lbf = function() {},
    get_lbf_history = function() return(private$lbf_history)
  ),
  private = list(
    L = NULL,
    SER_models = NULL, # Single effect regression models
    elbo = NULL, # Evidence lower bound
    kl = NULL, # KL divergence
    null_index = NULL, # index of null effect intentially added
    niter = NULL,
    pip_history = NULL, # keep track of pip
    lbf_history = NULL, # keep track of lbf
    tol = NULL, # tolerance level for convergence
    temp_erss = NULL,
    is_converged = function() {
        n = length(private$elbo)
        if (n<=1) return (FALSE)
        else return ((private$elbo[n]-private$elbo[n-1]) < private$tol)
    },
    compute_objective = function(d) {
        # update ERSS here
    },
    estimate_residual_variance = function(n) { 
        private$sigma2 = private$temp_erss / (n - 1)
    },
    save_history = function() {
        private$exit()
    },
    exit = function() {
        warning("Not yet implemented")
    }
  )
)