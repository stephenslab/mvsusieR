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
        private$kl = vector()
        private$niter = max_iter
        private$tol = tol
        if (track_pip) private$pip_history = list()
        if (track_lbf) private$lbf_history = list()
        # initialize single effect regression models
        private$L = L
        private$workers = lapply(1:private$J, function(j) m$clone(deep=T))
        self$set_prior_b(scaled_prior_variance * residual_variance)
    },
    set_prior = function(prior=NULL) {
        # set prior
        if (is.null(prior)) private$prior = rep(1/private$J,private$J)
        else private$prior = prior
    },
    fit = function(d) {
        for(i in 1:private$niter) {
            private$save_history()
            private$update_each_effect(d)
            private$compute_objective(d)
            private$estimate_residual_variance(d)
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
    workers = NULL, # Single effect regression models
    elbo = NULL, # Evidence lower bound
    kl = NULL, # KL divergence
    null_index = NULL, # index of null effect intentially added
    niter = NULL,
    pip_history = NULL, # keep track of pip
    lbf_history = NULL, # keep track of lbf
    tol = NULL, # tolerance level for convergence
    is_converged = function() {
        n = length(private$elbo)
        if (n<=1) return (FALSE)
        else return ((private$elbo[n]-private$elbo[n-1]) < private$tol)
    },
    update_each_effect = function(d) {
    },
    compute_objective = function() {
    },
    estimate_residual_variance = function() { 
    },
    save_history = function() {
        private$exit()
    },
    exit = function() {
        warning("Not yet implemented")
    }
  )
)