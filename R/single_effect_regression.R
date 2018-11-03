# Single effect regression
SingleEffectRegression <- R6Class("SingleEffectRegression",
  public = list(
    pip = NULL, # posterior inclusion probability, alpha
    posterior_b1 = NULL, # posterior first moment, alpha * posterior_b1_reg
    posterior_b2 = NULL, # posterior second moment, alpha * posterior_b2_reg
    lbf = NULL, # log Bayes factors from base regression computation
    initialize = function(m,J,prior=NULL) {
        private$J = J
        private$models = lapply(1:private$J, function(j) m$clone(deep=T))
        self$set_prior(prior)
    },
    set_prior = function(prior=NULL) {
        # set prior
        if (is.null(prior)) private$prior = rep(1/private$J,private$J)
        else private$prior = prior
    },
    fit = function(d) {
        private$fit_jth_model(d)
        private$comp_bf()
        private$comp_pip()
        private$comp_posterior_b()
    },
    # some get private numbers functions
    get_prior = function() { private$prior }
  ),
  private = list(
    models = NULL, # Bayesian regression for the j-th element
    J = NULL,
    prior = NULL, # prior on gamma
    fit_jth_model = function(d) {
        invisible(lapply(1:private$J, function(j) private$models[j]$fit(d,j)))
    }
    comp_lbf = function() {
        # compute Bayes factor
        self$lbf = sapply(1:private$J, function(j) private$models[j]$lbf)
    },
    comp_pip = function() { 
        # compute posterior inclusion probability
        self$pip = safe_comp_weight(self$lbfs, private$prior, log = TRUE)
    },
    comp_posterior_b = function() {
        # compute posterior on b
        self$posterior_b1 = lapply(1:private$J, function(j) self$pip[j] * private$models[j]$posterior_b1)
        self$posterior_b2 = lapply(1:private$J, function(j) self$pip[j] * private$models[j]$posterior_b2)
    },
    exit = function() {
        warning("Not yet implemented")
    }
  )
)