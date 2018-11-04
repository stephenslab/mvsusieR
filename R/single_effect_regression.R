# Single effect regression
SingleEffectRegression <- R6Class("SingleEffectRegression",
  public = list(
    pip = NULL, # posterior inclusion probability, alpha
    lbf = NULL, # log Bayes Factor for the single effect regression model
    initialize = function(m,J,prior_gamma=NULL) {
        private$J = J
        private$BR_models = lapply(1:private$J, function(j) m$clone(deep=T))
        self$set_prior_gamma(prior_gamma)
    },
    set_prior_gamma = function(prior_gamma=NULL) {
        # set prior inclusion probability
        if (is.null(prior_gamma)) private$prior_gamma = rep(1/private$J,private$J)
        else private$prior_gamma = prior_gamma
    },
    set_prior_b = function(prior=NULL) {
        # set prior
        if (is.null(prior)) private$prior_b = rep(1/private$J,private$J)
        else private$prior_b = prior
    },
    fit = function(d) {
        private$fit_jth_model(d)
        private$comp_pip()
    },
    # some get private numbers functions
    get_pi = function() { private$prior_gamma },
    get_BR_lbf = function() {
        return(sapply(1:private$J, function(j) private$BR_models[j]$lbf))
    },
    get_posterior_b1 = function() {
        # posterior first moment, alpha * posterior_b1_reg
        self$posterior_b1 = lapply(1:private$J, function(j) self$pip[j] * private$BR_models[j]$posterior_b1)
    },
    get_posterior_b2 = function() {
        # posterior second moment, alpha * posterior_b2_reg
        self$posterior_b2 = lapply(1:private$J, function(j) self$pip[j] * private$BR_models[j]$posterior_b2)
    }
  ),
  private = list(
    BR_models = NULL, # Bayesian regression for the j-th element
    prior_gamma = NULL, # prior on gamma
    prior_b = NULL, # prior on effect size
    J = NULL,
    fit_jth_model = function(d) {
        invisible(lapply(1:private$J, function(j) private$BR_models[j]$fit(d,j)))
    }
    comp_pip = function() { 
        # compute posterior inclusion probability
        self$pip = safe_comp_weight(self$get_lbf, private$prior_gamma, log = TRUE)
    },
    exit = function() {
        warning("Not yet implemented")
    }
  )
)