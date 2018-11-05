# Single effect regression
SingleEffectRegression <- R6Class("SingleEffectRegression",
  public = list(
    pip = NULL, # posterior inclusion probability, alpha
    lbf = NULL, # logBF for SER model: sum of logBF of all J effects
    loglik = NULL,
    kl = NULL,
    initialize = function(m,J,prior_gamma=NULL) {
        private$J = J
        private$BMR = m$clone(deep=TRUE)
        if (is.null(prior_gamma)) private$prior_gamma = rep(1/private$J,private$J)
        else private$prior_gamma = prior_gamma
    },
    fit = function(d, residual_variance) {
        private$BMR$fit(d, use_residual = TRUE, prior_weights = private$prior_gamma)
        ws = safe_comp_weight(private$BMR$lbf, private$prior_gamma, log = TRUE)
        self$pip = ws$alpha
        self$lbf = ws$log_total
        if (!is.null(d$Y)) self$loglik = self$lbf + sum(dnorm(d$Y,0,sqrt(residual_variance),log=TRUE))
        # keep track of prior for this SER
        # because it might have changed due to re-estimates in BMR
        private$prior_b = BMR$get_prior()
    },
    predict = function(d) {
        d$compute_Xb(self$get_posterior_b1())
    }
    get_posterior_b1 = function() {
        # posterior first moment, alpha * posterior_b1_reg
        self$pip * private$BMR$posterior_b1
    },
    get_posterior_b2 = function() {
        # posterior second moment, alpha * posterior_b2_reg
        self$pip * private$BMR$posterior_b2
    },
    comp_kl = function(d, residual_variance) {
        # compute KL divergence
        pp_eloglik = comp_expected_loglik_partial(d, residual_variance, 
                                                  self$get_posterior_b1(),
                                                  self$get_posterior_b2())
        self$kl = pp_eloglik - self$lbf
    }
  ),
  private = list(
    BMR = NULL, # Bayesian multiple regression model
    prior_gamma = NULL, # prior on gamma
    prior_b = NULL,
    J = NULL,
    exit = function() {
        warning("Not yet implemented")
    }
  )
)

comp_expected_loglik_partial = function(d, s2, Eb1, Eb2) {
    if (inherits(d, c("DenseData", "SSData", "SparseData"))) {
        return(- (0.5/s2) * (- 2*sum(Eb1*d$get_XtR()) + sum(d$d*as.vector(Eb2))))
    } else {
        stop("comp_expected_loglik_partial not implemented for given data type")
    }
}