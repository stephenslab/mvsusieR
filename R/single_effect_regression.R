#' @title Single effect regression object
#  It is meant to be dynamically inheriated from any regression model
#' @importFrom R6 R6Class
#' @keywords internal
SingleEffectRegression <- function(base)
    R6Class("SingleEffectRegression",
    inherit = base,
    public = list(
        initialize = function(J, prior_variance, estimate_prior_variance, prior_weights=NULL) {
            super$initialize(prior_variance, estimate_prior_variance)
            private$J = J
            if (is.null(prior_weights)) private$prior_weights = rep(1/private$J,private$J)
            else private$prior_weights = prior_weights
        },
        fit = function(d) {
            super$fit(d, use_residual = TRUE, prior_weights = private$prior_weights)
            ws = safe_comp_weight(private$get_lbf(), private$prior_weights, log = TRUE)
            private$pip = ws$alpha
            private$lbf_single_effect = ws$log_total
            if (!is.null(d$Y)) private$mloglik_single_effect = private$lbf_single_effect + sum(dnorm(d$Y,0,sqrt(private$get_residual_variance()),log=TRUE))
        },
        predict = function(d) {
            d$compute_Xb(self$get_posterior_b1())
        },
        get_posterior_b1 = function() {
            # posterior first moment, alpha * posterior_b1_reg
            private$pip * super$get_posterior_b1()
        },
        get_posterior_b2 = function() {
            # posterior second moment, alpha * posterior_b2_reg
            private$pip * super$get_posterior_b2()
        },
        comp_kl = function(d) {
            # compute KL divergence
            pp_eloglik = comp_expected_loglik_partial(d, self$get_residual_variance(), 
                                                    self$get_posterior_b1(),
                                                    self$get_posterior_b2())
            private$kl = pp_eloglik - private$lbf_single_effect
        },
        get_kl = function() private$kl,
        get_pip = function() private$pip,
        get_lbf_single_effect = function() private$lbf_single_effect
    ),
    private = list(
        prior_weights = NULL, # prior on gamma
        J = NULL,
        mloglik_single_effect = NULL,
        pip = NULL, # posterior inclusion probability, alpha
        lbf_single_effect = NULL, # logBF for SER model: sum of logBF of all J effects
        kl = NULL,
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