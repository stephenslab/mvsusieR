#' @title Single effect regression object
#  It is meant to be dynamically inheriated from any regression model
#' @importFrom R6 R6Class
#' @keywords internal
SingleEffectRegression <- function(base)
    R6Class("SingleEffectRegression",
    inherit = base,
    public = list(
        initialize = function(J, residual_variance, prior_variance, estimate_prior_variance=FALSE, prior_weights=NULL) {
            super$initialize(J, residual_variance, prior_variance, estimate_prior_variance)
            if (is.null(prior_weights)) private$prior_weights = rep(1/private$J,private$J)
            else private$prior_weights = prior_weights
            private$.pip = rep(0, J)
        },
        fit = function(d) {
            super$fit(d, use_residual = TRUE, prior_weights = private$prior_weights)
            ws = safe_compute_weight(private$.lbf, private$prior_weights, log = TRUE)
            private$.pip = ws$alpha
            private$.lbf_single_effect = ws$log_total
            self$compute_loglik_null(d)
            private$.mloglik_single_effect = private$.lbf_single_effect + private$.loglik_null
        },
        predict = function(d) {
            d$compute_Xb(self$posterior_b1)
        },
        compute_kl = function(d) {
            # compute KL divergence
            pp_eloglik = compute_expected_loglik_partial(d, private$.residual_variance,
                                                    self$posterior_b1,
                                                    self$posterior_b2)
            private$.kl = pp_eloglik - private$.lbf_single_effect
        }
    ),
    private = list(
        prior_weights = NULL, # prior on gamma
        .mloglik_single_effect = NULL,
        .pip = NULL, # posterior inclusion probability, alpha
        .lbf_single_effect = NULL, # logBF for SER model: sum of logBF of all J effects
        .kl = NULL
    ),
    active = list(
        # user accessible interface
        posterior_b1 = function(v) {
            # posterior first moment, alpha * posterior_b1_reg
            if (missing(v)) private$.pip * private$.posterior_b1
            else private$denied('posterior_b1')
        },
        posterior_b2 = function(v) {
            # posterior first moment, alpha * posterior_b2_reg
            if (missing(v)) {
                if (length(dim(private$.posterior_b2)) == 3)
                    b2 = t(apply(private$.posterior_b2, 3, diag))
                else
                    b2 = private$.posterior_b2
                private$.pip * b2
            } else { private$denied('posterior_b2') }
        },
        lfsr = function(v) {
            if (missing(v)) private$.lfsr
            else private$denied('lfsr')
        },
        kl = function(v) {
            if (missing(v)) private$.kl
            else private$.kl = v
        },
        pip = function(v) {
            if (missing(v)) private$.pip
            else private$.pip = v
        },
        lbf_single_effect = function(v) {
            if (missing(v)) private$.lbf_single_effect
            else private$.lbf_single_effect = v
        },
        mloglik_single_effect = function(v) {
            if (missing(v)) private$.mloglik_single_effect
            else private$denied('mloglik_single_effect')
        },
        mixture_posterior_weights = function(v) {
            if (missing(v)) private$.mixture_posterior_weights
            else private$denied('mixture_posterior_weights')
        }
    )
  )

compute_expected_loglik_partial = function(d, s2, Eb1, Eb2) {
    if (inherits(d, c("DenseData", "SSData", "SparseData"))) {
        return(- (0.5/s2) * (- 2*sum(Eb1*d$XtR) + sum(d$d*as.vector(Eb2))))
    } else {
        stop("compute_expected_loglik_partial not implemented for given data type")
    }
}