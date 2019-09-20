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
        },
        predict = function(d) {
            d$compute_Xb(self$posterior_b1)
        },
        compute_kl = function(d) {
            # compute KL divergence
            pp_eloglik = private$compute_expected_loglik_partial(d)
            private$.kl = pp_eloglik - private$.lbf_single_effect
        }
    ),
    private = list(
        prior_weights = NULL, # prior on gamma
        .pip = NULL, # posterior inclusion probability, alpha
        .lbf_single_effect = NULL, # logBF for SER model: sum of logBF of all J effects
        .kl = NULL,
        .vbxxb = NULL,
        # This is the expected loglik minus the loglik_null, that is, N(R|B,V) - N(R|0,V)
        compute_expected_loglik_partial = function(d) {
            if (inherits(d, c("DenseData", "SSData", "SparseData"))) {
                if (is.null(private$.residual_variance_inv)) {
                    return(- (0.5/private$.residual_variance) * (- 2*sum(self$posterior_b1*d$XtR) + sum(d$d*as.vector(self$posterior_b2))))
                } else {
                    # FIXME: Currently this computation does not work for case with missing data
                    E1 = tr(private$.residual_variance_inv %*% t(self$posterior_b1) %*% d$XtR)
                    # posterior variance covariance matrix, weighted by PIP
                    if (length(dim(private$.posterior_b2)) == 3) {
                        S = lapply(1:nrow(private$.posterior_b1), function(j) private$.pip[j] * private$.posterior_b2[,,j] - tcrossprod(private$.pip[j] * private$.posterior_b1[j,]))
                    } else {
                        S = lapply(1:nrow(private$.posterior_b1), function(j) private$.pip[j] * matrix(private$.posterior_b2[j,]) - tcrossprod(private$.pip[j] * private$.posterior_b1[j,]))
                    }
                    private$.vbxxb = sum(d$d * sapply(1:length(S), function(j) tr(private$.residual_variance_inv %*% S[[j]])))
                    E2 = sum(d$d * sapply(1:length(S), function(j) t(self$posterior_b1[j,]) %*% private$.residual_variance_inv %*% self$posterior_b1[j,])) + private$.vbxxb
                    return(E1 - E2 / 2)
                }
            } else {
                stop("compute_expected_loglik_partial not implemented for given data type")
            }
        }
    ),
    active = list(
        # user accessible interface
        mu = function(v) {
            if (missing(v)) private$.posterior_b1
            else private$.posterior_b1 = v
        },
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
            # allow for initializing it
            else private$.pip = v
        },
        lbf_single_effect = function(v) {
            if (missing(v)) private$.lbf_single_effect
            else private$.lbf_single_effect = v
        },
        mixture_posterior_weights = function(v) {
            if (missing(v)) private$.mixture_posterior_weights
            else private$denied('mixture_posterior_weights')
        },
        vbxxb = function(v) {
            if (missing(v)) private$.vbxxb
            else private$denied('vbxxb')
        }
    )
  )