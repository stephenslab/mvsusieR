#' @title Single effect model object
#  It is meant to be dynamically inheriated from any regression model
#' @importFrom R6 R6Class
#' @keywords internal
SingleEffectModel <- function(base)
    R6Class("SingleEffectModel",
    inherit = base,
    public = list(
        initialize = function(J, residual_variance, prior_variance) {
            super$initialize(J, residual_variance, prior_variance)
            private$.pip = rep(0, J)
        },
        fit = function(d, prior_weights=NULL, estimate_prior_variance_method=NULL) {
            if (is.null(prior_weights)) prior_weights = rep(1/private$J, private$J)
            super$fit(d, use_residual = TRUE, prior_weights = prior_weights, estimate_prior_variance_method=estimate_prior_variance_method)
            ws = compute_weighted_sum(private$.lbf, prior_weights, log = TRUE)
            private$.pip = ws$weights
            private$lbf_single_effect = ws$log_sum
        },
        predict = function(d) {
            d$compute_Xb(self$posterior_b1)
        },
        compute_kl = function(d) {
            # compute KL divergence
            pp_eloglik = private$compute_expected_loglik_partial(d)
            private$.kl = pp_eloglik - private$lbf_single_effect
        }
    ),
    private = list(
        .pip = NULL, # posterior inclusion probability, alpha
        lbf_single_effect = NULL, # logBF for SER model: sum of logBF of all J effects
        .kl = NULL,
        .vbxxb = NULL,
        # This is the expected loglik minus the loglik_null, that is, N(R|B,V) - N(R|0,V)
        compute_expected_loglik_partial = function(d) {
            if (inherits(d, c("DenseData", "SSData", "SparseData"))) {
                if (is.null(private$.residual_variance_inv)) {
                    return(- (0.5/private$.residual_variance) * (- 2*sum(self$posterior_b1*d$XtR) + sum(d$X2_sum*as.vector(self$posterior_b2))))
                } else {
                    # FIXME: Currently this computation does not work for case with missing data
                    E1 = tr(private$.residual_variance_inv %*% t(self$posterior_b1) %*% d$XtR)
                    # posterior variance covariance matrix, weighted by PIP
                    if (length(dim(private$.posterior_b2)) == 3) {
                        S = lapply(1:nrow(private$.posterior_b1), function(j) private$.pip[j] * private$.posterior_b2[,,j] - tcrossprod(private$.pip[j] * private$.posterior_b1[j,]))
                    } else {
                        S = lapply(1:nrow(private$.posterior_b1), function(j) private$.pip[j] * matrix(private$.posterior_b2[j,]) - tcrossprod(private$.pip[j] * private$.posterior_b1[j,]))
                    }
                    private$.vbxxb = sum(d$X2_sum * sapply(1:length(S), function(j) tr(private$.residual_variance_inv %*% S[[j]])))
                    private$.vbxxb = sum(d$X2_sum * sapply(1:length(S), function(j) t(self$posterior_b1[j,]) %*% private$.residual_variance_inv %*% self$posterior_b1[j,])) + private$.vbxxb
                    return(E1 - private$.vbxxb / 2)
                }
            } else {
                stop("compute_expected_loglik_partial not implemented for given data type")
            }
        }
    ),
    active = list(
        # allow for initialization of coefficients
        mu = function(v) {
            if (missing(v)) private$.posterior_b1
            else private$.posterior_b1 = v
        },
        # posterior first moment, alpha * posterior_b1_reg
        posterior_b1 = function() private$.pip * private$.posterior_b1,
        # posterior first moment, alpha * posterior_b2_reg
        posterior_b2 = function() {
            if (length(dim(private$.posterior_b2)) == 3) b2 = t(apply(private$.posterior_b2, 3, diag))
            else b2 = private$.posterior_b2
            return(private$.pip * b2)
        },
        pip = function() private$.pip,
        lbf = function() private$lbf_single_effect,
        kl = function() private$.kl,
        vbxxb = function() private$.vbxxb
    )
  )