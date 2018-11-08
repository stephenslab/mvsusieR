#' @title SUm of SIngle Effect (SuSiE) regression object
#' @importFrom R6 R6Class
#' @keywords internal
SuSiE <- R6Class("SuSiE",
  public = list(
    initialize = function(SER,L,estimate_residual_variance,
    max_iter=100,tol=1e-3,track_pip=FALSE,track_lbf=FALSE) 
    {
        # initialize single effect regression models
        private$L = L
        private$to_estimate_residual_variance = estimate_residual_variance
        private$SER = lapply(1:private$L, function(l) SER$clone(deep=T))
        private$sigma2 = SER$residual_variance
        private$elbo = vector()
        private$niter = max_iter
        private$tol = tol
        if (track_pip) private$pip_history = list()
        if (track_lbf) private$lbf_history = list()
    },
    fit = function(d) {
        for (i in 1:private$niter) {
            private$save_history()
            for (l in 1:private$L) {
                d$remove_from_fitted(private$SER[[l]]$predict(d))
                d$compute_residual()
                private$SER[[l]]$residual_variance = private$sigma2
                private$SER[[l]]$fit(d)
                private$SER[[l]]$compute_kl(d)
                d$add_back_fitted(private$SER[[l]]$predict(d))
            }
            # assign these variables for performance consideration
            # because these active bindings involve some computation
            b1 = self$posterior_b1
            b2 = self$posterior_b2
            private$estimate_residual_variance(d,b1,b2)
            private$elbo[i] = private$compute_objective(d,b1,b2)
            if (private$is_converged()) {
                private$save_history()
                private$niter = i
                break
            }
        }
    },
    predict = function(x) {},
    coef = function(d) {
        d$rescale_coef(Reduce(`+`, self$posterior_b1))
    },
    get_objective = function(dump = FALSE) {
        if (!all(diff(private$elbo) >= 0)) {
            warning('Objective is not non-decreasing')
            dump = TRUE
        }
        if (dump) return(private$elbo)
        else return(private$elbo[private$niter])
    },
    get_niter = function() private$niter,
    get_pip_history = function() private$pip_history,
    get_lbf_history = function() private$lbf_history
  ),
  private = list(
    to_estimate_residual_variance = NULL,
    L = NULL,
    SER = NULL, # Single effect regression models
    elbo = NULL, # Evidence lower bound
    null_index = NULL, # index of null effect intentially added
    niter = NULL,
    pip_history = NULL, # keep track of pip
    lbf_history = NULL, # keep track of lbf
    tol = NULL, # tolerance level for convergence
    sigma2 = NULL, # residual variance
    essr = NULL,
    is_converged = function() {
        n = length(private$elbo)
        if (n<=1) return (FALSE)
        else return ((private$elbo[n]-private$elbo[n-1]) < private$tol)
    },
    compute_objective = function(d,b1,b2) {
        if (is.null(private$essr)) {
            essr = compute_expected_sum_squared_residuals(d,b1,b2) 
        } else {
            essr = private$essr
        }
        expected_loglik = compute_expected_loglik(d$n_sample, private$sigma2, essr)
        return(expected_loglik - sum(self$kl))
    },
    estimate_residual_variance = function(d,b1,b2) { 
        if (private$to_estimate_residual_variance) {
            private$essr = compute_expected_sum_squared_residuals(d,b1,b2)
            # FIXME: should we bother with it being N - 1 (or N - 2)?
            private$sigma2 = private$essr / d$n_sample
        }
    },
    save_history = function() {
        warning('save_histroy not yet implemented')
    },
    denied = function(v) stop(paste0('$', v, ' is read-only'), call. = FALSE) 
  ),
  active = list(
    prior_variance = function(v) { 
        # get prior effect size, because it might be updated during iterations
        if (missing(v)) sapply(1:private$L, function(l) private$SER[[l]]$prior_variance)
        else private$denied('prior_variance')
    },
    residual_variance = function(v) {
        if (missing(v)) private$sigma2
        else private$denied('residual_variance')
    },
    kl = function(v) {
        if (missing(v)) sapply(1:private$L, function(l) private$SER[[l]]$kl)
        else private$denied('kl')
    },
    pip = function(v) {
        if (missing(v)) do.call(cbind, lapply(1:private$L, function(l) private$SER[[l]]$pip)) # posterior inclusion probability, J by L matrix
        else private$denied('pip')
    },
    lbf = function(v) {
        if (missing(v)) sapply(1:private$L, function(l) private$SER[[l]]$lbf_single_effect)
        else priviate$denied('lbf')
    },
    posterior_b1 = function(v) {
        if (missing(v)) lapply(1:private$L, function(l) private$SER[[l]]$posterior_b1)
        else priviate$denied('posterior_b1')
    },
    posterior_b2 = function(v) {
        if (missing(v)) lapply(1:private$L, function(l) private$SER[[l]]$posterior_b2)
        else priviate$denied('posterior_b2')
    }
  )
)

# expected squared residuals
compute_expected_sum_squared_residuals = function(d, Eb1, Eb2) {
    if (inherits(d, c("DenseData","SparseData", "SSData"))) {
        Eb1 = t(do.call(cbind, Eb1))
        Eb2 = t(do.call(cbind, Eb2))
    }
    if (inherits(d, c("DenseData","SparseData"))) {
        Xr = d$compute_MXt(Eb1)
        Xrsum = colSums(Xr)
        return(sum((d$Y-Xrsum)^2) - sum(Xr^2) + sum(d$d*t(Eb2)))
    } else {
        XB2 = sum((Eb1%*%d$XtX) * Eb1)
        betabar = colSums(Eb1)
        return(d$YtY - 2*sum(betabar * d$XtY) + sum(betabar * (d$XtX %*% betabar)) -
           XB2 + sum(d$d*t(Eb2)))
    }
}

# expected loglikelihood for a susie fit
compute_expected_loglik = function(n, residual_variance, essr) {
    -(n/2) * log(2*pi* residual_variance) - (1/(2*residual_variance)) * essr
}