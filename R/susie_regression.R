# SUm of SIngle Effect (SuSiE) regression
SuSiE <- R6Class("SuSiE",
  public = list(
    initialize = function(SER,L,residual_variance,estimate_residual_variance, 
        max_iter=100,tol=1e-3,track_pip=FALSE,track_lbf=FALSE) 
    {
        # initialize single effect regression models
        private$L = L
        private$SER = lapply(1:private$L, function(l) SER$clone(deep=T))
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
                fitted_l = private$SER[l].predict(d)
                d.remove_from_fitted(fitted_l)
                d.comp_residual()
                private$SER[l]$set_residual_variance(private$sigma2)
                private$SER[l]$fit(d)
                private$SER[l]$comp_kl(d)
                fitted_l = private$SER[l].predict(d)
                d.add_back_fitted(fitted_l)
            }
            b1 = self$get_posterior_b1()
            b2 = self$get_posterior_b2()
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
    coef = function() {
        d$rescale_coef(do.call(sum, get_posterior_b1()))
    },
    # some get private numbers functions
    get_prior = function() { 
        # get prior effect size, because it might be updated during iterations
        lapply(1:private$L, function(i) private$SER[l].get_prior())
    },
    get_residual_variance = function() private$sigma2,
    get_kl = function() {
        sapply(1:private$L, function(l) private$SER[l]$kl)
    },
    get_objective = function(dump = FALSE) {
        if (!all(diff(private$elbo) >= 0)) {
            warning('Objective is not non-decreasing')
            dump = TRUE
        }
        if (dump) return(private$elbo)
        else return(private$elbo[private$niter])
    },
    get_pip = function() {}, # posterior inclusion probability, L by J matrix
    get_pip_history = function() return(private$pip_history)
    get_lbf = function() {},
    get_lbf_history = function() return(private$lbf_history),
    get_posterior_b1 = function() {
        lapply(1:private$L, function(l) private$SER[l].get_posterior_b1())
    },
    get_posterior_b2 = function() {
        lapply(1:private$L, function(l) private$SER[l].get_posterior_b2())
    }
  ),
  private = list(
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
            essr = comp_expected_sum_squared_residuals(d,b1,b2) 
        } else {
            essr = private$essr
        }
        expected_loglik = comp_expected_loglik(d$get_n_sample(), private$sigma2, essr)
        return(expected_loglik - sum(self$get_kl()))
    },
    estimate_residual_variance = function(d,b1,b2) { 
        private$essr = comp_expected_sum_squared_residuals(d,b1,b2)
        private$sigma2 = private$essr / (d$get_n_sample()-1)
    },
    save_history = function() {
        private$exit()
    },
    exit = function() {
        warning("Not yet implemented")
    }
  )
)

# expected squared residuals
comp_expected_sum_squared_residuals = function(d, Eb1, Eb2) {
    if (inherits(d, c("DenseData","SparseData", "SSData")) {
        Eb1 = do.call(rbind, Eb1)
        Eb2 = do.call(rbind, Eb2)
    }
    if (inherits(d, c("DenseData","SparseData")) {
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
comp_expected_loglik = function(n, residual_variance, essr){
    -(n/2) * log(2*pi* residual_variance) - (1/(2*residual_variance)) * essr
}