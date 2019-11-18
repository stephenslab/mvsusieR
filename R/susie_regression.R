#' @title SUm of SIngle Effect (SuSiE) regression object
#' @importFrom R6 R6Class
#' @importFrom progress progress_bar
#' @keywords internal
SuSiE <- R6Class("SuSiE",
  public = list(
    initialize = function(SER,L,estimate_residual_variance,
    compute_objective = TRUE,max_iter=100,tol=1e-3,track_pip=FALSE,track_lbf=FALSE)
    {
        if (!compute_objective) {
            track_pip = TRUE
            estimate_residual_variance = FALSE
        }
        # initialize single effect regression models
        private$L = L
        private$to_estimate_residual_variance = estimate_residual_variance
        private$to_compute_objective = compute_objective
        private$SER = lapply(1:private$L, function(l) SER$clone(deep=T))
        private$sigma2 = SER$residual_variance
        private$elbo = vector()
        private$niter = max_iter
        private$tol = tol

        if (track_pip) private$pip_history = list()
        if (track_lbf) private$lbf_history = list()
    },
    init_coef = function(coef_index, coef_value, p, r) {
        L = length(coef_index)
        if (L <= 0) stop("Need at least one non-zero effect")
        if (L > private$L) stop("Cannot initialize more effects than the current model allows")
        if (any(which(apply(coef_value,1,sum)==0))) stop("Input coef_value must be at least one non-zero item per row")
        if (L != nrow(coef_value)) stop("Inputs coef_index and coef_value must of the same length")
        if (max(coef_index)>p) stop("Input coef_index exceeds the boundary of p")
        for (i in 1:L) {
            alpha = rep(0, p)
            mu = matrix(0, p, r)
            alpha[coef_index[i]] = 1
            mu[coef_index[i], ] = coef_value[i,]
            private$SER[[i]]$pip = alpha
            private$SER[[i]]$mu = mu
        }
    },
    fit = function(d, verbose=TRUE) {
        if (verbose) pb = progress_bar$new(format = "[:spin] Iteration :iteration (diff = :delta) :elapsed",
                                    clear = TRUE, total = private$niter, show_after = .5)
        else pb = null_progress_bar$new()
        for (i in 1:private$niter) {
            private$save_history()
            fitted = d$compute_Xb(Reduce(`+`, lapply(1:private$L, function(l) private$SER[[l]]$posterior_b1)))
            d$compute_residual(fitted)
            for (l in 1:private$L) {
                d$add_to_residual(private$SER[[l]]$predict(d))
                private$SER[[l]]$residual_variance = private$sigma2
                private$SER[[l]]$fit(d)
                if (private$to_compute_objective) private$SER[[l]]$compute_kl(d)
                d$remove_from_residual(private$SER[[l]]$predict(d))
            }
            # FIXME: different logic for univariate and multivariate cases
            if (private$to_estimate_residual_variance) 
                private$estimate_residual_variance(d)
            if (private$to_compute_objective)
                private$compute_objective(d)
            private$convergence = private$check_convergence(i)
            if (private$convergence$converged) {
                private$save_history()
                pb$tick(private$niter)
                private$niter = i
                break
            }
            pb$tick(tokens = list(delta=sprintf(private$convergence$delta, fmt = '%#.1e'), iteration=i))
        }
    },
    predict = function(x) {},
    get_objective = function(dump = FALSE) {
        if (length(private$elbo) == 0) return(NULL)
        if (!all(diff(private$elbo) >= 0)) {
            warning('Objective is not non-decreasing')
            dump = TRUE
        }
        if (dump) return(private$elbo)
        else return(private$elbo[private$niter])
    },
    get_niter = function() private$niter,
    get_convergence = function() private$convergence,
    get_pip_history = function() private$pip_history,
    get_lbf_history = function() private$lbf_history
  ),
  private = list(
    to_estimate_residual_variance = NULL,
    to_compute_objective = NULL,
    L = NULL,
    SER = NULL, # Single effect regression models
    elbo = NULL, # Evidence lower bound
    null_index = NULL, # index of null effect intentially added
    niter = NULL,
    convergence = NULL,
    pip_history = NULL, # keep track of pip
    lbf_history = NULL, # keep track of lbf
    tol = NULL, # tolerance level for convergence
    sigma2 = NULL, # residual variance
    essr = NULL,
    check_convergence = function(n) {
        if (n<=1) {
            return (list(delta=Inf, converged=FALSE))
        } else {
            if (private$to_compute_objective)
                delta = private$elbo[n] - private$elbo[n-1]
            else
                delta = max(abs(private$pip_history[[n]] - private$pip_history[[n-1]]))
            return (list(delta=delta, converged=(delta < private$tol)))
        }
    },
    compute_objective = function(d) {
        if (private$to_compute_objective) {
            v_inv = private$SER[[1]]$residual_variance_inv
            if (is.null(v_inv)) {
                # univeriate case
                # FIXME: should improve the way to identify univariate vs multivariate
                if (is.null(private$essr)) {
                    essr = compute_expected_sum_squared_residuals(d,self$posterior_b1,self$posterior_b2)
                } else {
                    essr = private$essr
                }
                expected_loglik = compute_expected_loglik(d$n_sample, private$sigma2, essr)
            } else {
                expected_loglik = -(d$n_sample * d$n_condition / 2) * log(2*pi) - d$n_sample / 2 * log(det(private$SER[[1]]$residual_variance))
                # a version not expanding the math
                resid = d$Y - d$X %*% Reduce('+', lapply(1:length(private$SER), function(l) private$SER[[l]]$posterior_b1))
                E1 = sapply(1:length(private$SER), function(l) tr(v_inv %*% t(private$SER[[l]]$posterior_b1) %*% d$XtX %*% private$SER[[l]]$posterior_b1))
                E1 = tr(v_inv%*%t(resid)%*%resid) - sum(E1)
                # After expanding the math
                #E1 = tr(v_inv%*%crossprod(d$Y, d$Y)) - 2 * tr(v_inv %*% t(d$Y) %*% d$X %*% Reduce('+', lapply(1:length(private$SER), function(l) private$SER[[l]]$posterior_b1)))
                #XtX = d$XtX
                #for (l1 in 1:length(private$SER)) {
                #    for (l2 in 1:length(private$SER)) {
                #        if (l1 != l2) {
                #            E1 = E1 + tr(v_inv %*% t(private$SER[[l1]]$posterior_b1) %*% XtX %*% private$SER[[l2]]$posterior_b1)
                #        }
                #    }
                #}
                scaled_essr = -0.5 * (E1 + Reduce('+', lapply(1:length(private$SER), function(l) private$SER[[l]]$vbxxb)))
                expected_loglik = expected_loglik + scaled_essr
            }
            elbo = expected_loglik - Reduce('+', self$kl)
        } else {
            elbo = NA
        }
        private$elbo = c(private$elbo, elbo)
    },
    estimate_residual_variance = function(d) {
        if (private$to_estimate_residual_variance) {
            private$essr = compute_expected_sum_squared_residuals(d,self$posterior_b1,self$posterior_b2)
            # FIXME: should we bother with it being N - 1 (or N - 2)?
            private$sigma2 = private$essr / d$n_sample
        }
    },
    save_history = function() {
        if (!is.null(private$pip_history)) {
            private$pip_history[[length(private$pip_history) + 1]] = self$pip
        }
        if (!is.null(private$lbf_history)) {
            private$lbf_history[[length(private$lbf_history) + 1]] = self$pip
        }
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
        if (!private$to_compute_objective) return(NULL)
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
    },
    lfsr = function(v) {
        if (missing(v)) lapply(1:private$L, function(l) private$SER[[l]]$lfsr)
        else priviate$denied('lfsr')
    },
    mixture_posterior_weights = function(v) {
        if (missing(v)) lapply(1:private$L, function(l) private$SER[[l]]$mixture_posterior_weights)
        else priviate$denied('mixture_posterior_weights')
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