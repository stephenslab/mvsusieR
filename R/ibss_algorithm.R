#' @title SUm of SIngle Effect (SuSiE) regression IBSS algorithm
#' @importFrom R6 R6Class
#' @importFrom progress progress_bar
#' @keywords internal
SuSiE <- R6Class("SuSiE",
  public = list(
    initialize = function(SER,L,estimate_residual_variance=FALSE,
                    compute_objective = TRUE,max_iter=100,tol=1e-3,
                    track_pip=FALSE,track_lbf=FALSE, track_prior_est=FALSE)
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
        private$elbo = vector()
        private$.niter = max_iter
        private$tol = tol
        if (track_pip) private$.pip_history = list()
        if (track_lbf) private$.lbf_history = list()
        if (track_prior_est) private$.prior_history = list()
    },
    init_from = function(model) {
        mu = model$b1
        if (is.null(dim(mu)) || length(dim(mu)) != 3)
            stop("Input coefficient should be a 3-D array for first moment estimate of each single effect.")
        L = dim(mu)[1]
        if (L != private$L)
            stop("Cannot initialize different number of effects than the current model allows")
        for (i in 1:L) {
            mu_l = mu[i,,]/model$alpha[i,]
            mu_l[which(is.nan(mu_l) | !is.finite(mu_l))] = 0
            private$SER[[i]]$mu = mu_l
            private$SER[[i]]$pip = model$alpha[i,]
            if (!is.null(model$V)) private$SER[[i]]$prior_variance = model$V[i]
        }
    },
    fit = function(d, prior_weights=NULL, estimate_prior_variance_method=NULL, check_null_threshold=0, verbosity=1) {
        if (verbosity==1) pb = progress_bar$new(format = "[:spin] Iteration :iteration (diff = :delta) :elapsed",
                                    clear = TRUE, total = private$.niter, show_after = .5)
        else pb = null_progress_bar$new()
        if (verbosity > 1) message("Running IBSS algorithm ...")
        for (i in 1:private$.niter) {
            private$save_history()
            fitted = d$compute_Xb(Reduce(`+`, lapply(1:private$L, function(l) private$SER[[l]]$posterior_b1)))
            d$compute_residual(fitted)
            for (l in 1:private$L) {
                d$add_to_residual(private$SER[[l]]$predict(d))
                # For the first 10 iterations, don't do the check with zero when EM updates are used to estimate prior variance
                # see https://github.com/stephenslab/mvsusieR/issues/26#issuecomment-612947198
                private$SER[[l]]$fit(d, prior_weights=prior_weights, estimate_prior_variance_method=estimate_prior_variance_method,
                                    check_null_threshold=ifelse(!is.null(estimate_prior_variance_method) && estimate_prior_variance_method == "EM" && i <= 10, NA, check_null_threshold))
                if (private$to_compute_objective) private$SER[[l]]$compute_kl(d)
                d$remove_from_residual(private$SER[[l]]$predict(d))
            }
            if (private$to_compute_objective)
              private$compute_objective(d)
            private$.convergence = private$check_convergence(i)
            if (private$.convergence$converged) {
                private$save_history()
                pb$tick(private$.niter)
                private$.niter = i
                private$add_back_zero_effects()
                break
            }
            if (private$to_estimate_residual_variance){
                d$set_residual_variance(private$estimate_residual_variance(d),
                                        quantities = c('residual_variance', 'effect_variance'))
            }
            if (!is.null(estimate_prior_variance_method) && estimate_prior_variance_method == "EM")
                private$trim_zero_effects()
            if (i == private$.niter) {
                warning(paste("IBSS failed to converge after", i, "iterations. Perhaps you should increase max_iter and try again."))
                private$add_back_zero_effects()
            }
            if (verbosity>1) {
               message(paste("Iteration", i, "delta =",
                            private$.convergence$delta))
            } else {
              pb$tick(tokens = list(delta=sprintf(private$.convergence$delta, fmt = '%#.1e'), iteration=i))
            }
        }
    },
    get_objective = function(dump = FALSE, warning_tol = 1E-6) {
        if (length(private$elbo) == 0) return(NA)
        if (!all(diff(private$elbo) >= (-1 * warning_tol))) {
            warning('Objective is not non-decreasing')
            dump = TRUE
        }
        if (dump) return(private$elbo)
        else return(private$elbo[private$.niter])
    }
  ),
  active = list(
    niter = function() private$.niter,
    convergence = function() private$.convergence,
    # get prior effect size, because it might be updated during iterations
    prior_variance = function() sapply(1:private$L, function(l) private$SER[[l]]$prior_variance),
    kl = function() {
        if (!private$to_compute_objective) NA
        else sapply(1:private$L, function(l) private$SER[[l]]$kl)
    },
    # posterior inclusion probability, J by L matrix
    pip = function() do.call(cbind, lapply(1:private$L, function(l) private$SER[[l]]$pip)),
    lbf = function() sapply(1:private$L, function(l) private$SER[[l]]$lbf),
    pip_history = function() lapply(ifelse(private$.niter>1, 2, 1):length(private$.pip_history), function(i) private$.pip_history[[i]]),
    lbf_history = function() lapply(ifelse(private$.niter>1, 2, 1):length(private$.lbf_history), function(i) private$.lbf_history[[i]]),
    prior_history = function() lapply(ifelse(private$.niter>1, 2, 1):length(private$.prior_history), function(i) private$.prior_history[[i]]),
    posterior_b1 = function() lapply(1:private$L, function(l) private$SER[[l]]$posterior_b1),
    posterior_b2 = function() lapply(1:private$L, function(l) private$SER[[l]]$posterior_b2),
    clfsr = function() lapply(1:private$L, function(l) private$SER[[l]]$lfsr),
    mixture_posterior_weights = function() lapply(1:private$L, function(l) private$SER[[l]]$mixture_posterior_weights)
  ),
  private = list(
    to_estimate_residual_variance = NULL,
    to_compute_objective = NULL,
    L = NULL,
    SER = NULL, # Single effect regression models
    SER_NULL = list(), # Single effect regression models removed
    elbo = NULL, # Evidence lower bound
    null_index = NULL, # index of null effect intentially added
    .niter = NULL,
    .convergence = NULL,
    .pip_history = NULL, # keep track of pip
    .lbf_history = NULL, # keep track of lbf
    .prior_history = NULL, # keep track of prior estimates
    tol = NULL, # tolerance level for convergence
    essr = NULL,
    check_convergence = function(n) {
        if (n<=1) {
            return (list(delta=Inf, converged=FALSE))
        } else {
            if (private$to_compute_objective)
                delta = private$elbo[n] - private$elbo[n-1]
            else
                # convergence check with marginal PIP
                delta = max(abs(apply(1 - private$.pip_history[[n]], 1, prod) - apply(1 - private$.pip_history[[n-1]], 1, prod)))
            return (list(delta=delta, converged=(delta < private$tol)))
        }
    },
    compute_objective = function(d) {
        if (private$to_compute_objective) {
            if (is.matrix(d$residual_variance)) {
                expected_loglik = private$compute_expected_loglik_multivariate(d)
            } else {
                expected_loglik = private$compute_expected_loglik_univariate(d)
            }
            elbo = expected_loglik - Reduce('+', self$kl)
        } else {
            elbo = NA
        }
        private$elbo = c(private$elbo, elbo)
    },
    # expected loglikelihood for SuSiE model
    compute_expected_loglik_univariate = function(d) {
        if(d$Y_has_missing){
          Y_missing_assign =  table(d$Y_missing_pattern_assign)
          expected_loglik = -0.5 * log(2*pi) * sum(sapply(d$residual_variance_eigenvalues, length) * Y_missing_assign) -
            0.5 * sum(sapply(d$residual_variance_eigenvalues, function(x) ifelse(length(x)>0,sum(log(x)),0)) * Y_missing_assign)
          essr = private$compute_expected_sum_squared_residuals_univariate(d)
          return(expected_loglik - 0.5 * essr)
        }else{
          n = d$n_sample
          residual_variance = d$residual_variance
          essr = private$compute_expected_sum_squared_residuals_univariate(d)
          return(-(n/2) * log(2*pi* residual_variance) - (1/(2*residual_variance)) * essr)
        }
    },
    compute_expected_loglik_multivariate = function(d) {
      if(d$Y_has_missing){
        Y_missing_assign =  table(d$Y_missing_pattern_assign)
        expected_loglik = -0.5 * log(2*pi) * sum(sapply(d$residual_variance_eigenvalues, length) * Y_missing_assign) -
          0.5 * sum(sapply(d$residual_variance_eigenvalues, function(x) ifelse(length(x)>0,sum(log(x)),0)) * Y_missing_assign)
        essr = private$compute_expected_sum_squared_residuals_multivariate(d)
      }else{
        expected_loglik = -(d$n_sample * d$n_condition / 2) * log(2*pi) - d$n_sample / 2 * log(det(d$residual_variance))
        essr = private$compute_expected_sum_squared_residuals_multivariate(d)
      }
      return(expected_loglik - 0.5 * essr)
    },
    estimate_residual_variance = function(d) {
        if (is.matrix(d$residual_variance)) {
            # FIXME: to implement estimating a vector of length R, or even a scalar
            E1 = lapply(1:length(private$SER), function(l) crossprod(private$SER[[l]]$posterior_b1, d$XtX %*% private$SER[[l]]$posterior_b1))
            E1 = crossprod(d$residual) - Reduce('+', E1)
            return((E1 + Reduce('+', lapply(1:length(private$SER), function(l) private$SER[[l]]$bxxb))) / d$n_sample)
        } else {
            return(private$compute_expected_sum_squared_residuals_univariate(d) / d$n_sample)
        }
    },
    # expected squared residuals
    compute_expected_sum_squared_residuals_univariate = function(d) {
      Eb1 = t(do.call(cbind, self$posterior_b1))
      Eb2 = t(do.call(cbind, self$posterior_b2))
      if (inherits(d, "RSSData")) {
        # RSSData is inherited from DenseData
        # actually code below will also work for DenseData
        # that is why there is no need to treat them separately in multivarite computation
        # However, computational complexity may be different (depending on the dimension of XtX)
        XB2 = sum((Eb1 %*% d$XtX) * Eb1)
        return(as.numeric(crossprod(d$residual) - XB2 + sum(d$X2_sum*t(Eb2))))
      } else if (inherits(d, "SSData")){
        XB2 = sum((Eb1 %*% d$XtX) * Eb1)
        EB = colSums(Eb1)
        return(as.numeric(d$YtY - 2*sum(EB * d$XtY) + sum(EB * d$compute_MXt(EB)) -
          XB2 + sum(d$X2_sum * t(Eb2))))
      }else {
        # full data, DenseData object
        if(d$Y_has_missing){
          resid_var_inv = unlist(d$residual_variance_inv)[d$Y_missing_pattern_assign]
          E1 = sapply(1:length(private$SER), function(l){
            Xb = d$compute_Xb(private$SER[[l]]$posterior_b1)
            sum(Xb^2 * resid_var_inv)
          })
          E1 = sum(d$residual^2 * resid_var_inv) - sum(E1)
          return(E1 + Reduce('+', lapply(1:length(private$SER), function(l) private$SER[[l]]$vbxxb)))
        }else{
          Xr = d$compute_MXt(Eb1)
          Xrsum = colSums(Xr)
          return(sum((d$Y-Xrsum)^2) - sum(Xr^2) + sum(d$X2_sum*t(Eb2)))
        }
      }
    },
    compute_expected_sum_squared_residuals_multivariate = function(d) {
      if(d$Y_has_missing){
        E1 = sapply(1:length(private$SER), function(l){
          Xb = d$compute_Xb(private$SER[[l]]$posterior_b1)
          sum(sapply(1:d$n_sample, function(i) crossprod(Xb[i,],
                                                         d$residual_variance_inv[[d$Y_missing_pattern_assign[i]]] %*% Xb[i, ])  ))
        })
        E1 = sum(sapply(1:d$n_sample, function(i) crossprod(d$residual[i,],
                                                            d$residual_variance_inv[[d$Y_missing_pattern_assign[i]]] %*% d$residual[i,]) )) - sum(E1)
      }else if(inherits(d, "SSData")){
        v_inv = d$residual_variance_inv
        E1 = sapply(1:length(private$SER), function(l) tr(v_inv %*% t(private$SER[[l]]$posterior_b1) %*% d$XtX %*% private$SER[[l]]$posterior_b1))
        Eb1 = aperm(abind::abind(lapply(1:private$L, function(l) private$SER[[l]]$posterior_b1),along=3), c(3,1,2))
        if(dim(Eb1)[1] == 1){
          Eb1 = Eb1[1,,]
        }else{
          Eb1 = do.call(cbind, lapply(1:dim(Eb1)[3], function(i) colSums(Eb1[,,i]))) # J by R
        }
        E2 = crossprod(Eb1, d$XtY)
        E3 = crossprod(Eb1,d$XtX) %*% Eb1
        E1 = tr(v_inv%*%(d$YtY - E2 - t(E2) + E3)) - sum(E1)
      }else{
        v_inv = d$residual_variance_inv
        E1 = sapply(1:length(private$SER), function(l) tr(v_inv %*% t(private$SER[[l]]$posterior_b1) %*% d$XtX %*% private$SER[[l]]$posterior_b1))
        E1 = tr(v_inv%*%crossprod(d$residual)) - sum(E1)
      }
      return(E1 + Reduce('+', lapply(1:length(private$SER), function(l) private$SER[[l]]$vbxxb)))
    },
    trim_zero_effects = function() {
        # remove single effect models where estimated prior is zero
        if (length(private$SER) > 1) {
            zero_idx = which(sapply(1:length(private$SER), function(i) private$SER[[i]]$prior_variance) == 0)
            if (length(zero_idx) == length(private$SER)) zero_idx = zero_idx[2:length(zero_idx)]
            if (length(zero_idx)) {
                private$SER_NULL = c(private$SER_NULL, private$SER[zero_idx])
                private$SER = private$SER[-zero_idx]
                private$L = length(private$SER)
            }
        }
    },
    add_back_zero_effects = function() {
        # keep output number of effects consistent with specified L
        if (length(private$SER_NULL)) {
            private$SER = c(private$SER, private$SER_NULL)
            private$L = length(private$SER)
        }
    },
    save_history = function() {
        if (!is.null(private$.pip_history)) {
            private$.pip_history[[length(private$.pip_history) + 1]] = self$pip
        }
        if (!is.null(private$.lbf_history)) {
            private$.lbf_history[[length(private$.lbf_history) + 1]] = self$lbf
        }
        if (!is.null(private$.prior_history)) {
            private$.prior_history[[length(private$.prior_history) + 1]] = sapply(1:private$L, function(l) ifelse(is.null(private$SER[[l]]$prior_variance_scalar), private$SER[[l]]$prior_variance, private$SER[[l]]$prior_variance_scalar))
        }
    }
  )
)
