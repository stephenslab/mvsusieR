#' @title Multiviate regression object
#' @importFrom R6 R6Class
#' @keywords internal
BayesianMultivariateRegression <- R6Class("BayesianMultivariateRegression",
  inherit = BayesianSimpleRegression,
  public = list(
    initialize = function(J, residual_variance, prior_variance) {
      private$J = J
      private$.prior_variance = prior_variance
      private$.posterior_b1 = matrix(0, J, nrow(prior_variance))
      private$prior_variance_scale = 1
      self$residual_variance = residual_variance
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE, save_var = FALSE, estimate_prior_variance_method = NULL, check_null_threshold=0) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      # X2_sum is either a length J vector or J by R by R array
      if(d$Y_has_missing){
        bhat = t(sapply(1:d$n_effect, function(j) solve(d$X2_sum[j,,], XtY[j,])))
        sbhat2 = lapply(1:d$n_effect, function(j) invert_via_chol(d$X2_sum[j,,]))
      }else{
        bhat = XtY / d$X2_sum
        sbhat2 = lapply(1:length(d$X2_sum), function(j) private$.residual_variance / d$X2_sum[j])
      }
      bhat[which(is.nan(bhat))] = 0
      for (j in 1:length(sbhat2)) {
        sbhat2[[j]][which(is.nan(sbhat2[[j]]) | is.infinite(sbhat2[[j]]))] = 1E6
      }
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sqrt(do.call(rbind, lapply(1:length(sbhat2), function(j) diag(sbhat2[[j]]))))
      }
      # deal with prior variance: can be "estimated" across effects
      if(!is.null(estimate_prior_variance_method)) {
        if (estimate_prior_variance_method == "EM") {
          private$cache = list(b=bhat, s=sbhat2, update_scale=T)
        } else {
          private$prior_variance_scale = private$estimate_prior_variance(bhat,sbhat2,prior_weights,method=estimate_prior_variance_method,check_null_threshold=check_null_threshold)
        }
      }
      # posterior
      post = multivariate_regression(bhat, sbhat2, private$.prior_variance * private$prior_variance_scale)
      private$.posterior_b1 = post$b1
      private$.posterior_b2 = post$b2
      if (save_var) private$.posterior_variance = post$cov
      private$.lbf = post$lbf
    }
  ),
  active = list(
    residual_variance_inv = function() private$.residual_variance_inv,
    residual_variance = function(v) {
      if (missing(v)) private$.residual_variance
      else {
        private$.residual_variance = v
        tryCatch({
          private$.residual_variance_inv = invert_via_chol(v)
        }, error = function(e) {
          warning(paste0('Cannot compute inverse for residual variance due to error:\n', e, '\nELBO computation will thus be skipped.'))
        })
      }
    },
    prior_variance = function(v) {
      if (missing(v)) private$prior_variance_scale
      else private$prior_variance_scale = v
    }
  ),
  private = list(
    .residual_variance_inv = NULL,
    .prior_variance_inv = NULL,
    prior_variance_scale = NULL,
    loglik = function(scalar, bhat, S, prior_weights) {
      U = private$.prior_variance * scalar
      lbf = multivariate_lbf(bhat, S, U)
      return(compute_softmax(lbf, prior_weights)$log_sum)
    },
    estimate_prior_variance_optim = function(betahat, shat2, prior_weights, ...) {
      # log(1) = 0
      lV = optim(par=0, fn=private$neg_loglik_logscale, betahat=betahat, shat2=shat2, prior_weights = prior_weights, ...)$par
      return(exp(lV))
    },
    estimate_prior_variance_em_direct_inv = function(pip, inv_function = invert_via_chol) {
      # Update directly using inverse of prior matrix
      # This is very similar to updating the univariate case via EM,
      # \sigma_0^2 = \mathrm{tr}(S_0^{-1} E[bb^T])/r
      # where S_0 is prior variance, E[bb^T] is 2nd moment of SER effect:
      # that is, E[bb^T] = \sum_j alpha_j * b_jb_j^T where b_j is posterior of j
      # Recall in univariate case it is \sigma_0^2 = E[bb^T] directly
      if (is.null(private$.prior_variance_inv)) private$.prior_variance_inv = inv_function(private$.prior_variance)
      if (length(dim(private$.posterior_b2)) == 3) {
        # when R > 1
        mu2 = Reduce("+", lapply(1:length(pip), function(j) pip[j] * private$.posterior_b2[,,j]))
      } else {
        # when R = 1 each post_b2 is a scalar.
        # Now make it a matrix to be compatable with later computations.
        if (ncol(private$.posterior_b2) != 1) stop("Data dimension is incorrect for posterior_b2")
        mu2 = matrix(sum(pip * private$.posterior_b2[,1]), 1,1)
      }
      V = sum(diag(private$.prior_variance_inv %*% mu2)) / nrow(private$.prior_variance)
      return(V)
    },
    estimate_prior_variance_em_inv_safe = function(pip) {
      # Instead of computing S_0^{-1} and E[bb^T] we compute them as one quantity to avoid explicit inverse
      # We need S_inv a J vector of R by R matrices (private$cache$s), bhat a J by R vector (private$cache$b),
      # the original prior matrix S_0 (private$.prior_variance)
      # and the scalar from previous update (private$prior_variance_scale)
      # U = \sigma_0 S_0
      U = private$prior_variance_scale * private$.prior_variance
      S_inv = lapply(1:private$J, function(j) invert_via_chol(private$cache$s[[j]]))
      # posterior covariance pre-multipled by U^{-1}
      post_cov_U = lapply(1:private$J, function(j) solve(diag(nrow(U)) + S_inv[[j]] %*% U))
      # posterior first moment pre-multipled by U^{-1}
      post_b1_U = lapply(1:private$J, function(j) post_cov_U[[j]] %*% (S_inv[[j]] %*% private$cache$b[j,]))
      # posterior 2nd moment pre-multiplied by S_0^{-1}
      b2_U = lapply(1:private$J, function(j) private$prior_variance_scale * (tcrossprod(post_b1_U[[j]]) %*% U + post_cov_U[[j]]))
      V = sum(diag(Reduce("+", lapply(1:private$J, function(j) pip[j] * b2_U[[j]])))) / nrow(U)
      return(V)
    },
    estimate_prior_variance_em = function(pip) {
      tryCatch({
          return(private$estimate_prior_variance_em_direct_inv(pip, inv_function = invert_via_chol))
        },
        error = function(e) {
          return(private$estimate_prior_variance_em_inv_safe(pip))
        })
    },
    estimate_prior_variance_simple = function() 1
  )
)

#' @title Multiviate regression calculations
#' @importFrom abind abind
#' @keywords internal
multivariate_regression = function(bhat, S, U) {
  # FIXME: this can be pre-computed to save some computations
  S_inv = lapply(1:length(S), function(j) invert_via_chol(S[[j]]))
  post_cov = lapply(1:length(S), function(j) U %*% solve(diag(nrow(U)) + S_inv[[j]] %*% U))
  lbf = sapply(1:length(S), function(j) 0.5 * (log(det(S[[j]])) - log(det(S[[j]]+U))) + 0.5*t(bhat[j,])%*%S_inv[[j]]%*%post_cov[[j]]%*%S_inv[[j]]%*%bhat[j,])
  lbf[which(is.nan(lbf))] = 0
  # lbf = multivariate_lbf(bhat, S, U)
  # using rbind here will end up with dimension issues for degenerated case on J; have to use t(...(cbind, )) instead
  post_b1 = t(do.call(cbind, lapply(1:length(S), function(j) post_cov[[j]] %*% (S_inv[[j]] %*% bhat[j,]))))
  post_b2 = lapply(1:length(post_cov), function(j) tcrossprod(post_b1[j,]) + post_cov[[j]])
  # deal with degerated case with 1 condition
  if (ncol(post_b1) == 1) {
    post_b2 = matrix(unlist(post_b2), length(post_b2), 1)
  } else {
    post_b2 = aperm(abind(post_b2, along = 3), c(2,1,3))
  }
  return(list(b1 = post_b1, b2 = post_b2, lbf = lbf, cov = post_cov))
}

#' @title Multiviate logBF
#' @importFrom mvtnorm dmvnorm
#' @keywords internal
multivariate_lbf = function(bhat, S, U) {
  lbf = sapply(1:length(S), function(j) dmvnorm(x = bhat[j,],sigma = S[[j]] + U,log = T) - dmvnorm(x = bhat[j,],sigma = S[[j]],log = T))
  lbf[which(is.nan(lbf))] = 0
  return(lbf)
}