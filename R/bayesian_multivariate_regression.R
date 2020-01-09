#' @title Multiviate regression object
#' @importFrom R6 R6Class
#' @keywords internal
BayesianMultivariateRegression <- R6Class("BayesianMultivariateRegression",
  inherit = BayesianSimpleRegression,
  public = list(
    initialize = function(J, residual_variance, prior_variance) {
      private$J = J
      private$.prior_variance = prior_variance
      private$.residual_variance = residual_variance
      private$.posterior_b1 = matrix(0, J, nrow(prior_variance))
      private$prior_variance_scale = 1
      tryCatch({
        private$.residual_variance_inv = invert_via_chol(residual_variance)
      }, error = function(e) {
        warning(paste0('Cannot compute inverse for residual variance due to error:\n', e, '\nELBO computation will thus be skipped.'))
      })
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE, estimate_prior_variance_method = NULL) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (d$Y_has_missing) stop("Cannot work with missing data in Bayesian Multivariate Regression module.")
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      bhat = XtY / d$X2_sum
      bhat[which(is.nan(bhat))] = 0
      if (d$Y_has_missing) sbhat2 = lapply(1:nrow(d$X2_sum), function(j) private$.residual_variance / d$X2_sum[j,])
      else sbhat2 = lapply(1:length(d$X2_sum), function(j) private$.residual_variance / d$X2_sum[j])
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sqrt(do.call(rbind, lapply(1:length(sbhat2), function(j) diag(sbhat2[[j]]))))
        private$.sbhat[which(is.nan(private$.sbhat) | is.infinite(private$.sbhat))] = 1E3
      }
      if (d$Y_has_missing) stop("Computation involving missing data in Y has not been implemented in BayesianMultivariateRegression method.")
      # deal with prior variance: can be "estimated" across effects
      if(!is.null(estimate_prior_variance_method)) {
          if (is.null(prior_weights)) prior_weights = rep(1/private$J, private$J)
        private$prior_variance_scale = private$estimate_prior_variance(bhat,sbhat2,prior_weights,method=estimate_prior_variance_method)
      }
      # posterior
      post = multivariate_regression(bhat, sbhat2, private$.prior_variance * private$prior_variance_scale)
      private$.posterior_b1 = post$b1
      private$.posterior_b2 = post$b2
      private$.lbf = post$lbf
    }
  ),
  active = list(
    residual_variance_inv = function() private$.residual_variance_inv,
    residual_variance = function(v) {
      if (missing(v)) private$.residual_variance
      else{
        private$.residual_variance = v
        private$.residual_variance_inv = invert_via_chol(v)
        }
    },
    prior_variance = function() private$prior_variance_scale
  ),
  private = list(
    .residual_variance_inv = NULL,
    prior_variance_scale = NULL,
    loglik = function(bhat,S,scalar,prior_weights) {
      U = private$.prior_variance * scalar
      lbf = multivariate_lbf(bhat, S, U)
      return(compute_weighted_sum(lbf, prior_weights)$log_sum)
    },
    neg_loglik_logscale = function(lV, bhat, S, prior_weights) {
      return(-1 * private$loglik(bhat, S, exp(lV), prior_weights))
    },
    estimate_prior_variance = function(bhat, sbhat2, prior_weights, method=c('optim','simple')) {
      if (method == 'optim') {
        # dont constrain on values of `lV` -- as a scalar it does not have to be between 0 and 1 (unlike the case with SuSiE)
        lV = optim(par=log(1), fn=private$neg_loglik_logscale, bhat=bhat, S=sbhat2, prior_weights = prior_weights, method='BFGS')$par
        V = exp(lV)
      } else {
        # just use 1 the default, and to be compared with 0 below
        V = 1
      }
      if(private$loglik(bhat, sbhat2, 0, prior_weights) >= private$loglik(bhat, sbhat2, V, prior_weights)) V=0 # set V exactly 0 if that beats the numerical value
      return(V)
    }
  )
)

#' @title Multiviate regression calculations
#' @importFrom abind abind
#' @keywords internal
multivariate_regression = function(bhat, S, U) {
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
  return(list(b1 = post_b1, b2 = post_b2, lbf = lbf))
}

#' @title Multiviate logBF
#' @importFrom mvtnorm dmvnorm
#' @keywords internal
multivariate_lbf = function(bhat, S, U) {
  lbf = sapply(1:length(S), function(j) dmvnorm(x = bhat[j,],sigma = S[[j]] + U,log = T) - dmvnorm(x = bhat[j,],sigma = S[[j]],log = T))
  lbf[which(is.nan(lbf))] = 0
  return(lbf)
}