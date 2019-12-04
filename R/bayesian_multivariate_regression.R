#' @title Multiviate regression object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm dmvnorm
#' @keywords internal
BayesianMultivariateRegression <- R6Class("BayesianMultivariateRegression",
  inherit = BayesianSimpleRegression,
  public = list(
    initialize = function(J, residual_variance, prior_variance, estimate_prior_variance = FALSE) {
      private$.prior_variance = prior_variance
      private$.residual_variance = residual_variance
      private$.posterior_b1 = matrix(0, J, nrow(prior_variance))
      private$to_estimate_prior_variance = estimate_prior_variance
      private$prior_variance_scale = 1
      tryCatch({
        private$.residual_variance_inv = solve(residual_variance)
      }, error = function(e) {
        warning(paste0('Cannot compute inverse for residual variance due to error:\n', e, '\nELBO computation will thus be skipped.'))
      })
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (d$Y_has_missing) stop("Cannot work with missing data in Bayesian Multivariate Regression module.")
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      bhat = XtY / d$X2_sum
      sbhat2 = lapply(1:length(d$X2_sum), function(j) private$.residual_variance / d$X2_sum[j])
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sqrt(do.call(cbind, lapply(1:length(sbhat2), function(j) diag(sbhat2[[j]]))))
      }
      # deal with prior variance: can be "estimated" across effects
      if(private$to_estimate_prior_variance) {
          if (is.null(prior_weights)) prior_weights = rep(1/private$J, private$J)
        private$prior_variance_scale = private$estimate_prior_variance(bhat,sbhat2,prior_weights)
      }
      # posterior
      post = multivariate_regression(bhat, sbhat2, private$.prior_variance * private$prior_variance_scale)
      private$.posterior_b1 = post$b1
      private$.posterior_b2 = post$b2
      private$.lbf = post$lbf
    },
    compute_loglik_null = function(d) {}
  ),
  active = list(
    residual_variance_inv = function() private$.residual_variance_inv,
    residual_variance = function(v) {
      if (missing(v)) private$.residual_variance
      else warning(paste("residual variance cannot be changed for object", "BayesianMultivariateRegression"))
    },
    prior_variance = function() private$prior_variance_scale
  ),
  private = list(
    .residual_variance_inv = NULL,
    prior_variance_scale = NULL,
    loglik = function(bhat,S,scalar,prior_weights) {
      U = private$.prior_variance * scalar
      lbf = sapply(1:length(S), function(j) dmvnorm(x = bhat[j,],sigma = S[[j]] + U,log = T) - dmvnorm(x = bhat[j,],sigma = S[[j]],log = T))
      lbf[which(is.nan(lbf))] = 0
      maxlbf = max(lbf)
      w = exp(lbf-maxlbf)
      w_weighted = w * prior_weights
      weighted_sum_w = sum(w_weighted)
      return(log(weighted_sum_w)+ maxlbf)
    },
    neg_loglik_logscale = function(lV, bhat, S, prior_weights) {
      return(-1 * private$loglik(bhat, S, exp(lV), prior_weights))
    },
    estimate_prior_variance = function(bhat, sbhat2, prior_weights) {
      # dont constrain on values of `lV` -- as a scalar it does not have to be between 0 and 1 (unlike the case with SuSiE)
      lV = optim(par=log(1), fn=private$neg_loglik_logscale, bhat=bhat, S=sbhat2, prior_weights = prior_weights, method='BFGS')$par
      V = exp(lV)
      if(private$loglik(bhat, sbhat2, 0, prior_weights) >= private$loglik(bhat, sbhat2, V, prior_weights)) V=0 # set V exactly 0 if that beats the numerical value
      return(V)
    }
  )
)

#' @title Multiviate regression calculations 
#' @importFrom abind abind
#' @keywords internal
multivariate_regression = function(bhat, S, U) {
  S_inv = lapply(1:length(S), function(j) solve(S[[j]]))
  post_cov = lapply(1:length(S), function(j) U %*% solve(diag(nrow(U)) + S_inv[[j]] %*% U))
  post_b1 = do.call(cbind, lapply(1:length(S), function(j) post_cov[[j]] %*% (S_inv[[j]] %*% bhat[j,])))
  #lbf = sapply(1:length(S), function(j) dmvnorm(x = bhat[j,],sigma = S[[j]] + U,log = T) - dmvnorm(x = bhat[j,],sigma = S[[j]],log = T))
  lbf = log(sapply(1:length(S), function(j) sqrt(det(S[[j]])/det(S[[j]]+U))*exp(0.5*t(bhat[j,])%*%S_inv[[j]]%*%post_cov[[j]]%*%S_inv[[j]]%*%bhat[j,])))
  lbf[which(is.nan(lbf))] = 0
  post_b2 = lapply(1:length(post_cov), function(j) tcrossprod(post_b1[,j]) + post_cov[[j]])
  return(list(b1 = t(post_b1), b2 = aperm(abind(post_b2, along = 3), c(2,1,3)), lbf = lbf))
}