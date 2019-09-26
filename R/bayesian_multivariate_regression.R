#' @title Multiviate regression object
#' @importFrom R6 R6Class
#' @keywords internal
BayesianMultivariateRegression <- R6Class("BayesianMultivariateRegression",
  inherit = BayesianMultipleRegression,
  public = list(
    initialize = function(J, residual_variance, prior_variance, estimate_prior_variance = FALSE) {
      private$.prior_variance = prior_variance
      private$.residual_variance = residual_variance
      private$.posterior_b1 = matrix(0, J, nrow(prior_variance))
      private$.posterior_b2 = matrix(0, J, nrow(prior_variance))
      # We dont estimate prior variance at this point
      private$estimate_prior_variance = FALSE
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
      if (d$Y_has_missing()) stop("Cannot work with missing data in Bayesian Multivariate Regression module.")
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      bhat = XtY / d$d
      sbhat2 = diag(private$.residual_variance) / d$d
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sqrt(diag(sbhat2))
      }
      # posterior
      post = multivariate_regression(bhat, 1/d$d, private$.residual_variance, private$.prior_variance)
      private$.posterior_b1 = post$b1
      private$.posterior_b2 = post$b2
      # Bayes factor
      private$.lbf = post$lbf
      private$.lbf[sbhat2==Inf] == 0
    },
    compute_loglik_null = function(d) {}
  ),
  private = list(
    .residual_variance_inv = NULL
  )
)

multivariate_regression = function(bhat, xtx_inv, V, U) {
  S = lapply(1:length(xtx_inv), function(j) xtx_inv[j] * V)
  S_inv = lapply(1:length(S), function(j) solve(S[[j]]))
  post_cov = lapply(1:length(S), function(j) U %*% solve(diag(nrow(U)) + S_inv[[j]] %*% U))
  post_b1 = do.call(rbind, lapply(1:length(S), function(j) post_cov[[j]] %*% (S_inv[[j]] %*% bhat[j,])))
  bf = sapply(1:length(S), function(j) sqrt(det(S[[j]])/det(S[[j]]+U))*exp(0.5*t(bhat[j,])%*%S_inv[[j]]%*%post_cov[[j]]%*%S_inv[[j]]%*%bhat[j,]))
  post_b2 = do.call(rbind, lapply(1:length(post_cov), function(j) tcrossprod(post_b1[j,]) + diag(post_cov[[j]])))
  return(list(b1 = post_b1, b2 = post_b2, lbf = log(bf)))
}