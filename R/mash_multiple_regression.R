#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @importFrom mashr calc_post_rcpp 
#' @keywords internal
MashMultipleRegression <- R6Class("MashMultipleRegression",
  inherit = BayesianMultipleRegression,
  public = list(
    initialize = function(J, residual_variance, mash_initializer, estimate_prior_variance = FALSE) {
      private$J = J
      private$.prior_variance = mash_initializer$prior_covariance
      private$.residual_variance = residual_variance
      if (is.null(mash_initializer$effect_correlation)) {
        private$effect_correlation = diag(mash_initializer$n_condition)
      } else {
        private$effect_correlation = mash_initializer$effect_correlation
      }
      private$.posterior_b1 = matrix(0, J, mash_initializer$n_condition)
      private$.posterior_b2 = matrix(0, J, mash_initializer$n_condition)
      # Though possible to estimate from MASH model on given variables
      # we insist that the information should be provided beforehand
      private$estimate_prior_variance = FALSE
      # We do not need to compute posterior variance
      # if we do not estimate residual covariance
      private$to_compute_posterior_b2 = mash_initializer$to_compute_posterior_b2
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      # FIXME: can this be done faster?
      bhat = diag(1/d$d) %*% XtY
      # sbhat is R by R
      sigma2 = diag(private$.residual_variance)
      sbhat = sqrt(do.call(rbind, lapply(1:private$J, function(j) sigma2 / d$d[j])))
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sbhat
      }
      # compute MASH posterior
      if (private$to_compute_posterior_b2) output_type = 2
      else output_type = 1
      mash_post = calc_post_rcpp(t(bhat),t(sbhat),private$effect_correlation,
                                 matrix(0,0,0), matrix(0,0,0),
                                 simplify2array(Ulist), t(posterior_weights),
                                 is_mat_common(sbhat), 
                                 output_type)
      # posterior
      private$.posterior_b1 = mash_post$post_mean
      if (private$to_compute_posterior_b2) {
        ## FIXME: we might not need to compute second moment at all if we do not need to estimate residual variance
        ## we can get away with checking for convergence by PIP not by ELBO
        if (ncol(private$.posterior_b1) == 1) {
          mash_post$post_cov = array(mash_post$post_cov, c(1, 1, private$J))
        } 
        m2 = simplify2array(lapply(1:private$J, function(i) tcrossprod(mash_post$post_mean[i,])))
        private$.posterior_b2 = aperm(mash_post$post_cov, c(3,1,2)) + m2
      }
      # likelihoods

      # Bayes factor
      private$.lbf = mobj$alt_loglik - mobj$null_loglik
      # loglik under the null
      private$.loglik_null = mobj$null_loglik
    },
    compute_loglik_null = function(d) {}
  ),
  private = list(
    effect_correlation = NULL,
    to_compute_posterior_b2 = NULL
  )
)

#' @title MASH initializer object
#' @importFrom R6 R6Class
#' @keywords internal
MashInitializer <- R6Class("MashInitializer",
  public = list(
      initialize = function(Ulist, grid, prior_weights, null_weight = 0, V = NULL) {
        # FIXME: check input
        private$R = nrow(Ulist[[1]])
        for (l in 1:length(Ulist)) {
            if (sum(Ulist[[l]]) == 0) 
            stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
        }
        private$mash_g = list(pi = c(null_weight, prior_weights), Ulist = Ulist, grid = grid, usepointmass=TRUE)
        private$V = V 
      }
  ),
  private = list(
      mash_g = NULL,
      R = NULL,
      V = NULL
  ),
  active = list(
      n_condition = function(v) private$R,
      prior_covariance = function() private$mash_g,
      effect_correlation = function() private$V
  )
)