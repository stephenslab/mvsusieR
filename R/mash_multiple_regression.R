#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @importFrom mashr mash_set_data mash
#' @keywords internal
MashMultipleRegression <- R6Class("MashMultipleRegression",
  inherit = BayesianMultipleRegression,
  public = list(
    initialize = function(J, residual_variance, prior_variance, effect_correlation = NULL, estimate_prior_variance = FALSE) {
      private$J = J
      private$.prior_variance = prior_variance
      private$.residual_variance = residual_variance
      if (is.null(effect_correlation)) {
        private$effect_correlation = diag(prior_variance$n_condition)
      } else {
        private$effect_correlation = effect_correlation
      }
      private$.posterior_b1 = matrix(0, J, prior_variance$n_condition)
      private$.posterior_b2 = matrix(0, J, prior_variance$n_condition)
      # Though possible to estimate from MASH model on given variables
      # we insist that the information should be provided beforehand
      private$estimate_prior_variance = FALSE
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
      # fit MASH model
      mash_data = mash_set_data(bhat, sbhat, V = private$effect_correlation)
      mobj = mash(mash_data, g = private$.prior_variance$dump(), fixg = TRUE, outputlevel = 3, verbose = FALSE, algorithm = 'R')
      # posterior
      private$.posterior_b1 = mobj$result$PosteriorMean
      ## FIXME: we might not need to compute second moment at all if we do not need to estimate residual variance
      ## we can get away with checking for convergence by PIP not by ELBO
      if (ncol(private$.posterior_b1) == 1) {
        mobj$result$PosteriorCov = array(mobj$result$PosteriorCov, c(1, 1, private$J))
      } 
      m2 = lapply(1:private$J, function(i) tcrossprod(mobj$result$PosteriorMean[i,])) 
      m2 = array(unlist(m2), dim = c(nrow(m2[[1]]), ncol(m2[[1]]), private$J))
      private$.posterior_b2 = mobj$result$PosteriorCov + m2
      # Bayes factor
      private$.lbf = mobj$alt_loglik - mobj$null_loglik
      # loglik under the null
      private$.loglik_null = mobj$null_loglik
    },
    compute_loglik_null = function(d) {}
  ),
  private = list(
    effect_correlation = NULL
  )
)

#' @title MASH prior object
#' @importFrom R6 R6Class
#' @keywords internal
MashPrior <- R6Class("MashPrior",
  public = list(
      initialize = function(Ulist, grid, prior_weights, null_weight = 0) {
        # FIXME: check input
        private$R = nrow(Ulist[[1]])
        for (l in 1:length(Ulist)) {
            if (sum(Ulist[[l]]) == 0) 
            stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
        }
        private$mash_g = list(pi = c(null_weight, prior_weights), Ulist = Ulist, grid = grid, usepointmass=TRUE)
      },
      dump = function() private$mash_g
  ),
  private = list(
      mash_g = NULL,
      R = NULL
  ),
  active = list(
      n_condition = function(v) private$R
  )
)