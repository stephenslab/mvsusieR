#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @keywords internal
MashMultipleRegression <- R6Class("MashMultipleRegression",
  inherit = BayesianMultipleRegression,
  public = list(
    initialize = function(J, residual_variance, mash_initializer, estimate_prior_variance = FALSE) {
      private$J = J
      private$.prior_variance = mash_initializer$prior_covariance
      private$.residual_variance = residual_variance
      if (is.null(mash_initializer$null_correlation)) {
        private$null_correlation = diag(mash_initializer$n_condition)
      } else {
        private$null_correlation = mash_initializer$null_correlation
      }
      private$alpha = mash_initializer$alpha
      private$.posterior_b1 = matrix(0, J, mash_initializer$n_condition)
      private$.posterior_b2 = matrix(0, J, mash_initializer$n_condition)
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
      if (private$alpha != 0 && !all(sbhat == 1)) {
        s_alpha = sbhat ^ private$alpha 
        bhat =  bhat / s_alpha
        sbhat = sbhat ^(1 - private$alpha)
      } else {
        s_alpha = matrix(0,0,0)
      }
      bhat[which(is.nan(bhat))] = 0
      # Fit MASH model
      is_common_cov = is_mat_common(sbhat)
      # 1.1 compute log-likelihood matrix given current estimates
      llik_mat = mashr:::calc_lik_rcpp(t(bhat), t(sbhat), private$null_correlation,
                             matrix(0,0,0), private$.prior_variance$xUlist, TRUE,
                             is_common_cov)$data
      # 1.2 give a warning if any columns have -Inf likelihoods.
      rows <- which(apply(llik_mat,2,function (x) any(is.infinite(x))))
      if (length(rows) > 0)
        warning(paste("Some mixture components result in non-finite likelihoods,",
                          "either\n","due to numerical underflow/overflow,",
                          "or due to invalid covariance matrices",
                          paste(rows,collapse=", "), "\n"))
      private$.loglik_null = llik_mat[,1]
      # 1.3 get relative loglik
      lfactors = apply(llik_mat,1,max)
      llik_mat = llik_mat - lfactors
      # 2. compute posterior weights
      posterior_weights = mashr:::compute_posterior_weights(private$.prior_variance$pi, exp(llik_mat))
      # 3. posterior
      ## FIXME: we might not need to compute second moment at all if we do not need to estimate residual variance
      ## we can get away with checking for convergence by PIP not by ELBO
      ## but let's set report_type = 4 and compute posterior covariance for now
      post = mashr:::calc_post_rcpp(t(bhat), t(sbhat), t(s_alpha), matrix(0,0,0), 
                            private$null_correlation,
                            matrix(0,0,0), matrix(0,0,0), 
                            private$.prior_variance$xUlist, t(posterior_weights),
                            is_common_cov, 4)
      private$.posterior_b1 = post$post_mean
      if (ncol(private$.posterior_b1) == 1) {
        post$post_cov = array(post$post_cov, c(1, 1, private$J))
      } 
      m2 = simplify2array(lapply(1:private$J, function(i) tcrossprod(post$post_mean[i,])))
      private$.posterior_b2 = aperm(post$post_cov, c(3,1,2)) + m2
      # 4. loglik under the alternative
      loglik_alt = log(exp(llik_mat[,-1,drop=FALSE]) %*% (private$.prior_variance$pi[-1]/(1-private$.prior_variance$pi[1]))) + lfactors
      # 5. adjust with alpha the EE vs EZ model
      if (nrow(s_alpha) > 0) {
        private$.loglik_null = private$.loglik_null - rowSums(log(s_alpha))
        loglik_alt = loglik_alt - rowSums(log(s_alpha))
      }
      # 6. Bayes factor
      private$.lbf = loglik_alt - private$.loglik_null
    },
    compute_loglik_null = function(d) {}
  ),
  private = list(
    null_correlation = NULL,
    alpha = NULL
  )
)

#' @title MASH initializer object
#' @importFrom R6 R6Class
#' @importFrom mashr expand_cov
#' @keywords internal
MashInitializer <- R6Class("MashInitializer",
  public = list(
      initialize = function(Ulist, grid, prior_weights, null_weight = 0, V = NULL) {
        # FIXME: need to check input
        private$R = nrow(Ulist[[1]])
        for (l in 1:length(Ulist)) {
            if (sum(Ulist[[l]]) == 0) 
            stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
        }
        xUlist = expand_cov(Ulist, grid, TRUE)
        weights = c(null_weight, prior_weights)
        which.comp = which(weights[-1] > 1e-10)
        which.comp = c(1, which.comp + 1)
        private$xU = list(pi = weights[which.comp], xUlist = simplify2array(xUlist[which.comp]))
        private$U = list(pi = weights, Ulist = Ulist, grid = grid, usepointmass = TRUE)
        private$V = V
        private$a = 0
      }
  ),
  private = list(
      R = NULL,
      V = NULL,
      U = NULL,
      xU = NULL,
      a = NULL
  ),
  active = list(
      n_condition = function(v) private$R,
      prior_covariance = function() private$xU,
      mash_prior = function() private$U,
      null_correlation = function() private$V,
      alpha = function(value) {
        if (missing(value)) return(private$a)
        else private$a = value
      }
  )
)