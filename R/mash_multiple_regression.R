#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @keywords internal
MashMultipleRegression <- R6Class("MashMultipleRegression",
  inherit = BayesianMultipleRegression,
  public = list(
    initialize = function(J, residual_variance, mash_initializer, estimate_prior_variance = FALSE) {
      private$J = J
      private$.prior_variance = mash_initializer$prior_covariance
      private$.prior_variance$xUlist = simplify2array(private$.prior_variance$xUlist)
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
    fit = function(d, precomputed_cov_matrices = NULL, use_residual = FALSE, save_summary_stats = FALSE) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      # FIXME: can this be done faster?
      bhat = diag(1/d$d) %*% XtY
      if (!is.null(precomputed_cov_matrices)) {
        sbhat = precomputed_cov_matrices$sbhat
      } else {
        # sbhat is R by R
        sigma2 = diag(private$.residual_variance)
        sbhat = sqrt(do.call(rbind, lapply(1:private$J, function(j) sigma2 / d$d[j])))
      }
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sbhat
      }
      if (private$alpha != 0 && !all(sbhat == 1)) {
        s_alpha = sbhat ^ private$alpha 
        bhat =  bhat / s_alpha
        sbhat = sbhat ^ (1 - private$alpha)
      } else {
        s_alpha = matrix(0,0,0)
      }
      bhat[which(is.nan(bhat))] = 0
      # Fit MASH model
      if (!is.null(precomputed_cov_matrices)) is_common_cov = precomputed_cov_matrices$common_sbhat
      else is_common_cov = is_mat_common(sbhat)
      # FIXME: add d$X_has_missing true or false
      # and if X has no missing data AND is_common_cov, we can do that faster version
      # where Vinv_mat and Vinv_cube are pre-specified inside MASH initializer
      # and MASH initializer should thus have a function to take input data d and residual variance
      # and output Vinv_mat and Vinv_cube
      # 1.1 compute log-likelihood matrix given current estimates
      if (is.null(precomputed_cov_matrices))
        llik_mat = mashr:::calc_lik_rcpp(t(bhat), t(sbhat), private$null_correlation, 
                                         matrix(0,0,0),
                                         private$.prior_variance$xUlist,
                                         TRUE, 
                                         is_common_cov)$data
      else
        llik_mat = mashr:::calc_lik_common_rcpp(t(bhat), 
                                         simplyfy2array(precomputed_cov_matrices$sigma_rooti),
                                         FALSE)

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
      private$.mixture_posterior_weights = mashr:::compute_posterior_weights(private$.prior_variance$pi, exp(llik_mat))
      # 3. posterior
      ## FIXME: we might not need to compute second moment at all if we do not need to estimate residual variance
      ## we can get away with checking for convergence by PIP not by ELBO
      ## but let's set report_type = 4 and compute posterior covariance for now
      if (is.null(precomputed_cov_matrices))
        post = mashr:::calc_post_rcpp(t(bhat), t(sbhat), t(s_alpha), matrix(0,0,0), 
                              private$null_correlation,
                              matrix(0,0,0), matrix(0,0,0), 
                              private$.prior_variance$xUlist,
                              t(private$.mixture_posterior_weights),
                              is_common_cov, 4)
      else
        post = mashr:::calc_post_precision_rcpp(t(bhat), t(sbhat), t(s_alpha), matrix(0,0,0), 
                              private$null_correlation,
                              matrix(0,0,0), matrix(0,0,0), 
                              precomputed_cov_matrices$Vinv,
                              precomputed_cov_matrices$U0,
                              t(private$.mixture_posterior_weights),
                              is_common_cov, 4)
      private$.posterior_b1 = post$post_mean
      # Format post_cov for degenerated case with R = 1
      # (no need for it)
      #if (ncol(private$.posterior_b1) == 1) {
      #  post$post_cov = array(post$post_cov, c(1, 1, private$J))
      #}
      private$.posterior_b2 = post$post_cov + simplify2array(lapply(1:private$J, function(i) tcrossprod(post$post_mean[i,])))
      # 4. lfsr
      private$.lfsr = ashr::compute_lfsr(post$post_neg, post$post_zero)
      # 5. loglik under the alternative
      loglik_alt = log(exp(llik_mat[,-1,drop=FALSE]) %*% (private$.prior_variance$pi[-1]/(1-private$.prior_variance$pi[1]))) + lfactors
      # 6. adjust with alpha the EE vs EZ model
      if (nrow(s_alpha) > 0) {
        private$.loglik_null = private$.loglik_null - rowSums(log(s_alpha))
        loglik_alt = loglik_alt - rowSums(log(s_alpha))
      }
      # 7. Bayes factor
      private$.lbf = loglik_alt - private$.loglik_null
    },
    compute_loglik_null = function(d) {}
  ),
  private = list(
    null_correlation = NULL,
    alpha = NULL,
    .mixture_posterior_weights = NULL,
    .lfsr = NULL
  )
)

#' @title MASH initializer object
#' @importFrom R6 R6Class
#' @importFrom mashr expand_cov
#' @keywords internal
MashInitializer <- R6Class("MashInitializer",
  public = list(
      initialize = function(Ulist, grid, prior_weights = NULL, null_weight = NULL, V = NULL, weights_tol = 1E-10, top_mixtures = 50) {
        # FIXME: need to check input
        private$R = nrow(Ulist[[1]])
        for (l in 1:length(Ulist)) {
            if (all(Ulist[[l]] == 0))
            stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
        }
        if (any(grid<=0)) stop("grid values should be greater than zero")
        xUlist = expand_cov(Ulist, grid, TRUE)
        plen = length(xUlist) - 1
        if (is.null(null_weight)) null_weight = 0
        if (is.null(prior_weights)) prior_weights = rep(1/plen, plen)
        if (length(prior_weights) != plen)
          stop(paste("Invalid prior_weights setting: expect length", plen, "but input is of length", length(prior_weights)))
        weights = c(null_weight, prior_weights)
        weights = weights / sum(weights)
        # Filter by weights lower bound
        which.comp = which(weights[-1] > weights_tol)
        which.comp = c(1, which.comp + 1)
        filtered_weights = weights[which.comp]
        xUlist = xUlist[which.comp]
        # Filter for top weights: we only keep top weights
        if (top_mixtures > 0 && top_mixtures < length(filtered_weights)) {
          which.comp = head(sort(filtered_weights[-1], index.return=T, decreasing=T)$ix, top_mixtures)
          which.comp = c(1, which.comp + 1)
          filtered_weights = weights[which.comp]
          xUlist = xUlist[which.comp]
        }
        private$xU = list(pi = filtered_weights / sum(filtered_weights), xUlist = xUlist)
        private$U = list(pi = weights, Ulist = Ulist, grid = grid, usepointmass = TRUE)
        if (is.null(V)) private$V = diag(private$R)
        else private$V = V
        private$a = 0
      },
    precompute_cov_matrices = function(d) {
      # computes constants (SVS + U)^{-1} and (SVS)^{-1} for posterior
      # and sigma_rooti for likelihooods
      # output of this function will provide input to `mashr`'s
      # functions calc_lik_common_rcpp() and
      # calc_post_precision_rcpp()
      # The input should be sbhat data matrix
      # FIXME: currently only allows for common sbhat (see issue #5)
      sigma2 = diag(private$.residual_variance)
      sbhat = sqrt(do.call(rbind, lapply(1:private$J, function(j) sigma2 / d$d[j])))
      if (!is_mat_common(sbhat))
        stop("Input summary statistics must have the same standard error for all variables (will be the case if the original X matrix is standardized)")
      svs = sbhat[1,] * t(private$V * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
      # this is in preparation for some constants used in dmvnrom() for likelihood calculations
      sigma_rooti = list()
      for (i in 1:length(private$xU$xUlist)) sigma_rooti[[i]] = backsolve(muffled_chol(svs + private$xU$xUlist[[i]], pivot=T), diag(nrow(svs)))
      # this is in prepartion for some constants used in posterior calculation
      Vinv = solve(svs)
      U0 = list()
      for (i in 1:length(private$xU$xUlist)) U0[[i]] = private$xU$xUlist[[i]] %*% solve(Vinv %*% private$xU$xUlist[[i]] + diag(nrow(private$xU$xUlist[[i]])))
      return(list(Vinv = Vinv, U0 = U0, sigma_rooti = sigma_rooti, sbhat = sbhat, common_sbhat = TRUE))
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