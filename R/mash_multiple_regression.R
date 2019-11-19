#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @keywords internal
MashMultipleRegression <- R6Class("MashMultipleRegression",
  inherit = BayesianMultipleRegression,
  public = list(
    initialize = function(J, residual_variance, mash_initializer, estimate_prior_variance = FALSE) {
      private$.prior_variance = mash_initializer$prior_covariance
      private$.prior_variance$xUlist = simplify2array(private$.prior_variance$xUlist)
      private$.residual_variance = residual_variance
      # FIXME: not sure if this is the best way to handle
      tryCatch({
        private$.residual_variance_inv = solve(residual_variance)
      }, error = function(e) {
        warning(paste0('Cannot compute inverse for residual variance due to error:\n', e, '\nELBO computation will thus be skipped.'))
      })
      if (is.matrix(residual_variance)) private$residual_correlation = cov2cor(residual_variance)
      else private$residual_correlation = diag(1) 
      private$alpha = mash_initializer$alpha
      private$precomputed_cov_matrices = mash_initializer$precomputed
      private$.posterior_b1 = matrix(0, J, mash_initializer$n_condition)
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
      bhat = XtY / d$d
      if (!is.null(private$precomputed_cov_matrices)) {
        sbhat = private$precomputed_cov_matrices$sbhat
      } else {
        # sbhat is R by R
        # for non-missing Y d$d is a J vector
        # for missing Y it is a J by R matrix
        sigma2 = diag(private$.residual_variance)
        if (d$Y_has_missing()) sbhat = sqrt(do.call(rbind, lapply(1:nrow(d$d), function(j) sigma2 / d$d[j,])))
        else sbhat = sqrt(do.call(rbind, lapply(1:length(d$d), function(j) sigma2 / d$d[j])))
        sbhat[which(is.nan(sbhat) | is.infinite(sbhat))] = 1E6
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
      if (!is.null(private$precomputed_cov_matrices)) is_common_cov = private$precomputed_cov_matrices$common_sbhat
      else is_common_cov = is_mat_common(sbhat)
      # 1.1 compute log-likelihood matrix given current estimates
      if (is.null(private$precomputed_cov_matrices) || ncol(bhat) == 1) {
        llik_mat = mashr:::calc_lik_rcpp(t(bhat), t(sbhat), private$residual_correlation, 
                                         matrix(0,0,0),
                                         private$.prior_variance$xUlist,
                                         TRUE,
                                         is_common_cov)$data
      } else {
        llik_mat = mashr:::calc_lik_rooti_rcpp(t(bhat), 
                                         private$precomputed_cov_matrices$sigma_rooti,
                                         TRUE,
                                         is_common_cov)$data
      }

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
      if (is.null(private$precomputed_cov_matrices) || ncol(bhat) == 1) {
        post = mashr:::calc_post_rcpp(t(bhat), t(sbhat), t(s_alpha), matrix(0,0,0), 
                              private$residual_correlation,
                              matrix(0,0,0), matrix(0,0,0), 
                              private$.prior_variance$xUlist,
                              t(private$.mixture_posterior_weights),
                              is_common_cov, 4)
      } else {
        post = mashr:::calc_post_precision_rcpp(t(bhat), t(sbhat), t(s_alpha), matrix(0,0,0), 
                              private$residual_correlation,
                              matrix(0,0,0), matrix(0,0,0), 
                              private$precomputed_cov_matrices$Vinv,
                              private$precomputed_cov_matrices$U0,
                              t(private$.mixture_posterior_weights), 
                              is_common_cov, 4)
      }
      private$.posterior_b1 = post$post_mean
      # Format post_cov for degenerated case with R = 1
      # (no need for it)
      #if (ncol(private$.posterior_b1) == 1) {
      #  post$post_cov = array(post$post_cov, c(1, 1, private$J))
      #}
      private$.posterior_b2 = post$post_cov + simplify2array(lapply(1:nrow(post$post_mean), function(i) tcrossprod(post$post_mean[i,])))
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
      # FIXME: NA situation
      private$.lbf[which(is.na(private$.lbf))] = 0
    },
    compute_loglik_null = function(d) {}
  ),
  private = list(
    residual_correlation = NULL,
    precomputed_cov_matrices = NULL,
    alpha = NULL,
    .mixture_posterior_weights = NULL,
    .lfsr = NULL,
    .residual_variance_inv = NULL
  ),
  active = list(
    residual_variance_inv = function(v) {
      if (missing(v)) private$.residual_variance_inv
      else private$.residual_variance_inv = v
    }
  )
)

#' @title MASH initializer object
#' @importFrom R6 R6Class
#' @importFrom mashr expand_cov
#' @keywords internal
MashInitializer <- R6Class("MashInitializer",
  public = list(
      initialize = function(Ulist, grid, prior_weights = NULL, null_weight = 0, alpha = 1, weights_tol = 1E-10, top_mixtures = 20, xUlist = NULL, include_conditions = NULL) {
        all_zeros = vector()
        if (is.null(xUlist)) {
          # FIXME: mashr::: namespace
          for (l in 1:length(Ulist)) {
              if (all(Ulist[[l]] == 0))
              stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
          }
          if (any(grid<=0)) stop("grid values should be greater than zero")
          private$U = list(pi = weights, Ulist = Ulist, grid = grid, usepointmass = TRUE)
          if (!is.null(include_conditions)) {
            for (l in 1:length(Ulist)) {
              Ulist[[l]] = Ulist[[l]][include_conditions, include_conditions]
              all_zeros[l] = all(Ulist[[l]] == 0)
            }
          }
          xUlist = expand_cov(Ulist, grid, usepointmass=TRUE)
        }
        plen = length(xUlist) - 1
        if (is.null(prior_weights)) prior_weights = rep(1/plen, plen)
        if (length(prior_weights) != plen)
          stop(paste("Invalid prior_weights setting: expect length", plen, "but input is of length", length(prior_weights)))
        weights = c(null_weight, prior_weights)
        weights = weights / sum(weights)
        # Filter by weights lower bound
        # Have to keep the first null component
        which.comp = which(weights[-1] > weights_tol)
        which.comp = c(1, which.comp + 1)
        filtered_weights = weights[which.comp]
        xUlist = xUlist[which.comp]
        # There are all zero priors, after some conditions are removed
        # we will have to adjust the prior weights based on it
        # This is a not very efficient yet safe and clear way to do it
        if (length(which(all_zeros))>0) {
          non_zeros = which(sapply(1:length(xUlist), function(l) !all(xUlist[[l]] == 0)))
          xUlist = xUlist[c(1, non_zeros[-1])]
          filtered_weights = filtered_weights[c(1, non_zeros[-1])]
        }
        # Filter for top weights: we only keep top weights
        if (top_mixtures > 0 && top_mixtures < length(filtered_weights)) {
          which.comp = head(sort(filtered_weights[-1], index.return=T, decreasing=T)$ix, top_mixtures)
          which.comp = c(1, which.comp + 1)
          filtered_weights = weights[which.comp]
          xUlist = xUlist[which.comp]
        }
        # Check on xUlist
        u_rows = vector()
        for (i in 1:length(xUlist)) {
          mashr:::check_covmat_basics(xUlist[[i]])
          u_rows[i] = nrow(xUlist[[i]])
        }
        if (length(unique(u_rows)) > 1) stop("Ulist contains matrices of different dimensions.")
        private$xU = list(pi = filtered_weights / sum(filtered_weights), xUlist = xUlist)
        private$a = alpha
      },
    precompute_cov_matrices = function(d, residual_covariance, algorithm = c('cpp', 'R')) {
      # computes constants (SVS + U)^{-1} and (SVS)^{-1} for posterior
      # and sigma_rooti for likelihooods
      # output of this function will provide input to `mashr`'s
      # functions calc_lik_common_rcpp() and
      # calc_post_precision_rcpp()
      # The input should be sbhat data matrix
      # d[j,] can be different for different conditions due to missing Y data
      # FIXME: did not use alpha information
      V = cov2cor(residual_covariance)
      sigma2 = diag(residual_covariance)
      if (d$Y_has_missing()) {
        res = get_sumstats_missing_data(d$X, d$Y, sigma2, V, private$a)
        svs = res$svs
        sbhat0 = res$sbhat0
        sbhat = sbhat0 ^ (1 - private$a)
        common_sbhat = is_mat_common(sbhat)
      } else {
        sbhat0 = sqrt(do.call(rbind, lapply(1:length(d$d), function(j) sigma2 / d$d[j])))
        sbhat0[which(is.nan(sbhat0) | is.infinite(sbhat0))] = 1E6
        sbhat = sbhat0 ^ (1 - private$a)
        common_sbhat = is_mat_common(sbhat)
        if (common_sbhat) svs = sbhat[1,] * t(V * sbhat[1,]) # faster than diag(s) %*% V %*% diag(s)
        else svs = lapply(1:nrow(sbhat), function(j) sbhat[j,] * t(V * sbhat[j,]))
      }
      # the `if` condition is used due to computational reasons: we can save RxRxP matrices but not RxRxPxJ
      # FIXME: compute this in parallel in the future
      algorithm = match.arg(algorithm)
      if (is.matrix(svs)) {
        # sigma_rooti is R * R * P
        # this is in preparation for some constants used in dmvnrom() for likelihood calculations
        sigma_rooti = list()
        for (i in 1:length(private$xU$xUlist)) {
          if (algorithm == 'R') sigma_rooti[[i]] = t(backsolve(muffled_chol(svs + private$xU$xUlist[[i]]), diag(nrow(svs))))
          else sigma_rooti[[i]] = mashr:::calc_rooti_rcpp(svs + private$xU$xUlist[[i]])$data
        }
        # this is in prepartion for some constants used in posterior calculation
        Vinv = list()
        Vinv[[1]] = solve(svs)
        U0 = list()
        for (i in 1:length(private$xU$xUlist)) U0[[i]] = private$xU$xUlist[[i]] %*% solve(Vinv[[1]] %*% private$xU$xUlist[[i]] + diag(nrow(private$xU$xUlist[[i]])))
      } else {
          # have to do this for every effect
          # sigma_rooti and U0 will be R * R * (J * P)
          # and Vinv will be a J list, not a matrix
          # this is in preparation for some constants used in dmvnrom() for likelihood calculations
          sigma_rooti = list()
          k = 1
          for (j in 1:length(svs)) {
            for (i in 1:length(private$xU$xUlist)) {
              if (algorithm == 'R') {
                sigma_rooti[[k]] = t(backsolve(muffled_chol(svs[[j]] + private$xU$xUlist[[i]]), diag(nrow(svs[[j]]))))
              } else {
                sigma_rooti[[k]] = mashr:::calc_rooti_rcpp(svs[[j]] + private$xU$xUlist[[i]])$data
              }
              k = k + 1
            }
          }
          Vinv = lapply(1:length(svs), function(i) solve(svs[[i]]))
          U0 = list()
          k = 1
          for (j in 1:length(svs)) {
            for (i in 1:length(private$xU$xUlist)) {
                U0[[k]] = private$xU$xUlist[[i]] %*% solve(Vinv[[j]] %*% private$xU$xUlist[[i]] + diag(nrow(private$xU$xUlist[[i]])))
                k = k + 1
            }
          }
      }
      private$inv_mats = list(Vinv = simplify2array(Vinv), U0 = simplify2array(U0),
                              sigma_rooti = simplify2array(sigma_rooti), sbhat = sbhat0, common_sbhat = common_sbhat)
    }
  ),
  private = list(
      U = NULL,
      xU = NULL,
      a = NULL,
      inv_mats = NULL
  ),
  active = list(
      n_condition = function(v) nrow(private$xU$xUlist[[1]]),
      prior_covariance = function() private$xU,
      mash_prior = function() private$U,
      precomputed = function() private$inv_mats,
      alpha = function(value) {
        if (missing(value)) return(private$a)
        else private$a = value
      }
  )
)