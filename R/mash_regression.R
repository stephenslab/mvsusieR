#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @importFrom ashr compute_lfsr
#' @keywords internal
MashRegression <- R6Class("MashRegression",
  inherit = BayesianSimpleRegression,
  public = list(
    initialize = function(J, residual_variance, mash_initializer) {
      private$J = J
      private$.prior_variance = mash_initializer$prior_variance
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
      private$prior_variance_scale = 1
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE, estimate_prior_variance_method = NULL) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      # bhat is J by R
      bhat = XtY / d$X2_sum
      bhat[which(is.nan(bhat))] = 0
      if (!is.null(private$precomputed_cov_matrices)) {
        sbhat = private$precomputed_cov_matrices$sbhat
      } else {
        # sbhat is R by R
        # for non-missing Y d$X2_sum is a J vector
        # for missing Y it is a J by R matrix
        sigma2 = diag(private$.residual_variance)
        if (d$Y_has_missing) sbhat = sqrt(do.call(rbind, lapply(1:nrow(d$X2_sum), function(j) sigma2 / d$X2_sum[j,])))
        else sbhat = sqrt(do.call(rbind, lapply(1:length(d$X2_sum), function(j) sigma2 / d$X2_sum[j])))
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
      # Fit MASH model
      if (!is.null(private$precomputed_cov_matrices)) is_common_cov = private$precomputed_cov_matrices$common_sbhat
      else is_common_cov = is_mat_common(sbhat)
      # 1.1 compute log-likelihood matrix given current estimates
      if (is.null(private$precomputed_cov_matrices) || ncol(bhat) == 1) {
        llik = mashr:::calc_lik_rcpp(t(bhat), t(sbhat), private$residual_correlation,
                                         matrix(0,0,0),
                                         private$.prior_variance$xUlist,
                                         TRUE,
                                         is_common_cov)$data
      } else {
        llik = mashr:::calc_lik_rooti_rcpp(t(bhat),
                                         private$precomputed_cov_matrices$sigma_rooti,
                                         TRUE,
                                         is_common_cov)$data
      }

      # 1.2 give a warning if any columns have -Inf likelihoods.
      rows <- which(apply(llik,2,function (x) any(is.infinite(x))))
      if (length(rows) > 0)
        warning(paste("Some mixture components result in non-finite likelihoods,",
                          "either\n","due to numerical underflow/overflow,",
                          "or due to invalid covariance matrices",
                          paste(rows,collapse=", "), "\n"))
      # 1.3 get relative loglik
      lfactors = apply(llik,1,max)
      llik = list(loglik_matrix=llik-lfactors, lfactors=lfactors)
      # 2. lbf
      lbf_obj = private$compute_lbf(llik, s_alpha)
      private$.lbf = lbf_obj$lbf
      private$.loglik_null = lbf_obj$loglik_null
      if (!is.null(estimate_prior_variance_method)) {
        if (estimate_prior_variance_method != 'simple')
          stop(paste("Estimate prior method", estimate_prior_variance_method, "is not available for MashRegression."))
        if (is.null(prior_weights))
          prior_weights = rep(1/private$J, private$J)
        # the commented out line is consistent with comparing loglik for V=0 vs current V
        # but is more agressive than checking if all `private$.lbf` are smaller than zero.
        #if (compute_weighted_sum(private$.lbf, prior_weights)$log_sum <= 0) {
        if (all(private$.lbf <= 0)) {
          private$prior_variance_scale = 0
          private$.lbf = private$.lbf * 0
        }
      }
      # 3. compute posterior weights
      private$.mixture_posterior_weights = mashr:::compute_posterior_weights(private$.prior_variance$pi, exp(llik$loglik_matrix))
      # 4. posterior
      ## FIXME: we might not need to compute second moment at all if we do not need to estimate residual variance
      ## we can get away with checking for convergence by PIP not by ELBO
      ## but let's set report_type = 4 and compute posterior covariance for now
      if (is.null(private$precomputed_cov_matrices) || ncol(bhat) == 1 || private$prior_variance_scale != 1) {
        if (private$prior_variance_scale != 1)
          xUlist = private$.prior_variance$xUlist * private$prior_variance_scale
        else
          xUlist = private$.prior_variance$xUlist
        post = mashr:::calc_post_rcpp(t(bhat), t(sbhat), t(s_alpha), matrix(0,0,0),
                              private$residual_correlation,
                              matrix(0,0,0), matrix(0,0,0),
                              xUlist,
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
      # 5. lfsr
      private$.lfsr = compute_lfsr(post$post_neg, post$post_zero)
    }
  ),
  active = list(
    mixture_posterior_weights = function() private$.mixture_posterior_weights,
    lfsr = function() private$.lfsr,
    residual_variance_inv = function() private$.residual_variance_inv,
    prior_variance = function() NA,
    residual_variance = function(v) {
      if (missing(v)) private$.residual_variance
      else warning(paste("residual variance cannot be changed for object", "MashRegression"))
    }
  ),
  private = list(
    residual_correlation = NULL,
    precomputed_cov_matrices = NULL,
    alpha = NULL,
    prior_variance_scale = NULL,
    .mixture_posterior_weights = NULL,
    .lfsr = NULL,
    .residual_variance_inv = NULL,
    compute_lbf = function(llik, s) {
      # using mashr functions have to ensure s_alpha has valid log and rowSums
      if (nrow(s) == 0) s = matrix(1,1,1)
        loglik_null = mashr:::compute_null_loglik_from_matrix(llik, s)
        loglik_alt = mashr:::compute_alt_loglik_from_matrix_and_pi(private$.prior_variance$pi, llik, s)
        lbf = loglik_alt - loglik_null
        if (!is.null(ncol(lbf)) && ncol(lbf) == 1)
          lbf = as.vector(lbf)
        # Inf - Inf above can cause NaN
        lbf[which(is.na(lbf))] = 0
        return(list(lbf=lbf, loglik_null=loglik_null))
    }
  ),
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
          for (l in 1:length(Ulist)) {
              if (all(Ulist[[l]] == 0))
              stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
          }
          if (any(grid<=0)) stop("grid values should be greater than zero")
          if (!is.null(include_conditions)) {
            for (l in 1:length(Ulist)) {
              Ulist[[l]] = Ulist[[l]][include_conditions, include_conditions]
              all_zeros[l] = all(Ulist[[l]] == 0)
            }
          }
          xUlist = expand_cov(Ulist, grid, usepointmass=TRUE)
        } else {
          if (!all(xUlist[[1]] == 0)) xUlist = c(list(matrix(0, nrow(xUlist[[1]]), ncol(xUlist[[1]]))), xUlist)
        }
        plen = length(xUlist) - 1
        if (is.null(prior_weights)) prior_weights = rep(1/plen, plen)
        if (length(prior_weights) != plen)
          stop(paste("Invalid prior_weights setting: expect length", plen, "but input is of length", length(prior_weights)))
        # Filter by weights lower bound
        # Have to keep the first null component
        if (weights_tol > 0) {
          which.comp = which(prior_weights > weights_tol)
          prior_weights = prior_weights[which.comp]
          xUlist = xUlist[c(1, which.comp + 1)]
        }
        # There are all zero priors, after some conditions are removed
        # we will have to adjust the prior weights based on it
        # This is a not very efficient yet safe and clear way to do it
        if (length(which(all_zeros))>0) {
          # must exclude first xUlist which is always null here
          which.comp = which(sapply(2:length(xUlist), function(l) !all(xUlist[[l]] == 0)))
          prior_weights = prior_weights[which.comp]
          xUlist = xUlist[c(1, which.comp + 1)]
        }
        # Filter for top weights: we only keep top weights
        if (top_mixtures > 0 && top_mixtures < length(prior_weights)) {
          which.comp = head(sort(prior_weights, index.return=T, decreasing=T)$ix, top_mixtures)
          prior_weights = prior_weights[which.comp]
          xUlist = xUlist[c(1, which.comp + 1)]
        }
        # Check on xUlist
        u_rows = vector()
        for (i in 1:length(xUlist)) {
          mashr:::check_covmat_basics(xUlist[[i]])
          u_rows[i] = nrow(xUlist[[i]])
        }
        if (length(unique(u_rows)) > 1) stop("Ulist contains matrices of different dimensions.")
        prior_weights = prior_weights / sum(prior_weights)
        private$xU = list(pi = c(null_weight, prior_weights * (1 - null_weight)), xUlist = xUlist)
        private$a = alpha
      },
    precompute_cov_matrices = function(d, residual_covariance, algorithm = c('R', 'cpp')) {
      # computes constants (SVS + U)^{-1} and (SVS)^{-1} for posterior
      # and sigma_rooti for likelihooods
      # output of this function will provide input to `mashr`'s
      # functions calc_lik_common_rcpp() and
      # calc_post_precision_rcpp()
      # The input should be sbhat data matrix
      # d[j,] can be different for different conditions due to missing Y data
      # FIXME: did not use alpha information
      V = cov2cor(residual_covariance)
      res = d$get_sumstats(diag(residual_covariance), V, private$a)
      svs = res$svs
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
                              sigma_rooti = simplify2array(sigma_rooti),
                              sbhat = res$sbhat0,
                              common_sbhat = res$is_common_sbhat)
    },
    remove_precomputed = function() private$inv_mats = NULL
  ),
  active = list(
      n_condition = function() nrow(private$xU$xUlist[[1]]),
      prior_variance = function() private$xU,
      precomputed = function() private$inv_mats,
      alpha = function() private$a
  ),
  private = list(
      U = NULL,
      xU = NULL,
      a = NULL,
      inv_mats = NULL
  )
)