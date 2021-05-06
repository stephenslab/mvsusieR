#' @title MASH multiple regression object
#' @importFrom R6 R6Class
#' @importFrom ashr compute_lfsr
#' @keywords internal
MashRegression <- R6Class("MashRegression",
  inherit = BayesianSimpleRegression,
  public = list(
    initialize = function(J, mash_initializer) {
      private$J = J
      private$.prior_variance = mash_initializer$prior_variance
      private$.prior_variance$xUlist = matlist2array(private$.prior_variance$xUlist)
      private$precomputed_cov_matrices = mash_initializer$precomputed
      if (is.null(private$.prior_variance$xUlist_inv))
        private$.prior_variance$xUlist_inv = 0
      private$.posterior_b1 = matrix(0, J, mash_initializer$n_condition)
      private$prior_variance_scale = 1
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE, save_var = FALSE, estimate_prior_variance_method = NULL, check_null_threshold = 0) {
      # When prior changes (private$prior_variance_scale != 1),
      # we can no longer use precomputed quantities
      # because the precomputed quantities will be wrong in scale.
      private$residual_correlation = d$residual_correlation
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      # bhat is J by R
      bhat = d$get_coef(use_residual)
      sbhat = d$sbhat
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sbhat
      }
      if (!is.null(estimate_prior_variance_method) && estimate_prior_variance_method != "EM") {
        if (estimate_prior_variance_method != 'simple')
          stop(paste("Estimate prior method", estimate_prior_variance_method, "is not available for MashRegression."))
        private$prior_variance_scale = private$estimate_prior_variance(bhat,sbhat,prior_weights,method=estimate_prior_variance_method,check_null_threshold=check_null_threshold)
      }
      if (!is.null(estimate_prior_variance_method) && estimate_prior_variance_method != 'simple' && !is.null(private$precomputed_cov_matrices$U0)) {
        # Cannot use precomputed quantities if prior variance scalar is being estimated
        # we set it to null so it will not be used later
        private$precomputed_cov_matrices$U0 = NULL
      }
      # Fit MASH model
      # 1. compute log-likelihood matrix given current estimates
      llik = private$compute_loglik_mat(private$prior_variance_scale, bhat, sbhat, d$svs, d$is_common_cov)
      # 2. lbf
      lbf_obj = private$compute_lbf(llik)
      private$.lbf = lbf_obj$lbf
      private$.loglik_null = lbf_obj$loglik_null
      # 3. compute posterior weights
      private$.mixture_posterior_weights = private$compute_mixture_posterior_weights(private$.prior_variance$pi, llik)
      if (!is.null(estimate_prior_variance_method) && estimate_prior_variance_method == 'EM') {
        variable_posterior_weights = private$compute_variable_posterior_weights(prior_weights, llik)
        private$cache = list(b=bhat, s=sbhat, update_scale=T)
      } else {
        variable_posterior_weights = matrix(0,0,0)
      }
      # 4. posterior
      ## FIXME: we do not need to compute second moment unless:
      # 1. need ELBO to check for convergence
      # 2. need ELBO to estimate residual variance
      # 3. need to update prior via EM
      # but let's compute it here anyways
      post = private$compute_posterior(bhat, sbhat, d$svs_inv, d$is_common_cov, private$.mixture_posterior_weights, variable_posterior_weights)
      private$.posterior_b1 = post$post_mean
      private$.posterior_b2 = post$post_cov + matlist2array(lapply(1:nrow(post$post_mean), function(i) tcrossprod(post$post_mean[i,])))
      if (save_var) private$.posterior_variance = post$post_cov
      # flatten posterior_b2 for degenerated case with R = 1
      if (ncol(private$.posterior_b1) == 1) {
        private$.posterior_b2 = as.matrix(apply(private$.posterior_b2,3,diag))
        if (!is.null(private$.posterior_variance)) private$.posterior_variance = as.matrix(apply(private$.posterior_variance,3,diag))
      }
      # 5. lfsr
      private$.lfsr = compute_lfsr(post$post_neg, post$post_zero)
      # 6. estimate prior via EM
      if (!is.null(estimate_prior_variance_method) && estimate_prior_variance_method == 'EM') {
        private$cache$SER_posterior_mixture_weights = private$get_SER_posterior_mixture_weights(llik, prior_weights, private$.prior_variance$pi)
        private$cache$mixture_prior_variance_scale = post$prior_scale_em_update
      }
      # 7. clean up workspace
      rm(post)
      rm(llik)
    }
  ),
  active = list(
    mixture_posterior_weights = function() private$.mixture_posterior_weights,
    lfsr = function() private$.lfsr,
    prior_variance = function(v) {
      if (missing(v)) private$prior_variance_scale
      else private$prior_variance_scale = v
    }
  ),
  private = list(
    precomputed_cov_matrices = NULL,
    prior_variance_scale = NULL,
    .mixture_posterior_weights = NULL,
    .lfsr = NULL,
    residual_correlation = NULL,
    compute_loglik_mat = function(scalar, bhat, sbhat, svs, is_common_cov) {
      if (is.null(private$precomputed_cov_matrices$sigma_rooti) || (scalar != 1 && scalar != 0)) {
        llik = mashr:::calc_lik_rcpp(t(bhat),
                                    # t(sbhat) and d$residual_correlation can both be empty (matrix(0,0,0)) if SVS is provided
                                    t(sbhat),
                                    private$residual_correlation,
                                    matrix(0,0,0),
                                    private$get_scaled_prior(scalar),
                                    # should be matlist2array(d$svs), if t(sbhat) and d$residual_correlation are not empty
                                    matlist2array(svs),
                                    TRUE,
                                    is_common_cov,
                                    private$n_thread)$data
      } else {
        # Here private$prior_variance_scale is either 0 or 1.
        # This line below assumes it is 1; will adjust it after for case of 0.
        llik = mashr:::calc_lik_precomputed_rcpp(t(bhat),
                                         private$precomputed_cov_matrices$sigma_rooti,
                                         TRUE,
                                         is_common_cov,
                                         private$n_thread)$data
        if (scalar == 0) {
          # The precomputed sigma_rooti is not correct
          # but the first column of llik is llik under the null anyways
          # that corresponds to scalar == 0
          # so we can simply set all columns of llik to the first column
          llik = replicate(ncol(llik), llik[,1])
        }
      }
      # give a warning if any columns have -Inf likelihoods.
      rows = which(apply(llik,2,function (x) any(is.infinite(x))))
      if (length(rows) > 0)
        warning(paste("Some mixture components result in non-finite likelihoods,",
                          "either\n","due to numerical underflow/overflow,",
                          "or due to invalid covariance matrices",
                          paste(rows,collapse=", "), "\n"))
      return(llik)
    },
    compute_posterior = function(bhat, sbhat, svs_inv, is_common_cov, mixture_posterior_weights, variable_posterior_weights) {
      if (is.null(private$precomputed_cov_matrices$U0) || (private$prior_variance_scale != 1 && private$prior_variance_scale != 0)) {
        post = mashr:::calc_sermix_rcpp(t(bhat),
                              # sbhat is not needed (can safely be replaced by matrix(0,0,0)) IF Vinv is provided
                              t(sbhat),
                              # residual correlation is not needed (can safely be replaced by matrix(0,0,0)) IF Vinv is provided
                              private$residual_correlation,
                              matlist2array(svs_inv),
                              private$get_scaled_prior(private$prior_variance_scale),
                              # because we define the scalar with respect to the original prior
                              # the inverse should always be the original.
                              private$.prior_variance$xUlist_inv,
                              0,
                              t(mixture_posterior_weights),
                              t(variable_posterior_weights),
                              is_common_cov,
                              private$n_thread)
      } else {
        # Use precomputed quantities
        # here private$prior_variance_scale is either 0 or 1
        post = mashr:::calc_sermix_rcpp(t(bhat),
                              # No need for sbhat and residual correlation when Vinv is precomputed
                              # So we just put in an empty matrix for them (matrix(0,0,0)).
                              matrix(0,0,0), matrix(0,0,0),
                              matlist2array(svs_inv),
                              private$get_scaled_prior(private$prior_variance_scale),
                              private$.prior_variance$xUlist_inv,
                              private$precomputed_cov_matrices$U0 * private$prior_variance_scale,
                              t(mixture_posterior_weights),
                              matrix(0,0,0),
                              is_common_cov,
                              private$n_thread)
      }
      return(post)
    },
    compute_mixture_posterior_weights = function(prior_mixture_weights, llik) {
      lfactors = apply(llik,1,max)
      d = t(prior_mixture_weights * t(exp(llik-lfactors)))
      return(d/rowSums(d))
    },
    compute_variable_posterior_weights = function(prior_variable_weights, llik) {
      lbf = t(llik - llik[,1])
      return(t(private$compute_mixture_posterior_weights(prior_variable_weights, lbf)))
    },
    compute_lbf = function(llik, s = NULL) {
      # get relative loglik
      lfactors = apply(llik,1,max)
      llik = list(loglik_matrix=llik-lfactors, lfactors=lfactors)
      # using mashr functions have to ensure input s_alpha parameter has valid log and rowSums
      if (is.null(s) || (is.matrix(s) && nrow(s) == 0)) s = matrix(1,1,1)
      loglik_null = mashr:::compute_null_loglik_from_matrix(llik, s)
      loglik_alt = mashr:::compute_alt_loglik_from_matrix_and_pi(private$.prior_variance$pi, llik, s)
      lbf = loglik_alt - loglik_null
      if (!is.null(ncol(lbf)) && ncol(lbf) == 1)
        lbf = as.vector(lbf)
      # Inf - Inf above can cause NaN
      lbf[which(is.na(lbf))] = 0
      return(list(lbf=lbf, loglik_null=loglik_null))
    },
    loglik = function(V,B,S,prior_weights) {
      llik = private$compute_loglik_mat(V,B,S)
      return(compute_softmax(private$compute_lbf(llik)$lbf, prior_weights)$log_sum)
    },
    get_SER_posterior_mixture_weights = function(llik, prior_weights, prior_mixture_weights) {
      # This function computes p(\gamma_p) in estimate_prior_variance_em() function
      lbf = llik - llik[,1]
      ser_lbf = apply(lbf, 2, function(x) compute_softmax(x, prior_weights)$log_sum)
      return(compute_softmax(ser_lbf, prior_mixture_weights)$weights)
    },
    estimate_prior_variance_em = function(pip) {
      # The EM update is
      # \sigma_0^2 = \sum_{p=1}^P p(\gamma_p) \mathrm{tr}(U_p^{-1} E[bb^T \,|\, \gamma_p])/r
      # where E[bb^T \,|\, \gamma_p] = \sum_j \alpha_{p,j} * mu2_mat_{p,j}
      # The trace(.) / r part has already been computed in function calc_sermix_rcpp()
      # the output is saved as private$cache$mixture_prior_variance_scale
      # The (\gamma_p) part has already been computed in function get_SER_posterior_mixture_weights()
      # the output is saved as private$cache$SER_posterior_mixture_weights
      # Here PIP is not used. The notion of PIP here has been reflected in
      # variable_posterior_weights an input to calc_sermix_rcpp()
      # this PIP is for per mixture component.
      V = sum(private$cache$SER_posterior_mixture_weights * private$cache$mixture_prior_variance_scale)/
        sum(private$cache$SER_posterior_mixture_weights * attr(private$.prior_variance$xUlist_inv, 'rank'))
      return(V)
    },
    estimate_prior_variance_simple = function() 1,
    get_scaled_prior = function(scalar) {
      # xUlist here is a 3D array
      if (scalar != 1) {
        return(private$.prior_variance$xUlist * scalar)
      } else {
        return(private$.prior_variance$xUlist)
      }
    }
  ),
)

#' @title MASH initializer object
#' @importFrom R6 R6Class
#' @keywords internal
MashInitializer <- R6Class("MashInitializer",
  public = list(
      initialize = function(Ulist, grid, prior_weights = NULL, null_weight = 0, weights_tol = 1E-10, null_tol = 5E-7, top_mixtures = 20, xUlist = NULL, include_conditions = NULL) {
        all_zeros = vector()
        if (is.null(xUlist)) {
          if (is.null(Ulist)) stop("Either xUlist or Ulist have to be non-null")
          for (l in 1:length(Ulist)) {
              if (all(abs(Ulist[[l]])< null_tol)) stop(paste("Prior covariance", l , "is zero matrix. This is not allowed."))
          }
          if (any(grid<=0)) stop("grid values should be greater than zero") 
          xUlist = mashr:::expand_cov(Ulist, grid, usepointmass=TRUE)
        } else {
          if (!all(xUlist[[1]] == 0)) xUlist = c(list(matrix(0, nrow(xUlist[[1]]), ncol(xUlist[[1]]))), xUlist)
        }
        if (!is.null(include_conditions)) {
            for (l in 1:length(xUlist)) {
              xUlist[[l]] = xUlist[[l]][include_conditions, include_conditions]
              if (l > 1) all_zeros[l-1] = all(abs(xUlist[[l]])< null_tol)
            }
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
        u_rows = vector(length=length(xUlist))
        for (i in 1:length(xUlist)) {
          mashr:::check_covmat_basics(xUlist[[i]])
          u_rows[i] = nrow(xUlist[[i]])
          if (!mashr:::issemidef(xUlist[[i]]))
            stop(paste0("The prior matrices ", i, " should be positive semi-definite."))
        }
        if (length(unique(u_rows)) > 1) stop("Ulist contains matrices of different dimensions.")
        prior_weights = prior_weights / sum(prior_weights)
        private$xU = list(pi = c(null_weight, prior_weights * (1 - null_weight)), xUlist = xUlist)
      },
    compute_prior_inv = function() {
      # compute pseudo inverse for prior matrices and divided by its rank
      # this is relevant to the EM update of prior variance scalar
      K = length(private$xU$xUlist)
      Uinv = vector("list", length = K)
      Urank = numeric(K)
      for(i in 1:K){
        uinv = pseudo_inverse(private$xU$xUlist[[i]])
        Uinv[[i]] = uinv$inv
        Urank[i] = uinv$rank
      }
      private$xU$xUlist_inv = matlist2array(Uinv)
      attr(private$xU$xUlist_inv, 'rank') = Urank
    },
    precompute_cov_matrices = function(d, algorithm = c('R', 'cpp')) {
      # computes constants (SVS + U)^{-1} and (SVS)^{-1} for posterior
      # and sigma_rooti for likelihooods
      # output of this function will provide input to `mashr`'s
      # functions calc_lik_common_rcpp() and
      # calc_post_precision_rcpp()
      # The input should be sbhat data matrix
      # d[j,] can be different for different conditions due to missing Y data
      # the `if` condition is used due to computational reasons: we can save RxRxP matrices but not RxRxPxJ
      # FIXME: compute this in parallel in the future
      algorithm = match.arg(algorithm)

      if (d$is_common_cov) {
        K = length(private$xU$xUlist)
        # sigma_rooti is R * R * P
        # this is in preparation for some constants used in dmvnrom() for likelihood calculations
        sigma_rooti = vector("list", length = K)
        # this is in prepartion for some constants used in posterior calculation
        U0 = vector("list", length = K)
        for (i in 1:K) {
          if (algorithm == 'R') sigma_rooti[[i]] = invert_chol_tri(d$svs[[1]] + private$xU$xUlist[[i]])$inv
          else sigma_rooti[[i]] = mashr:::inv_chol_tri_rcpp(d$svs[[1]] + private$xU$xUlist[[i]])$data
          U0[[i]] = private$xU$xUlist[[i]] %*% solve(d$svs_inv[[1]] %*% private$xU$xUlist[[i]] + diag(nrow(private$xU$xUlist[[i]])))
        }
      } else {
        # have to do this for every effect
        # sigma_rooti and U0 will be R * R * (J * P)
        # and Vinv will be a J list, not a matrix
        # this is in preparation for some constants used in dmvnrom() for likelihood calculations
        K = length(private$xU$xUlist) * d$n_effect
        sigma_rooti = vector("list", length = K)
        U0 = vector("list", length = K)
        k = 1
        for (j in 1:d$n_effect) {
          for (i in 1:length(private$xU$xUlist)) {
            if (algorithm == 'R') {
              sigma_rooti[[k]] = invert_chol_tri(d$svs[[j]] + private$xU$xUlist[[i]])$inv
            } else {
              sigma_rooti[[k]] = mashr:::inv_chol_tri_rcpp(d$svs[[j]] + private$xU$xUlist[[i]])$data
            }
            U0[[k]] = private$xU$xUlist[[i]] %*% solve(d$svs_inv[[j]] %*% private$xU$xUlist[[i]] + diag(nrow(private$xU$xUlist[[i]])))
            k = k + 1
          }
        }
      }
      private$inv_mats = list(U0 = matlist2array(U0),
                              sigma_rooti = matlist2array(sigma_rooti))
    },
    remove_precomputed = function() private$inv_mats = NULL,
    scale_prior_variance = function(sigma) {
      private$xU$xUlist = lapply(1:length(private$xU$xUlist), function(i) scale_covariance(private$xU$xUlist[[i]], sigma))
    }
  ),
  active = list(
      n_condition = function() nrow(private$xU$xUlist[[1]]),
      n_component = function() length(private$xU$xUlist),
      prior_variance = function() private$xU,
      precomputed = function() private$inv_mats
  ),
  private = list(
      U = NULL,
      xU = NULL,
      inv_mats = NULL
  )
)