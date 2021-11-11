#' @title Compute condition specific posterior inclusion probability.
#' @description This is only relevant when canonical priors are used
#' @param m M&M model
#' @param prior_obj prior mixture object
#' @return P by R matrix of PIP per condition
#' @keywords internal
mvsusie_get_pip_per_condition = function(m, prior_obj) {
  condition_pip = mvsusie_get_alpha_per_condition(m, prior_obj)
  return(do.call(cbind, lapply(1:dim(condition_pip)[3], function(r) apply(condition_pip[,,r], 2, function(x) 1-prod(1-x)))))
}

#' @title Compute condition specific posterior inclusion probability per effect
#' @keywords internal
mvsusie_get_alpha_per_condition = function(m, prior_obj) {
  condition_indicator = do.call(rbind, lapply(1:length(prior_obj$prior_variance$xUlist), function(i) as.integer(diag(prior_obj$prior_variance$xUlist[[i]]) != 0)))
  condition_pip = array(0, dim=dim(m$b1))
  for (r in 1:dim(condition_pip)[3]) {
    for (p in 1:length(condition_indicator[,r])) {
        condition_pip[,,r] = condition_pip[,,r] + m$mixture_weights[,,p] * condition_indicator[p,r]
    }
    condition_pip[,,r] = condition_pip[,,r] * m$alpha
  }
  return(condition_pip)
}

#' @title Local false sign rate (lfsr) for single effects
#' @details This computes the lfsr of single effects for each condition.
#' @param alpha L by P matrix
#' @param clfsr L by P by R conditonal lfsr
#' @return a L by R matrix of lfsr
#' @export
mvsusie_single_effect_lfsr = function(clfsr, alpha) {
  if(!is.array(clfsr) && is.na(clfsr)){
    return(NA)
  }else{
    return(do.call(cbind, lapply(1:dim(clfsr)[3], function(r){
      clfsrr = clfsr[,,r]
      if (is.null(nrow(clfsrr))) clfsrr = matrix(clfsrr, 1, length(clfsrr))
      pmax(0, rowSums(alpha * clfsrr))
    })))
  }
}

#' @title Local false sign rate (lfsr) for variables
#' @details This computes the lfsr of variables for each condition.
#' @param alpha L by P matrix
#' @param clfsr L by P by R conditonal lfsr
#' @param weighted TRUE to weight lfsr by PIP; FALSE otherwise.
#' @return a P by R matrix of lfsr
#' @export
mvsusie_get_lfsr = function(clfsr, alpha, weighted = TRUE) {
  if(!is.array(clfsr) && is.na(clfsr)){
    return(NA)
  }else{
    if (weighted) alpha = alpha
    else alpha = matrix(1, nrow(alpha), ncol(alpha))
    return(do.call(cbind, lapply(1:dim(clfsr)[3], function(r){
      true_sign_mat = alpha * (1 - clfsr[,,r])
      pmax(1e-20, 1 - apply(true_sign_mat, 2, max))
    })))
  }
}

# SuSiE model extractor
# 
#' @importFrom abind abind
report_susie_model = function(d, m, estimate_prior_variance = TRUE) {
    if (length(dim(m$posterior_b1[[1]])) < 2) {
      # univariate case
      b1 = t(do.call(cbind, m$posterior_b1))
      b2 = t(do.call(cbind, m$posterior_b2))
      b = colSums(b1)
    } else {
      b1 = aperm(abind::abind(m$posterior_b1,along=3), c(3,1,2))
      b2 = aperm(abind::abind(m$posterior_b2,along=3), c(3,1,2))
      if (dim(b1)[1] == 1) {
        # only one effect specified or left
        b = do.call(cbind, lapply(1:dim(b1)[3], function(i) b1[,,i]))
      } else {
        # multiple effects
        b = do.call(cbind, lapply(1:dim(b1)[3], function(i) colSums(b1[,,i])))
      }
      if (dim(b)[2] == 1) {
        b1 = b1[,,1]
        b2 = b2[,,1]
        b = as.vector(b)
      }
    }
    if (is.null(m$mixture_posterior_weights[[1]])) mixture_weights = NA
    else {
      if (length(dim(m$mixture_posterior_weights[[1]])) < 2) mixture_weights = t(do.call(cbind, m$mixture_posterior_weights))
      else mixture_weights = aperm(abind::abind(m$mixture_posterior_weights,along=3), c(3,1,2))
    }
    if (is.null(m$clfsr[[1]])) clfsr = NA
    else {
      if (length(dim(m$clfsr[[1]])) < 2) clfsr = t(do.call(cbind, m$clfsr))
      else clfsr = aperm(abind::abind(m$clfsr,along=3), c(3,1,2))
    }
    s = list(
        alpha = t(m$pip),
        b1 = b1,
        b2 = b2,
        KL = m$kl,
        lbf = m$lbf,
        V = m$prior_variance,
        sigma2 = d$residual_variance,
        elbo = m$get_objective(dump=TRUE),
        niter = m$niter,
        convergence = m$convergence,
        coef = d$rescale_coef(b),
        mixture_weights = mixture_weights,
        conditional_lfsr = clfsr,
        lfsr = mvsusie_get_lfsr(clfsr, t(m$pip)),
        single_effect_lfsr = mvsusie_single_effect_lfsr(clfsr, t(m$pip))
        )
    if (!is.null(m$pip_history)) s$alpha_history = m$pip_history
    if (!is.null(m$lbf_history)) s$lbf_history = m$lbf_history
    if (!is.null(m$prior_history)) s$prior_history = m$prior_history
    # FIXME: unit test for scaling issue for the fitted
    if(inherits(d, "RSSData")){
      s$fitted = d$XtX %*% b
    }else{
      s$fitted = d$compute_Xb(b)
    }
    if (is.null(dim(s$coef))) s$intercept = s$coef[1]
    else s$intercept = s$coef[1,]
    if (estimate_prior_variance) s$V = m$prior_variance
    class(s) = 'susie'
    return(s)
}
