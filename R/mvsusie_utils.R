#' @title Local false sign rate (lfsr) for single effects
#' 
#' @details This function returns the lfsr for identifying nonzero
#'   single effects, separately for each condition.
#' 
#' @param alpha L x J matrix.
#' 
#' @param clfsr L x J x R conditonal lfsr.
#' 
#' @return L x R matrix of lfsr
#' 
#' @export
#' 
mvsusie_single_effect_lfsr = function (clfsr, alpha) {
  if (!is.array(clfsr) && is.na(clfsr))
    return(as.numeric(NA))
  else
    return(do.call(cbind,
                   lapply(1:dim(clfsr)[3],
                          function(r) {
                            clfsrr = clfsr[,,r]
                            if (is.null(nrow(clfsrr)))
                              clfsrr = matrix(clfsrr,1,length(clfsrr))
                            return(pmax(0,rowSums(alpha * clfsrr)))
                          })))
}

#' @title Local false sign rate (lfsr) for variables.
#' 
#' @details This function returns the lfsr for identifying nonzero
#'   effects for each condition.
#' 
#' @param alpha L x J matrix.
#' 
#' @param clfsr L x J x R conditonal lfsr.
#' 
#' @param weighted Set \code{weighted = TRUE} to weight lfsr by PIP;
#'   otherwise set \code{weighted = FALSE}.
#' 
#' @return J x R lfsr matrix.
#' 
#' @export
#' 
mvsusie_get_lfsr = function (clfsr, alpha, weighted = TRUE) {
  if (!is.array(clfsr) && is.na(clfsr))
    return(as.numeric(NA))
  else {
    if (weighted)
      alpha = alpha
    else
      alpha = matrix(1,nrow(alpha),ncol(alpha))
    return(do.call(cbind,
                   lapply(1:dim(clfsr)[3],
                          function (r) {
                            true_sign_mat = alpha * (1 - clfsr[,,r])
                            pmax(1e-20,1 - apply(true_sign_mat, 2, max))
                          })))
  }
}

# SuSiE model extractor
# 
#' @importFrom abind abind
report_susie_model = function (d, m, estimate_prior_variance = TRUE) {
  if (length(dim(m$posterior_b1[[1]])) < 2) {
      
    # univariate case
    b1 = t(do.call(cbind,m$posterior_b1))
    b2 = t(do.call(cbind,m$posterior_b2))
    b = colSums(b1)
  } else {
    b1 = aperm(abind(m$posterior_b1,along = 3),c(3,1,2)) # L x J x R
    b2 = aperm(abind(m$posterior_b2,along = 3),c(3,1,2)) # L x J x R
    if (dim(b1)[1] == 1)
        
      # only one effect specified or left
      b = do.call(cbind, lapply(1:dim(b1)[3],function(i) b1[,,i])) 
    else 
        
      # Multiple effects.
      b = do.call(cbind,lapply(1:dim(b1)[3],function(i) colSums(b1[,,i]))) # J x R

    if (dim(b)[2] == 1) {
      b1 = b1[,,1]
      b2 = b2[,,1]
      b  = as.vector(b)
    }
  }
  if (is.null(m$mixture_posterior_weights[[1]]))
    mixture_weights = as.numeric(NA)
  else {
    if (length(dim(m$mixture_posterior_weights[[1]])) < 2)
      mixture_weights = t(do.call(cbind, m$mixture_posterior_weights))
    else
      mixture_weights = aperm(abind(m$mixture_posterior_weights,along = 3),
                                    c(3,1,2)) # L x J x R
  }
  if (is.null(m$clfsr[[1]]))
    clfsr = as.numeric(NA)
  else {
    if (length(dim(m$clfsr[[1]])) < 2)
      clfsr = t(do.call(cbind, m$clfsr))
    else
      clfsr = aperm(abind(m$clfsr,along = 3),c(3,1,2)) # L x J x R
  }
  s = list(
        alpha  = t(m$pip),
        b1     = b1,
        b2     = b2,
        KL     = m$kl,
        lbf    = m$lbf,
        V      = m$prior_variance,
        sigma2 = d$residual_variance,
        elbo   = m$get_objective(dump = TRUE),
        niter  = m$niter,
        convergence        = m$convergence,
        coef               = d$rescale_coef(b),
        b1_rescaled        = rescale_single_effects(b1,d$rescale_coef),
        mixture_weights    = mixture_weights,
        conditional_lfsr   = clfsr,
        lfsr               = mvsusie_get_lfsr(clfsr, t(m$pip)),
        single_effect_lfsr = mvsusie_single_effect_lfsr(clfsr, t(m$pip)))
    if (!is.null(m$pip_history))
      s$alpha_history = m$pip_history
    if (!is.null(m$lbf_history))
      s$lbf_history = m$lbf_history
    if (!is.null(m$prior_history))
      s$prior_history = m$prior_history
  
    # FIXME: unit test for scaling issue for the fitted
    if (inherits(d, "RSSData"))
      s$fitted = d$XtX %*% b
    else
      s$fitted = d$compute_Xb(b)
    if (is.null(dim(s$coef)))
      s$intercept = s$coef[1]
    else
      s$intercept = s$coef[1,]
    if (estimate_prior_variance)
      s$V = m$prior_variance
    class(s) = "susie"
    return(s)
}

# This is used by report_susie_model to rescale the L x J x R matrix
# of single effect estimates.
rescale_single_effects <- function (b1, rescale_coef) {
  L <- dim(b1)[1]
  J <- dim(b1)[2]
  if (is.matrix(b1))
    R <- 1
  else
    R <- dim(b1)[3]
  out <- array(0,c(L,J + 1,R))
  for (l in 1:L)
    if (is.matrix(b1))
      out[l,,] = rescale_coef(b1[l,])
    else
      out[l,,] = rescale_coef(b1[l,,])
  return(drop(out))
}
