#' @title Bayesian multiple regression object
#' @importFrom R6 R6Class
#' @keywords internal
BayesianSimpleRegression <- R6Class("BayesianSimpleRegression",
  public = list(
    initialize = function(J, residual_variance, prior_variance, estimate_prior_variance = FALSE) {
      private$J = J
      private$.prior_variance = prior_variance
      private$.residual_variance = residual_variance
      private$.posterior_b1 = matrix(0, J, 1)
      private$estimate_prior_variance = estimate_prior_variance
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      bhat = 1/d$X2_sum * XtY
      sbhat2 = private$.residual_variance / d$X2_sum
      if (save_summary_stats) {
        private$.bhat = bhat
        private$.sbhat = sqrt(sbhat2)
      }
      # deal with prior variance: can be "estimated" across effects
      if(private$estimate_prior_variance) {
          if (is.null(prior_weights)) prior_weights = rep(1/private$J, private$J)
        private$.prior_variance = est.prior.variance(bhat,sbhat2,prior_weights,method='optim')
      }
      # posterior
      post_var = (1/private$.prior_variance + d$X2_sum/private$.residual_variance)^(-1) # posterior variance
      private$.posterior_b1 = (1/private$.residual_variance) * post_var * XtY
      private$.posterior_b2 = post_var + private$.posterior_b1^2 # second moment
      # Bayes factor
      private$.lbf = dnorm(bhat,0,sqrt(private$.prior_variance+sbhat2),log=TRUE) - dnorm(bhat,0,sqrt(sbhat2),log=TRUE)
      private$.lbf[sbhat2==Inf] == 0
    },
    compute_loglik_null = function(d) {
      if (inherits(d, "DenseData")) {
        private$.loglik_null = dnorm(d$Y,0,sqrt(private$.residual_variance),log=TRUE)
      } else {
        private$.loglik_null = NA
      }
    }
  ),
  active = list(
    loglik_null = function() private$.loglik_null,
    posterior_b1 = function() private$.posterior_b1,
    posterior_b2 = function() private$.posterior_b2,
    lbf = function() private$.lbf,
    bhat = function() private$.bhat,
    sbhat = function() private$.sbhat,
    prior_variance = function() private$.prior_variance,
    residual_variance = function(v) {
      if (missing(v)) private$.residual_variance
      else private$.residual_variance = v
    }
  ),
  private = list(
    estimate_prior_variance = NULL,
    J = NULL,
    .bhat = NULL,
    .sbhat = NULL,
    .prior_variance = NULL, # prior on effect size
    .residual_variance = NULL,
    .loglik_null = NULL,
    .lbf = NULL, # log Bayes factor
    .posterior_b1 = NULL, # posterior first moment
    .posterior_b2 = NULL # posterior second moment
  )
)

loglik = function(V,betahat,shat2,prior_weights) {
  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  #log(bf) on each SNP
  lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)
  maxlbf = max(lbf)
  w = exp(lbf-maxlbf) # w =BF/BFmax
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  return(log(weighted_sum_w)+ maxlbf)
}

neg.loglik.logscale = function(lV,betahat,shat2,prior_weights){
  return(-loglik(exp(lV),betahat,shat2,prior_weights))
}

# vector of gradients of logBF_j for each j, with respect to prior variance V
lbf.grad = function(V,sbhat2,T2){
  l = 0.5* (1/(V+sbhat2)) * ((sbhat2/(V+sbhat2))*T2-1)
  l[is.nan(l)] = 0
  return(l)
}

loglik.grad = function(V,bhat,sbhat2,prior_weights) {
  #log(bf) on each effect
  lbf = dnorm(bhat,0,sqrt(V+sbhat2),log=TRUE) - dnorm(bhat,0,sqrt(sbhat2),log=TRUE)
  lbf[sbhat2==Inf] = 0 # deal with special case of infinite sbhat2 (eg happens if X does not vary)
  alpha = safe_compute_weight(lbf, prior_weights)$alpha
  sum(alpha*lbf.grad(V,sbhat2,bhat^2/sbhat2))
}
# define gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale = function(lV,betahat,shat2,prior_weights) {
  -exp(lV)*loglik.grad(exp(lV),betahat,shat2,prior_weights)
}

est.prior.variance = function(betahat,shat2,prior_weights, method=c('optim', 'uniroot')) {
  if(method=="optim"){
    lV = optim(par=log(max(c(betahat^2-shat2, 1), na.rm = TRUE)), fn=neg.loglik.logscale, betahat=betahat, shat2=shat2, prior_weights = prior_weights, method='Brent', lower = -30, upper = 15)$par
    V = exp(lV)
  } else {
    V.u = uniroot(negloglik.grad.logscale,c(-10,10),extendInt = "upX",betahat=betahat,shat2=shat2,prior_weights=prior_weights)
    V = exp(V.u$root)
  }
  if(loglik(0,betahat,shat2,prior_weights) >= loglik(V,betahat,shat2,prior_weights)) V=0 # set V exactly 0 if that beats the numerical value
  return(V)
}