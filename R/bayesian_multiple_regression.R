#' @title Bayesian multiple regression object
#' @importFrom R6 R6Class
#' @keywords internal
BayesianMultipleRegression <- R6Class("BayesianMultipleRegression",
  public = list(
    initialize = function(J, residual_variance, prior_variance, estimate_prior_variance =FALSE) {
      private$J = J
      private$.prior_variance = prior_variance
      private$.residual_variance = residual_variance
      private$.posterior_b1 = matrix(0, J, 1)
      private$.posterior_b2 = matrix(0, J, 1)
      private$estimate_prior_variance = estimate_prior_variance

    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE) {
      # d: data object
      # use_residual: fit with residual instead of with Y, 
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d$XtR
      else XtY = d$XtY
      # OLS estimates
      betahat = 1/d$d * XtY
      shat2 = private$.residual_variance / d$d
      # deal with prior variance: can be "estimated" across effects
      if(private$estimate_prior_variance) {
        private$.prior_variance = est.prior.variance(betahat,shat2,prior_weights)
      }
      # posterior
      post_var = (1/private$.prior_variance + d$d/private$.residual_variance)^(-1) # posterior variance
      private$.posterior_b1 = (1/private$.residual_variance) * post_var * XtY
      private$.posterior_b2 = post_var + private$.posterior_b1^2 # second moment
      # Bayes factor
      private$.lbf = dnorm(betahat,0,sqrt(private$.prior_variance+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
      private$.lbf[shat2==Inf] == 0
    },
    compute_loglik_null = function(d) {
      if (inherits(d, "DenseData")) {
        private$.loglik_null = dnorm(d$Y,0,sqrt(private$.residual_variance),log=TRUE)
      } else {
        private$.loglik_null = NA
      }
    }
  ),
  private = list(
    estimate_prior_variance = FALSE,
    J = NULL,
    .prior_variance = NULL, # prior on effect size
    .residual_variance = NULL,
    .loglik_null = NULL,
    .lbf = NULL, # log Bayes factor
    .posterior_b1 = NULL, # posterior first moment
    .posterior_b2 = NULL, # posterior second moment
    denied = function(v) stop(paste0('$', v, ' is read-only'), call. = FALSE) 
  ),
  active = list(
    loglik_null = function(v) {
      if (missing(v)) private$.loglik_null
      else private$denied('loglik_null')
    },
    posterior_b1 = function(v) { 
      if (missing(v)) private$.posterior_b1
      else private$denied('posterior_b1')
    },
    posterior_b2 = function(v) { 
      if (missing(v)) private$.posterior_b2
      else private$denied('posterior_b2')
    },
    lbf = function(v) {
      if (missing(v)) private$.lbf
      else private$denied('lbf')
    }, 
    prior_variance = function(v) {
      if (missing(v)) {
        private$.prior_variance
      } else {
        if (length(v) != 2) stop('need argument of length 2')
        private$.prior_variance = v[1] 
        private$estimate_prior_variance = v[2]
      }
    },
    residual_variance = function(v) {
      if (missing(v)) private$.residual_variance
      else private$.residual_variance = v
    }
  )
)

# vector of gradients of logBF_j for each j, with respect to prior variance V
lbf.grad = function(V,shat2,T2){
  l = 0.5* (1/(V+shat2)) * ((shat2/(V+shat2))*T2-1)
  l[is.nan(l)] = 0
  return(l)
}

loglik.grad = function(V,betahat,shat2,prior_weights) {
  #log(bf) on each effect 
  lbf = dnorm(betahat,0,sqrt(V+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
  lbf[shat2==Inf] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)
  alpha = safe_comp_weight(lbf, prior_weights)$alpha
  sum(alpha*lbf.grad(V,shat2,betahat^2/shat2))
}

# define gradient as function of lV:=log(V)
# to improve numerical optimization
negloglik.grad.logscale = function(lV,betahat,shat2,prior_weights) {
  -exp(lV)*loglik.grad(exp(lV),betahat,shat2,prior_weights)
}

est.prior.variance = function(betahat,shat2,prior_weights) {
  if(loglik.grad(0,betahat,shat2,prior_weights)<0){
    return(0)
  } else {
    ##V.o = optim(par=log(V),fn=negloglik.logscale,gr = negloglik.grad.logscale,betahat=betahat,shat2=shat2,prior_weights=prior_weights,method="BFGS")
    ##if(V.o$convergence!=0){
    ##  warning("optimization over prior variance failed to converge")
    ##}
    V.u=uniroot(negloglik.grad.logscale,c(-10,10),extendInt = "upX",betahat=betahat,shat2=shat2,prior_weights=prior_weights)
    return(exp(V.u$root))
  }
}