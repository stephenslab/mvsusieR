# Base class for regression model
# Compare each vector of Y and X matrix at a time
BaseBayesianRegression <- R6Class("BaseBayesianRegression",
  public = list(
    lbf = NULL, # log Bayes factor
    posterior_b1 = NULL, # posterior first moment
    posterior_b2 = NULL, # posterior second moment
    initialize = function(prior=NULL) {
      self$set_prior(prior)
    },
    set_prior = function(prior=NULL) {
        # set prior
        if (is.null(prior)) {
            private$prior = -9
            private$estimate_prior = TRUE
        } else {
            private$prior = prior
            private$estimate_prior = FALSE
        }
    },
    set_residual_variance = function(r) private$residual_variance = r,
    get_residual_variance = function() private$residual_variance,
    get_prior = function() private$prior,
    fit = function(d, prior_weights = NULL, use_residual = FALSE) {
      # d: data object
      # use_residual: fit with residual instead of with Y, 
      # a special feature for when used with SuSiE algorithm
      if (use_residual) XtY = d.get_XtR()
      else XtY = d.get_XtY()
      # OLS estimates
      betahat = 1/d$d * XtY
      shat2 = private$residual_variance / d$d
      # deal with prior: can be "estimated" across effects
      if(private$estimate_prior) {
        private$prior = estimate.prior(betahat,shat2,prior_weights)
      }
      # posterior
      post_var = (1/private$prior + d$d/private$residual_variance)^(-1) # posterior variance
      self$posterior_b1 = (1/private$residual_variance) * post_var * XtY
      self$posterior_b2 = post_var + self$posterior_b1^2 # second moment
      # Bayes factor
      self$lbf = dnorm(betahat,0,sqrt(prior+shat2),log=TRUE) - dnorm(betahat,0,sqrt(shat2),log=TRUE)
      self$lbf[shat2==Inf] == 0
    }
  ),
  private = list(
    prior = NULL, # prior on effect size
    residual_variance = NULL,
    estimate_prior = FALSE,
    exit = function() {
      warning("Not implemented.")
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

estimate.prior = function(betahat,shat2,prior_weights) {
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