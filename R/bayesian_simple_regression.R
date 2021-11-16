# Bayesian multiple regression object.
#
#' @importFrom R6 R6Class
#' @importFrom stats dnorm
#' @importFrom stats uniroot
#' @importFrom stats optim
BayesianSimpleRegression = R6Class("BayesianSimpleRegression",
  public = list(
    initialize = function (J, prior_variance) {
      private$J                     = J
      private$prior_variance_scalar = prior_variance
      private$.posterior_b1         = matrix(0,J,1)
    },
      
    fit = function (d, prior_weights = NULL, use_residual = FALSE,
                    save_summary_stats = FALSE, save_var = FALSE,
                    estimate_prior_variance_method = NULL,
                    check_null_threshold = 0) {
      # d: data object
      # use_residual: fit with residual instead of with Y,
      # a special feature for when used with SuSiE algorithm
        
      # bhat is J by R
      bhat = d$get_coef(use_residual)
      if(is.numeric(d$svs))
        # X2_sum is a length-J vector.
        sbhat2 = d$sbhat^2
      else
        sbhat2 = matrix(unlist(d$svs),ncol = 1)
      sbhat2[which(is.nan(sbhat2) | is.infinite(sbhat2))] = 1e6
      if (save_summary_stats) {
        private$.bhat  = bhat
        private$.sbhat = sqrt(sbhat2)
      }
      
      # Deal with prior variance: can be "estimated" across effects.
      if (!is.null(estimate_prior_variance_method)) {
        if (estimate_prior_variance_method == "EM")
          private$cache = list(b=bhat, s=sbhat2)
        else
          private$prior_variance_scalar =
            private$estimate_prior_variance(bhat,sbhat2,prior_weights,
              method = estimate_prior_variance_method,
              check_null_threshold = check_null_threshold)
       }
      
      # Posterior calculations.
      post_var = 1/(1/private$prior_variance_scalar + 1/sbhat2)
      if (save_var)
        private$.posterior_variance = post_var
      private$.posterior_b1 = post_var * bhat/sbhat2
      private$.posterior_b2 = post_var + private$.posterior_b1^2
      
      # Bayes factor.
      private$.lbf = dnorm(bhat,0,sqrt(private$prior_variance_scalar+sbhat2),
                       log = TRUE) - dnorm(bhat,0,sqrt(sbhat2),log = TRUE)
      if (!is.null(ncol(private$.lbf)) && ncol(private$.lbf) == 1)
        private$.lbf = as.vector(private$.lbf)
      private$.lbf[is.infinite(sbhat2)] = 0
    },
    set_thread = function  (value) private$n_thread = value
  ),
    
  active = list(
    loglik_null  = function() private$.loglik_null,
    posterior_b1 = function() private$.posterior_b1,
    posterior_b2 = function() private$.posterior_b2,
    lbf          = function() private$.lbf,
    bhat         = function() private$.bhat,
    sbhat        = function() private$.sbhat,
    prior_variance = function (v) {
      if (missing(v))
        private$prior_variance_scalar
      else
        private$prior_variance_scalar = v
    },
    posterior_variance = function() private$.posterior_variance
  ),
    
  private = list(
    J                     = NULL,
    .bhat                 = NULL,
    .sbhat                = NULL,
    prior_variance_scalar = NULL, # prior on effect size
    .loglik_null          = NULL,
    .lbf                  = NULL, # log Bayes factor
    .posterior_b1         = NULL, # posterior first moment
    .posterior_b2         = NULL, # posterior second moment
    .posterior_variance   = NULL, # posterior second moment
    cache                 = NULL, # some cached data
    n_thread = 4,
      
    loglik = function (V, betahat, shat2, prior_weights) {
          
      # log(bf) of each SNP.
      lbf = dnorm(betahat,0,sqrt(V + shat2),log = TRUE) -
          dnorm(betahat,0,sqrt(shat2),log = TRUE)

      # Deal with special case of infinite shat2 (e.g., happens if X
      # does not vary).
      lbf[is.infinite(shat2)] = 0
      return(compute_softmax(lbf,prior_weights)$log_sum)
    },
      
    neg_loglik_logscale = function (lV, betahat, shat2, prior_weights)
      -private$loglik(exp(lV),betahat,shat2,prior_weights),
      
    # Vector of gradients of logBF_j for each j, with respect to prior
    # variance V.
    lbf_grad = function (V, sbhat2, T2) {
      l = (1/(V + sbhat2)) * ((sbhat2/(V + sbhat2))*T2 - 1)/2
      l[is.nan(l)] = 0
      return(l)
    },
      
    loglik_grad = function (V, bhat, sbhat2, prior_weights) {
          
      # Log(bf) on each effect.
      lbf = dnorm(bhat,0,sqrt(V + sbhat2),log = TRUE) -
            dnorm(bhat,0,sqrt(sbhat2),log = TRUE)

      # deal with special case of infinite sbhat2 (eg happens if X
      # does not vary)
      lbf[is.infinite(sbhat2)] = 0 
      alpha = compute_softmax(lbf, prior_weights)$weights
      return(sum(alpha * private$lbf_grad(V,sbhat2,bhat^2/sbhat2)))
    },
      
    # Define gradient as function of lV := log(V) to improve numerical
    # stability.
    negloglik_grad_logscale = function (lV, betahat, shat2, prior_weights)
      -exp(lV) * private$loglik_grad(exp(lV),betahat,shat2,prior_weights),

    estimate_prior_variance = function (betahat, shat2, prior_weights,
                                        method = c("optim","uniroot","simple"),
                                        check_null_threshold = 0) {
      if (is.null(prior_weights))
        prior_weights = rep(1/private$J,private$J)
      method = match.arg(method)
      if (method == "optim")
          
        # Method BFGS is 1.5x slower than Brent with upper 15 lower -15
        # although it does not require specifying upper/lower.
        V = private$estimate_prior_variance_optim(betahat,shat2,prior_weights,
                                                  method = "Brent",lower = -30,
                                                  upper = 15)
      else if (method == "uniroot") {
        V.u = uniroot(private$negloglik_grad_logscale,c(-10,10),
                      extendInt = "upX",betahat = betahat,shat2 = shat2,
                      prior_weights = prior_weights)
        V = exp(V.u$root)
      } else if (method == "simple")
        V = private$estimate_prior_variance_simple()
      else
        stop("Optimization method not supported")
      
      # Set V exactly to zero if that beats the numerical value by a
      # loglik factor of 1 + check_null_threshold.
      if (private$loglik(0,betahat,shat2,prior_weights) +
          check_null_threshold >= private$loglik(V,betahat,shat2,
                                                 prior_weights))
        V = 0
      return(V)
    },
      
    estimate_prior_variance_optim = function (betahat, shat2, prior_weights,
                                              ...) {
      lV = optim(par = log(max(c(betahat^2 - shat2,1),na.rm = TRUE)),
                 fn = private$neg_loglik_logscale,betahat = betahat,
                 shat2 = shat2,prior_weights = prior_weights, ...)$par
      if (private$neg_loglik_logscale(log(private$prior_variance_scalar),
                                      betahat,shat2,prior_weights) 
          < private$neg_loglik_logscale(lV,betahat,shat2,prior_weights))
        lV = log(private$prior_variance_scalar)
      return(exp(lV))
    },
      
    estimate_prior_variance_em = function (pip)
      sum(pip*private$.posterior_b2),
      
    estimate_prior_variance_simple = function() private$prior_variance_scalar
  )
)
