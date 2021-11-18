# @title Compute condition specific posterior inclusion probability.
# @description This is only relevant when canonical priors are used
# @param m M&M model
# @param prior_obj prior mixture object
# @return P by R matrix of PIP per condition
# @keywords internal
mvsusie_get_pip_per_condition = function (m, prior_obj) {
  condition_pip = mvsusie_get_alpha_per_condition(m,prior_obj)
  return(do.call(cbind,lapply(1:dim(condition_pip)[3],
           function(r) apply(condition_pip[,,r],2,function(x) 1-prod(1-x)))))
}

# Compute condition specific posterior inclusion probability per
# effect.
mvsusie_get_alpha_per_condition = function (m, prior_obj) {
  condition_indicator = do.call(rbind,
    lapply(1:length(prior_obj$prior_variance$xUlist),
      function(i) as.integer(diag(prior_obj$prior_variance$xUlist[[i]]) != 0)))
  condition_pip = array(0,dim = dim(m$b1))
  for (r in 1:dim(condition_pip)[3]) {
    for (p in 1:length(condition_indicator[,r]))
      condition_pip[,,r] =
        condition_pip[,,r] + m$mixture_weights[,,p] * condition_indicator[p,r]
    condition_pip[,,r] = condition_pip[,,r] * m$alpha
  }
  return(condition_pip)
}

# Multiviate regression object with NIG-MG prior
#
#' @importFrom R6 R6Class
#' @importFrom stats model.matrix
#' @importFrom stats lm
#' @importFrom stats rgamma
#' 
NIGMGRegression <- R6Class("NIGMGRegression",
  inherit = BayesianMultivariateRegression,
  public = list(
    # FIXME: have to make a prior variance object for this prior
    initialize = function(J, R) {
      private$J = J
      # FIXME: need to initialize with prior variance which will contain info on R
      private$.posterior_b1 = matrix(0, J, R)
    },
    fit = function(d, prior_weights = NULL, use_residual = FALSE, save_summary_stats = FALSE, save_var = FALSE, estimate_prior_variance_method = NULL, check_null_threshold=0) {
      # FIXME: need data object for DWT 
      if (use_residual) Y = list(C=rep(0, nrow(d$Y)), D=d$residual[,-1], J = 7, filter.number=10,family="DaubLeAsymm")
      else Y = list(C=rep(0, nrow(d$Y)), D=d$Y[,-1], J = 7, filter.number=10,family="DaubLeAsymm") 
      class(Y) = "DWT"
      ans0   <- NULL # grove::Denoise(Y) 
      m   <- apply(d$X, 2, function(X) fit_MS(X,Y)) 
      private$.posterior_b1 = do.call(rbind,lapply(1:private$J, function(j) m[[j]]$effect))
      post_cov = matlist2array(lapply(1:private$J, function(j) m[[j]]$cov))
      private$.posterior_b2 = post_cov + matlist2array(lapply(1:nrow(private$.posterior_b1), function(i) tcrossprod(private$.posterior_b1[i,])))
      private$.lbf = unlist(lapply(1:P, function(j) m[[j]]$marginal_likelihood)) - ans0$marginal_likelihood
      # FIXME: need to handle baseline
      private$.baseline = do.call(cbind,lapply(1:P, function(j) m[[j]]$baseline))
    }
  ),
  private = list(
    .baseline = NULL
  )
)

  ###########
  #Description
  #function to obtain the fitted effect from a Single effect Wavelet regression
  #Depend on more parameters than the one currently in place, these parameters stay fix during a SuSIe.
  
  #Input
  #grove.obj an output from FAnova (from grove package)
  #n_samples  the number of Monte carlo simuation considered
  #Output: list
  #effect: vector of length 2^J (first parameter= scaling parameter), the other=  wavelet coeff
  #baseline: vector of length 2^J (first parameter= scaling parameter), the other=  wavelet coeff
  #cov_effect= variance covariance matrix of the estimated effect, not sure to double chec
  
  get_post <- function(grove.obj, n_samples, m)
  {
    D <- grove.obj$samples$mean
    sub_D <- t(D[,,1]) #sub_D[,2]  as wavelet coeff effect, sub_D[,1] has the baseline
    cov_D <- D[2,,1]
    for ( i in 2:n_samples)#n_samples being the number of sample of the posterior
    {
      sub_D  <- sub_D + t(D[,,i])
      cov_D  <- cbind(cov_D ,D[2,,i])
    }
    y <- grove.obj$data$W$C
    X <- model.matrix(grove.obj$data$formula, grove.obj$data$X)
    s20 <- summary(lm(y ~ X))$sigma
    nu0 <- 10
    N <- 1
    tempC <- matrix(as.numeric(NA),ncol = 1,nrow = N) 
    for (i in 1:N) {
      #average C coefficient for reg 2
      x <- c(1, 0)
      g <- length(y)
      S <- 1
      n <- dim(X)[1]
      p <- dim(X)[2]
      Hg <- (g/(g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
      SSRg <- t(y) %*% (diag(1, nrow = n) - Hg) %*% y
      s2 <- 1/rgamma(S, (nu0 + n)/2, (nu0 * s20 + SSRg)/2)
      Vb <- g * solve(t(X) %*% X)/(g + 1)
      Eb <- Vb %*% t(X) %*% y
      E <- matrix(rnorm(S * p, 0, sqrt(s2)), S, p)
      C <- t(t(E %*% chol(Vb)) + c(Eb))
      tempC[i,1]<- sum(C * x)
    }
    varC  <- var(tempC)#pas sur de ça
    
    #Sub_D first colum is b0 for wavelet coeff, and second is the regression effect of the covariate
    sub_D <- sub_D/n_samples
    
    # SNP effect
    SNP_effect <- sub_D[,2]
    Baseline_effect <- sub_D[,1]
    
    
    effectC <-grove.obj$C_hat[2]#scaling effect
    baselineC <- grove.obj$C_hat[1]
    # SNP effect
    
    
    #Building a matrix corresponding for the first position to the scaling coefficient free from the other d coefficients
    temp         <- diag(1, m)
    temp [-1,-1] <- cov(t(cov_D  ))
    temp[1,1]    <- varC
    cov_effect   <- temp
    
    
    effect <-  c(effectC,SNP_effect) 
    baseline <-  c(baselineC,Baseline_effect)
    
    out <- list( effect     = effect, 
                 baseline   = baseline,
                 cov_effect = cov_effect)
    return(out)
    
    
  }
  
  
  ###########
  #get_baseline
  ###########
  #Description
  #function to obtain the fitted baseline from a Single effect Wavelet regression
  #Depend on more parameters than the one currently in place, these parameters stay fix during a SuSIe.
  
  #Input
  #grove.obj an output from FAnova (from grove package)
  #n_samples  the number of Monte carlo simuation considered
  #Output:
  # vector of length 2^J (first parameter= scaling parameter), the other=  wavelet coeff
  get_baseline <- function(grove.obj,  n_samples=n_samples,N)
  {
    D <- grove.obj$samples$mean
    sub_D <- t(D[,,1])
    for ( i in 2:n_samples)#n_samples being the number of sample of the posterior
    {
      sub_D  <- sub_D + t(D[,,i])
    }
    #Sub_D first colum is b0 for wavelet coeff, and second is the regression effect of the covariate
    sub_D <- sub_D/n_samples
    y <- grove.obj$data$W$C
    X <- model.matrix(grove.obj$data$formula, grove.obj$data$X)
    s20 <- summary(lm(y ~ X))$sigma
    nu0 <- 10
    baseline <- sub_D[,1]
    tempC <- rep(NA, N)
    for (i in 1:N) {
      #average C coefficient for reg 2
      x <- c(1, 0)
      g <- length(y)
      S <- 1
      n <- dim(X)[1]
      p <- dim(X)[2]
      Hg <- (g/(g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)
      SSRg <- t(y) %*% (diag(1, nrow = n) - Hg) %*% y
      s2 <- 1/rgamma(S, (nu0 + n)/2, (nu0 * s20 + SSRg)/2)
      Vb <- g * solve(t(X) %*% X)/(g + 1)
      Eb <- Vb %*% t(X) %*% y
      E <- matrix(rnorm(S * p, 0, sqrt(s2)), S, p)
      C <- t(t(E %*% chol(Vb)) + c(Eb))
      tempC[i]<- sum(C * x)
    }
    tempC <- mean(tempC)#/2
    # SNP effect
    
    out <- c(tempC,baseline)
    return(out)
    
  }


      fit_MS <- function(X, Y)
      {
        
        temp_X <- data.frame ( factorA  = as.factor(X))
        temp_res <- NULL #  FAnova(Y, #regression
                         #    temp_X,
                         #    ~ 1 + factorA,
                         #    transition.mode = "Markov",
                         #    is.kappa.fixed=FALSE
                         # )                     

        # FIXME: n_samples is a variable to change
        fitted_effect <- get_post(temp_res, n_samples = nrow(Y$D), m=ncol(Y$D)+1)
        out <- list(  marginal_likelihood = temp_res$marginal_likelihood , 
                      baseline            = fitted_effect$baseline,# retriving baseline,
                      effect              = fitted_effect$effect,# retriving estimated effect
                      cov                 = fitted_effect$cov_effect
        )
        
        
        return(out)
      }
