#' @title SUm of Single Effect (SuSiE) Regression of multivariate Y on X
#' @details Performs Bayesian multiple linear regression of Y on X.
#' That is, this function
#' fits the regression model Y= sum_l Xb_l + e, where elements of e are iid N(0,residual_variance) and the
#' sum_l b_l is a p vector of effects to be estimated.
#' The assumption is that each b_l has exactly one non-zero element, with all elements
#' equally likely to be non-zero.
#' @param X an n by p matrix of covariates
#' @param Y an n vector, or n by r matrix of response variables
#' @param L maximum number of non-zero effects
#' @param prior_variance Can be 1) a vector of length L, or a scalar, for scaled prior variance when Y is univariate (equivalent to `susieR::susie`);
#' 2) a matrix for simple Multivariate regression or 3) a MASH fit that contains an array of prior covariance matrices and their weights
#' @param residual_variance the residual variance (defaults to sample variance of Y)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting.
#' If you do not standardize you may need
#' to think more carefully about specifying the scale of prior variance.
#' Whatever the value of standardize, the coefficients (returned from `coef`) are for X on the original input scale.
#' Any column of X that has zero variance is not standardized, but left as is.
#' @param intercept Should intercept be fitted (default=TRUE) or set to zero (FALSE). The latter is generally not recommended.
#' @param estimate_residual_variance indicates whether to estimate residual variance (currently only works for univariate Y input)
#' @param estimate_prior_variance indicates whether to estimate prior (currently only works for univariate Y and for multivariate Y when prior is a single matrix)
#' @param estimate_prior_method the method used for estimating prior variance: "optim", "uniroot" and "em" for univariate Y, "optim" and "simple" for multivariate Y.
#' @param check_null_threshold when prior variance is estimated, compare the estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example, you can set it to 0.1 to nudge the estimate towards zero. When used with "EM" method
#' setting \code{check_null_threshold=NA} will skip the check and instead relying solely on EM to update this parameter.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it.
#' @param precompute_covariances if TRUE, precomputes various covariance quantities to speed up computations at the cost of increased memory usage
#' @param s_init a previous model fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param compute_univariate_zscore if TRUE outputs z-score from per variable univariate regression
#' @param n_thread maximum number of threads to use for parallel computation (only applicable to mixture prior)
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param verbose if TRUE outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves some current quantities of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{b1}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{b2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{KL}{an L vector of KL divergence}
#' \item{lbf}{an L vector of logBF}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{niter}{number of iterations took for convergence}
#' \item{convergence}{convergence status}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{walltime}{records runtime of the fitting algorithm}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' res = msusie(X,y,L=10)
#'
#' @importFrom stats var
#' @importFrom susieR susie_get_cs 
#' @export
msusie = function(X,Y,L=10,
                 prior_variance=0.2,
                 residual_variance=NULL,
                 prior_weights=NULL,
                 standardize=TRUE,intercept=TRUE,
                 estimate_residual_variance=FALSE,
                 estimate_prior_variance=TRUE,
                 estimate_prior_method='EM',
                 check_null_threshold=0, prior_tol=1E-9,
                 compute_objective=FALSE,
                 s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                 compute_univariate_zscore = FALSE,
                 precompute_covariances = FALSE,
                 n_thread=1,max_iter=100,tol=1e-3,
                 verbose=TRUE,track_fit=FALSE,approximate=FALSE) {
  # adjust prior weights
  if (is.null(prior_weights)) prior_weights = c(rep(1/ncol(X), ncol(X)))
  else prior_weights = prior_weights / sum(prior_weights)
  if(inherits(prior_variance, 'MashInitializer')){
    prior_variance = prior_variance$clone(deep=T)
  }
  # adjust prior effects
  is_numeric_prior = !(is.matrix(prior_variance) || inherits(prior_variance, 'MashInitializer'))
  if (standardize && !is_numeric_prior) {
    # Scale prior variance
    # https://github.com/stephenslab/mmbr/blob/master/inst/prototypes/prior_matrices_scale.ipynb
    sigma = sapply(1:ncol(Y), function(i) sd(Y[,i], na.rm=T))
    n = sapply(1:ncol(Y), function(i) length(which(!is.na(Y[,1]))))
    sigma = sigma / sqrt(n)
    if (is.matrix(prior_variance)) prior_variance = scale_covariance(prior_variance, sigma)
    else prior_variance$scale_prior_variance(sigma)
  }
  # set data object
  if (any(is.na(Y))) {
    # When the residual variance is a diagonal matrix,
    # the approximate version has the same result as the exact version,
    # and it is faster, so we aet approximate = T
    if(isDiagonal(residual_variance)){
      approximate = TRUE
    }
    data = DenseDataYMissing$new(X, Y, approximate)
    estimate_residual_variance = FALSE
  } else {
    data = DenseData$new(X, Y)
  }
  # include residual variance in data
  data$set_residual_variance(residual_variance, numeric = is_numeric_prior, quantities = 'residual_variance')
  data$standardize(intercept, standardize)
  data$set_residual_variance(quantities='effect_variance')
  #
  s = mmbr_core(data, s_init, L, prior_variance, prior_weights,
            estimate_residual_variance, estimate_prior_variance, estimate_prior_method, check_null_threshold,
            precompute_covariances, compute_objective, n_thread, max_iter, tol, prior_tol, track_fit, verbose)
  # CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage=coverage, X=X, min_abs_corr=min_abs_corr)
  }
  # report z-scores from univariate regression
  if (compute_univariate_zscore) {
    s$z = susieR:::calc_z(X,Y,center=intercept,scale=standardize)
  }
  return(s)
}

#' @title SUm of Single Effect (SuSiE) Regression using Sufficient Statistics XtX, XtY, YtY, N
#' @param XtX a J by J matrix
#' @param XtY a J by R matrix
#' @param YtY an R by R matrix
#' @param N sample size
#' @param L maximum number of non-zero effects
#' @param prior_variance Can be 1) a vector of length L, or a scalar, for scaled prior variance when Y is univariate (equivalent to `susieR::susie`); 2) a matrix for simple Multivariate regression or 3) a MASH fit that contains an array of prior covariance matrices and their weights
#' @param residual_variance the residual variance (defaults to 1)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param standardize logical flag (default=TRUE) for whether to standardize columns of X to unit variance prior to fitting
#' @param estimate_residual_variance indicates whether to estimate residual variance (currently only works for univariate Y input)
#' @param estimate_prior_variance indicates whether to estimate prior (currently only works for univariate Y and for multivariate Y when prior is a single matrix)
#' @param estimate_prior_method the method used for estimating prior variance: "optim", "uniroot" and "em" for univariate Y, "optim" and "simple" for multivariate Y.
#' @param check_null_threshold when prior variance is estimated, compare the estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example, you can set it to 0.1 to nudge the estimate towards zero. When used with "EM" method
#' setting \code{check_null_threshold=NA} will skip the check and instead relying solely on EM to update this parameter.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it.
#' @param s_init a previous model fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param n_thread maximum number of threads to use for parallel computation (only applicable to mixture prior)
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param verbose if TRUE outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves some current quantities of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{b1}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{b2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{KL}{an L vector of KL divergence}
#' \item{lbf}{an L vector of logBF}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{niter}{number of iterations took for convergence}
#' \item{convergence}{convergence status}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{walltime}{records runtime of the fitting algorithm}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' X = t(t(X) - colMeans(X))
#' y = t(t(y) - colMeans(y))
#' XtX = crossprod(X)
#' XtY = crossprod(X, y)
#' YtY = crossprod(y)
#' res = msusie_ss(XtX,XtY,YtY,n,L=10)
#'
#' @importFrom stats var
#' @importFrom susieR susie_get_cs 
#' @export
msusie_suff_stat = function(XtX, XtY, YtY, N, L=10,
                            prior_variance=50,
                            residual_variance=NULL,
                            prior_weights=NULL,
                            standardize=TRUE,
                            estimate_residual_variance=FALSE,
                            estimate_prior_variance=TRUE,
                            estimate_prior_method='EM',
                            check_null_threshold=0, prior_tol=1E-9,
                            compute_objective=FALSE,
                            precompute_covariances = FALSE,
                            s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                            n_thread=1, max_iter=100,tol=1e-3,
                            verbose=TRUE,track_fit=FALSE) {
  if (is.null(prior_weights)) prior_weights = c(rep(1/ncol(XtX), ncol(XtX)))
  else prior_weights = prior_weights / sum(prior_weights)
  
  if(inherits(prior_variance, 'MashInitializer')){
    prior_variance = prior_variance$clone(deep=T)
  }
  is_numeric_prior = !(is.matrix(prior_variance) || inherits(prior_variance, 'MashInitializer'))
  if (standardize && !is_numeric_prior) {
    # Scale prior variance
    # https://github.com/stephenslab/mmbr/blob/master/inst/prototypes/prior_matrices_scale.ipynb
    sigma = sqrt(diag(YtY)/(N-1))
    sigma = sigma / sqrt(N)
    if (is.matrix(prior_variance)) prior_variance = scale_covariance(prior_variance, sigma)
    else prior_variance$scale_prior_variance(sigma)
  }
  
  data = SSData$new(XtX, XtY, YtY, N)
  # include residual variance in data
  data$set_residual_variance(residual_variance, numeric = is_numeric_prior,
                             quantities = 'residual_variance')
  data$standardize(standardize)
  data$set_residual_variance(quantities='effect_variance')
  #
  s = mmbr_core(data, s_init, L, prior_variance, prior_weights,
                estimate_residual_variance, estimate_prior_variance, estimate_prior_method, check_null_threshold,
                precompute_covariances, compute_objective, n_thread, max_iter, tol, prior_tol, track_fit, verbose)
  # CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage=coverage, Xcorr=cov2cor(XtX), min_abs_corr=min_abs_corr)
  }
  return(s)
}


#' @title SUm of Single Effect (SuSiE) Regression using Summary Statistics Z and R
#' @param Z a J by R matrix of z scores
#' @param R a J by J LD matrix
#' @param eigenR a list containing eigenvalues (values) and eigenvectors (vectors) of R. One of R and eigenR must be specified.
#' @param L maximum number of non-zero effects
#' @param prior_variance Can be 1) a vector of length L, or a scalar, for scaled prior variance when Y is univariate (equivalent to `susieR::susie`); 2) a matrix for simple Multivariate regression or 3) a MASH fit that contains an array of prior covariance matrices and their weights
#' @param residual_variance the residual variance (defaults to 1)
#' @param prior_weights a p vector of prior probability that each element is non-zero
#' @param estimate_residual_variance indicates whether to estimate residual variance (currently only works for univariate Y input)
#' @param estimate_prior_variance indicates whether to estimate prior (currently only works for univariate Y and for multivariate Y when prior is a single matrix)
#' @param estimate_prior_method the method used for estimating prior variance: "optim", "uniroot" and "em" for univariate Y, "optim" and "simple" for multivariate Y.
#' @param check_null_threshold when prior variance is estimated, compare the estimate with the null and set prior variance to null (zero) unless the log-likelihood
#' using the estimate is larger than that of null by this threshold. For example, you can set it to 0.1 to nudge the estimate towards zero. When used with "EM" method
#' setting \code{check_null_threshold=NA} will skip the check and instead relying solely on EM to update this parameter.
#' @param prior_tol when prior variance is estimated, compare the estimated value to this tol at the end of
#' the analysis and exclude a single effect from PIP computation if the estimated prior variance is smaller than it.
#' @param s_init a previous model fit with which to initialize
#' @param coverage coverage of confident sets. Default to 0.95 for 95\% credible interval.
#' @param min_abs_corr minimum of absolute value of correlation allowed in a credible set.
#' Default set to 0.5 to correspond to squared correlation of 0.25,
#' a commonly used threshold for genotype data in genetics studies.
#' @param n_thread maximum number of threads to use for parallel computation (only applicable to mixture prior)
#' @param max_iter maximum number of iterations to perform
#' @param tol convergence tolerance
#' @param z_thresh the z score threshold below which to call an effect null
#' @param verbose if TRUE outputs some progress messages
#' @param track_fit add an attribute \code{trace} to output that saves some current quantities of all iterations
#' @return a susie fit, which is a list with some or all of the following elements\cr
#' \item{alpha}{an L by p matrix of posterior inclusion probabilites}
#' \item{b1}{an L by p matrix of posterior means (conditional on inclusion)}
#' \item{b2}{an L by p matrix of posterior second moments (conditional on inclusion)}
#' \item{KL}{an L vector of KL divergence}
#' \item{lbf}{an L vector of logBF}
#' \item{sigma2}{residual variance}
#' \item{V}{prior variance}
#' \item{elbo}{a vector of values of elbo achieved (objective function)}
#' \item{niter}{number of iterations took for convergence}
#' \item{convergence}{convergence status}
#' \item{sets}{a list of `cs`, `purity` and selected `cs_index`}
#' \item{pip}{a vector of posterior inclusion probability}
#' \item{walltime}{records runtime of the fitting algorithm}
#' \item{z}{a vector of univariate z-scores}
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' R = t(X) %*% X
#' z = susieR:::calc_z(X,y)
#' res = msusie_rss(z,R,L=10)
#'
#' @importFrom stats var
#' @importFrom susieR susie_get_cs 
#' @export
msusie_rss = function(Z,R=NULL,eigenR=NULL,L=10,r_tol = 1e-08,
                      prior_variance=50,
                      residual_variance=NULL,
                      prior_weights=NULL,
                      estimate_residual_variance=FALSE,
                      estimate_prior_variance=TRUE,
                      estimate_prior_method='EM',
                      check_null_threshold=0, prior_tol=1E-9,
                      compute_objective=FALSE,
                      precompute_covariances = FALSE,
                      s_init = NULL,coverage=0.95,min_abs_corr=0.5,
                      n_thread=1, max_iter=100,tol=1e-3,z_thresh = 2,
                      verbose=TRUE,track_fit=FALSE) {
  
  if (is.null(prior_weights)){
    if(is.null(dim(Z))){
      prior_weights = c(rep(1/length(Z), length(Z)))
    }else{
      prior_weights = c(rep(1/nrow(Z), nrow(Z)))
    }
  }
  else prior_weights = prior_weights / sum(prior_weights)

  if(inherits(prior_variance, 'MashInitializer')){
    prior_variance = prior_variance$clone(deep=T)
  }
  
  data = RSSData$new(Z, R, eigenR, r_tol)
  if (is.null(residual_variance)) {
    residual_variance = diag(data$n_condition)
  }
  #
  data$set_residual_variance(residual_variance, numeric = !(is.matrix(prior_variance) || inherits(prior_variance, 'MashInitializer')))
  s = mmbr_core(data, s_init, L, prior_variance, prior_weights,
                estimate_residual_variance, estimate_prior_variance, estimate_prior_method, check_null_threshold,
                precompute_covariances, compute_objective, n_thread, max_iter, tol, prior_tol, track_fit, verbose)
  # CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage=coverage, Xcorr=data$XtX, min_abs_corr=min_abs_corr)
  }
  return(s)
}

#' @title Core MMBR code
#' @importFrom susieR susie_get_pip 
#' @keywords internal
mmbr_core = function(data, s_init, L, prior_variance, prior_weights,
            estimate_residual_variance, estimate_prior_variance, estimate_prior_method, check_null_threshold,
            precompute_covariances, compute_objective, n_thread, max_iter, tol, prior_tol, track_fit, verbose) {
  start_time = proc.time()
  # for now the type of prior_variance controls the type of regression
  if (is.numeric(prior_variance)) {
    n_thread = NULL
    if (data$n_condition > 1 && !is.matrix(prior_variance))
      stop(paste("prior variance cannot be a number for multivariate analysis with", data$n_condition, "response variables."))
    if (is.matrix(prior_variance)) {
      base = BayesianMultivariateRegression
    } else {
      base = BayesianSimpleRegression
      prior_variance = prior_variance * data$residual_variance
    }
  } else {
    if (inherits(prior_variance, 'MashInitializer')) {
      if (prior_variance$n_condition != data$n_condition) stop("Dimension mismatch between input prior covariance and response variable data.")
      base = MashRegression
      if (!precompute_covariances)
        warning("precompute_covariances option is set to FALSE by default to save memory usage with MASH prior. The computation can be a lot slower as a result. It is recommended that you try setting it to TRUE, see if there is a memory usage issue and only switch back if it is a problem.")
      if (precompute_covariances)
        prior_variance$precompute_cov_matrices(data)
      if (estimate_prior_variance && !is.null(estimate_prior_method) && estimate_prior_method == 'EM')
        prior_variance$compute_prior_inv()
    } else if (is.null(prior_variance)) {
      # FIXME: a temporary interface for MS for methylation data
      base = NIGMGRegression
      prior_variance = data$n_condition 
    } else {
      stop("prior_variance should be a scalar for univariate response, or a matrix or MashInitializer object for multivariate response.")
    }
  }
  if (!estimate_prior_variance) estimate_prior_method = NULL
  # Below are the core computations
  SER_model = SingleEffectModel(base)$new(data$n_effect, prior_variance)
  if (!is.null(n_thread)) {
    if (n_thread<1) stop("number of threads cannot be smaller than 1")
    SER_model$set_thread(floor(n_thread))
  }
  SuSiE_model = SuSiE$new(SER_model, L, estimate_residual_variance, compute_objective, max_iter, tol, track_pip=track_fit, track_lbf=track_fit, track_prior=track_fit)
  if (!is.null(s_init)) SuSiE_model$init_from(s_init)
  SuSiE_model$fit(data, prior_weights, estimate_prior_method, check_null_threshold, verbose)
  s = report_susie_model(data, SuSiE_model, estimate_prior_variance)
  s$pip = susie_get_pip(s, prior_tol=prior_tol)
  ## clean up prior object
  if (inherits(prior_variance, "MashInitializer")) prior_variance$remove_precomputed()
  s$walltime = proc.time() - start_time
  return(s)
}