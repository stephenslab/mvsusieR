#' @rdname mvsusie
#' 
#' @title Multivariate SUm of Single Effect (SuSiE) Regression
#' 
#' @description Performs a Bayesian multiple linear regression of Y on X.
#'   That is, this function fits the regression model \deqn{Y = \sum_l X
#'   b_l + e}, where the elements of \eqn{e} are \emph{i.i.d.} normal
#'   with zero mean and variance \code{residual_variance}, and the sum
#'   \eqn{\sum_l b_l} is a vector of p effects to be estimated. The
#'   SuSiE assumption is that each \eqn{b_l} has exactly one non-zero
#'   element.
#' 
#' @param X n by p matrix of covariates.
#' 
#' @param Y Vector of length n, or n by r matrix of response
#'   variables.
#' 
#' @param L Maximum number of non-zero effects.
#' 
#' @param prior_variance Can be either (1) a vector of length L, or a
#'   scalar, for scaled prior variance when Y is univariate (which
#'   should then be equivalent to \code{\link[susieR]{susie}}); or (2) a
#'   matrix for a simple multivariate regression; or (3) a MASH fit that
#'   contains an array of prior covariance matrices and their weights.
#' 
#' @param residual_variance The residual variance (defaults to the
#'   sample variance of Y).
#' 
#' @param prior_weights A vector of length p giving the prior
#'   probability that each element is non-zero. Note that the prior
#'   weights need to be non-negative but do not need to sum to 1; they
#'   will automatically be normalized to sum to 1 so that they represent
#'   probabilities. The default setting is that the prior weights are
#'   the same for all variables.
#' 
#' @param standardize Logical flag specifying whether to standardize
#'   columns of X to unit variance prior to fitting. If you do not
#'   standardize you may need to think more carefully about specifying
#'   the scale of the prior variance. Whatever the value of
#'   standardize, the coefficients (returned by \code{\link{coef}}) are
#'   for X on the original input scale. Note that any column of X with
#'   zero variance is not standardized, but left as is.
#' 
#' @param intercept Should intercept be fitted or set to zero. Setting
#'   \code{intercept = FALSE} is generally not recommended.
#' 
#' @param approximate Specifies whether to use approximate computation
#'   for the intercept when there are missing values in Y. The
#'   approximation saves some computational effort. Note that when the
#'   residual_variance is a diagonal matrix, running \code{mvsusie} with
#'   \code{approximate = TRUE} will give same result as
#'   \code{approximate = FALSE}, but with less running time. This
#'   setting is only relevant when there are missing values in Y and
#'   \code{intercept} = TRUE.
#' 
#' @param estimate_residual_variance When
#'   \code{estimate_residual_variance = TRUE}, the residual variance is
#'   estimated; otherwise it is fixed. Currently
#'   \code{estimate_residual_variance = TRUE} only works for univariate Y.
#' 
#' @param estimate_prior_variance When \code{estimate_prior_variance =
#'   TRUE}, the prior variance is estimated; otherwise it is
#'   fixed. Currently \code{estimate_prior_variance = TRUE} only works
#'   for univariate Y, or for multivariate Y when the prior variance is
#'   a matrix).
#' 
#' @param estimate_prior_method The method used for estimating the
#'   prior variance; valid choices are \code{"optim"}, \code{"uniroot"}
#'   or \code{"em"} for univariate Y; and \code{"optim"},
#'   \code{"simple"} for multivariate Y.
#' 
#' @param check_null_threshold When the prior variance is estimated,
#'   the estimate is compared against the null, and the prior variance
#'   is set to zero unless the log-likelihood using the estimate is
#'   larger than that of null by this threshold. For example, setting
#'   \code{check_null_threshold = 0.1} will \dQuote{nudge} the estimate
#'   towards zero. When used with \code{estimate_prior_method = "EM"},
#'   setting \code{check_null_threshold = NA} will skip this check.
#' 
#' @param prior_tol When the prior variance is estimated, compare the
#'   estimated value to this value at the end of the analysis and
#'   exclude a single effect from PIP computation if the estimated prior
#'   variance is smaller than it.
#'
#' @param compute_objective Add description of "compute_objective"
#'   input argument here.
#' 
#' @param precompute_covariances If \code{precompute_covariances =
#'   TRUE}, precomputes various covariance quantities to speed up
#'   computations at the cost of increased memory usage..
#' 
#' @param s_init A previous model fit with which to initialize.
#' 
#' @param coverage Coverage of confident sets.
#' 
#' @param min_abs_corr Minimum of absolute value of correlation
#'   allowed in a credible set. The setting \code{min_abs_corr = 0.5}
#'   corresponds to squared correlation of 0.25, which is a commonly
#'   used threshold for genotype data in genetics studies.
#' 
#' @param compute_univariate_zscore When
#'   \code{compute_univariate_zscore = TRUE}, the z-scores from the
#'   per-variable univariate regressions are outputted.
#' 
#' @param n_thread Maximum number of threads to use for parallel
#'   computation (only applicable when a mixture prior is used).
#' 
#' @param max_iter Maximum number of iterations to perform.
#' 
#' @param tol The model fitting will terminate when the increase in
#'   ELBOs between two successive iterations is less than \code{tol}.
#' 
#' @param verbosity Set \code{verbosity = 0} for no messages;
#'   \code{verbosity = 1} for a progress bar; and \code{verbosity = 2}
#'   for more detailed information about the algorithm's progress at the
#'   end of each iteration.
#' 
#' @param track_fit Add attribute \code{trace} to the return value
#'   which records the algorithm's progress at each iteration.
#' 
#' @return A multivariate susie fit, which is a list with some or all
#' of the following elements:
#' 
#' \item{alpha}{L by p matrix of posterior inclusion probabilites.}
#' 
#' \item{b1}{L by p matrix of posterior means (conditional on inclusion).}
#' 
#' \item{b2}{L by p matrix of posterior second moments (conditional on
#'   inclusion).}
#' 
#' \item{KL}{Vector of single-ffect KL divergences}
#' 
#' \item{lbf}{Vector of single-effect log-Bayes factors.}
#' 
#' \item{sigma2}{Residual variance.}
#' 
#' \item{V}{Prior variance.}
#' 
#' \item{elbo}{Vector storing the the evidence lower bound, or
#'   \dQuote{ELBO}, achieved at each iteration of the model fitting
#'   algorithm, which attempts to maximize the ELBO.}
#' 
#' \item{niter}{Number of iterations performed.}
#' 
#' \item{convergence}{Convergence status.}
#' 
#' \item{sets}{Estimated credible sets.}
#' 
#' \item{pip}{Vector of posterior inclusion probabilities.}
#' 
#' \item{walltime}{Records runtime of the model fitting algorithm.}
#' 
#' \item{z}{Vector of univariate z-scores.}
#' 
#' @examples
#' # Example with one response.
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow = n,ncol = p)
#' Y = X %*% beta + rnorm(n)
#' fit = mvsusie(X,Y,L = 10)
#'
#' # Example with three responses.
#' n = 500
#' p = 1000
#' true_eff = 2
#' X = matrix(sample(c(0,1,2),size = n*p,replace = TRUE),nrow = n,ncol = p)
#' beta1 = rep(0,p)
#' beta2 = rep(0,p)
#' beta3 = rep(0,p)
#' beta1[1:true_eff] = runif(true_eff)
#' beta2[1:true_eff] = runif(true_eff)
#' beta3[1:true_eff] = runif(true_eff)
#' y1 = X %*% beta1 + rnorm(n)
#' y2 = X %*% beta2 + rnorm(n)
#' y3 = X %*% beta3 + rnorm(n)
#' Y = cbind(y1,y2,y3)
#' prior = create_mash_prior(max_mixture_len = -1,
#'                           sample_data = list(X = X,Y = Y,
#'                                              residual_variance = cov(Y)))
#' fit = mvsusie(X,Y,prior_variance = prior)
#'
#' @importFrom Matrix isDiagonal
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom susieR susie_get_cs
#' 
#' @export
#' 
mvsusie = function (X, Y, L = 10, prior_variance = 0.2,
                    residual_variance = NULL, prior_weights = NULL,
                    standardize = TRUE,intercept = TRUE, approximate = FALSE,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method = "EM",
                    check_null_threshold = 0, prior_tol = 1e-9,
                    compute_objective = TRUE, s_init = NULL,
                    coverage = 0.95, min_abs_corr = 0.5,
                    compute_univariate_zscore = FALSE,
                    precompute_covariances = FALSE, n_thread = 1,
                    max_iter = 100, tol = 1e-3, verbosity = 2,
                    track_fit = FALSE) {
    
  # Adjust prior weights.
  if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X),ncol(X))
  else
    prior_weights = prior_weights/sum(prior_weights)
  
  # Check and process prior variance.
  if (inherits(prior_variance,"MashInitializer"))
    prior_variance = prior_variance$clone(deep = TRUE)
  is_numeric_prior = !(is.matrix(prior_variance) ||
                       inherits(prior_variance,"MashInitializer"))
  if (!is.null(dim(Y)) && ncol(Y) > 1 && is_numeric_prior)
    stop("Please specify prior variance for the multivariate response Y")
  if (standardize && !is_numeric_prior) {
      
    # Scale prior variance; see
    # https://github.com/stephenslab/mvsusieR/blob/master/
    # inst/prototypes/prior_matrices_scale.ipynb
    sigma = sapply(1:ncol(Y),function(i) sd(Y[,i],na.rm = TRUE))
    n     = sapply(1:ncol(Y),function(i) length(which(!is.na(Y[,1]))))
    sigma = sigma/sqrt(n)
    
    # Make sigma numerically more robust against extreme values.
    if (estimate_prior_variance)
      sigma = sigma/max(sigma)
    if (is.matrix(prior_variance))
      prior_variance = scale_covariance(prior_variance,sigma)
    else
      prior_variance$scale_prior_variance(sigma)
  }
  if (verbosity > 1) {
    message("Initializing data object...")
    message(paste("Dimension of X matrix:",nrow(X),ncol(X)))
    message(paste("Dimension of Y matrix:",nrow(Y),ncol(Y)))
  }
  
  # Set data object.
  if (any(is.na(Y))) {
      
    # When the residual variance is a diagonal matrix, the approximate
    # version has the same result as the exact version, and it is
    # faster, so we set approximate = TRUE.
    if (isDiagonal(residual_variance))
      approximate = TRUE
    data = DenseDataYMissing$new(X,Y,approximate)
    estimate_residual_variance = FALSE
  } else
    data = DenseData$new(X,Y)

  # Include residual variance in the data object.
  data$set_residual_variance(residual_variance,numeric = is_numeric_prior,
                             quantities = "residual_variance")
  data$standardize(intercept,standardize)
  data$set_residual_variance(quantities = "effect_variance")

  # Fit the susie model.
  s = mvsusie_core(data,s_init,L,prior_variance,prior_weights,
                   estimate_residual_variance,estimate_prior_variance,
                   estimate_prior_method,check_null_threshold,
                   precompute_covariances,compute_objective,n_thread,
                   max_iter,tol,prior_tol,track_fit,verbosity)
  
  # Compute CSs and PIPs.
  if (!is.null(coverage) && !is.null(min_abs_corr))
    s$sets = susie_get_cs(s,coverage = coverage,X = X,
                          min_abs_corr = min_abs_corr)
  
  # Report z-scores from univariate regression.
  if (compute_univariate_zscore)
    s$z = susieR:::calc_z(X,Y,center = intercept,scale = standardize)
  
  # FIXME: need a better approach.
  s$variable_names = colnames(X)
  s$condition_names = colnames(Y)
  return(s)
}
