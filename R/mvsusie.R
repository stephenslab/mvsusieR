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
#'   probabilities. TO DO: Give default setting for
#'   \code{prior_weights}.
#' 
#' @param standardize Logical flag specifying whether to standardize
#'   columns of X to unit variance prior to fitting. If you do not
#'   standardize you may need to think more carefully about specifying
#'   the scale of the prior variance.  Whatever the value of
#'   standardize, the coefficients (returned by \code{\link{coef}}) are
#'   for X on the original input scale. Note that any column of X with
#'   zero variance is not standardized, but left as is.
#' 
#' @param intercept Should intercept be fitted or set to zero. Setting
#'   \code{intercept = FALSE} is generally not recommended.
#' 
#' @param approximate Specifies whether to use approximate computation
#'   for the intercept when there are missing values in Y.  The
#'   approximation saves some computational effort. Note that when the
#'   residual_variance is a diagonal matrix, running \code{mvsusie} with
#'   \code{approximate = TRUE} will give same result as
#'   \code{approximate = FALSE}, but with less running time. This
#'   setting is only relevant when there are missing values in Y and
#'   \code{intercept} = TRUE.
#' 
#' @param estimate_residual_variance indicates whether to estimate
#' residual variance (currently only works for univariate Y input)
#' 
#' @param estimate_prior_variance indicates whether to estimate prior
#' (currently only works for univariate Y and for multivariate Y when
#' prior is a single matrix)
#' 
#' @param estimate_prior_method the method used for estimating prior
#' variance: "optim", "uniroot" and "em" for univariate Y, "optim" and
#' "simple" for multivariate Y.
#' 
#' @param check_null_threshold when prior variance is estimated,
#' compare the estimate with the null and set prior variance to null
#' (zero) unless the log-likelihood using the estimate is larger than
#' that of null by this threshold. For example, you can set it to 0.1
#' to nudge the estimate towards zero. When used with "EM" method
#' setting \code{check_null_threshold=NA} will skip the check and
#' instead relying solely on EM to update this parameter.
#' 
#' @param prior_tol when prior variance is estimated, compare the
#' estimated value to this tol at the end of the analysis and exclude
#' a single effect from PIP computation if the estimated prior
#' variance is smaller than it.
#' 
#' @param precompute_covariances if TRUE, precomputes various
#' covariance quantities to speed up computations at the cost of
#' increased memory usage
#' 
#' @param s_init a previous model fit with which to initialize
#' 
#' @param coverage coverage of confident sets. Default to 0.95 for
#' 95\% credible interval.
#' 
#' @param min_abs_corr minimum of absolute value of correlation
#' allowed in a credible set.  Default set to 0.5 to correspond to
#' squared correlation of 0.25, a commonly used threshold for genotype
#' data in genetics studies.
#' 
#' @param compute_univariate_zscore if TRUE outputs z-score from per
#' variable univariate regression
#' 
#' @param n_thread maximum number of threads to use for parallel
#' computation (only applicable to mixture prior)
#' 
#' @param max_iter Maximum number of iterations to perform.
#' 
#' @param tol Convergence tolerance.
#' 
#' @param verbosity Set \code{verbosity = 0} for no messages;
#' \code{verbosity = 1} for a progress bar; and \code{verbosity = 2}
#' for more detailed information about the algorithm's progress at the
#' end of each iteration.
#' 
#' @param track_fit add an attribute \code{trace} to output that saves some current quantities of all iterations
#' 
#' @return a susie fit, which is a list with some or all of the
#' following elements
#' 
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
#' 
#' @examples
#' set.seed(1)
#' n = 1000
#' p = 1000
#' beta = rep(0,p)
#' beta[1:4] = 1
#' X = matrix(rnorm(n*p),nrow=n,ncol=p)
#' y = X %*% beta + rnorm(n)
#' res = mvsusie(X,y,L=10)
#'
#' @importFrom stats var
#' @importFrom susieR susie_get_cs
#' 
#' @export
mvsusie = function (X, Y, L = 10, prior_variance = 0.2,
                    residual_variance = NULL, prior_weights=NULL,
                    standardize = TRUE,intercept = TRUE, approximate = FALSE,
                    estimate_residual_variance = FALSE,
                    estimate_prior_variance = TRUE,
                    estimate_prior_method = "EM",
                    check_null_threshold = 0, prior_tol = 1e-9,
                    compute_objective = TRUE, s_init = NULL,
                    coverage = 0.95, min_abs_corr = 0.5,
                    compute_univariate_zscore = FALSE,
                    precompute_covariances = FALSE,
                    n_thread = 1, max_iter = 100, tol = 1e-3,
                    verbosity = 2, track_fit = FALSE) {
    
  # Adjust prior weights.
  if (is.null(prior_weights))
    prior_weights = rep(1/ncol(X),ncol(X))
  else
    prior_weights = prior_weights / sum(prior_weights)
  if(inherits(prior_variance, 'MashInitializer')){
    prior_variance = prior_variance$clone(deep=T)
  }
  # adjust prior effects
  is_numeric_prior = !(is.matrix(prior_variance) || inherits(prior_variance, 'MashInitializer'))
  if (!is.null(dim(Y)) && ncol(Y) > 1 && is_numeric_prior) stop("Please specify prior variance for the multivariate response Y")
  if (standardize && !is_numeric_prior) {
    # Scale prior variance
    # https://github.com/stephenslab/mvsusieR/blob/master/inst/prototypes/prior_matrices_scale.ipynb
    sigma = sapply(1:ncol(Y), function(i) sd(Y[,i], na.rm=T))
    n = sapply(1:ncol(Y), function(i) length(which(!is.na(Y[,1]))))
    sigma = sigma / sqrt(n)
    # Make sigma numerically more robust against extreme values
    if (estimate_prior_variance) sigma = sigma / max(sigma)
    if (is.matrix(prior_variance)) prior_variance = scale_covariance(prior_variance, sigma)
    else prior_variance$scale_prior_variance(sigma)
  }
  if (verbosity>1) {
    message("Initializing data object ...")
    message(paste("Dimension of X matrix:", nrow(X), ncol(X)))
    message(paste("Dimension of Y matrix:", nrow(Y), ncol(Y)))
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
  s = mvsusie_core(data, s_init, L, prior_variance, prior_weights,
            estimate_residual_variance, estimate_prior_variance, estimate_prior_method, check_null_threshold,
            precompute_covariances, compute_objective, n_thread, max_iter, tol, prior_tol, track_fit, verbosity)
  # CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
    s$sets = susie_get_cs(s, coverage=coverage, X=X, min_abs_corr=min_abs_corr)
  }
  
  # Report z-scores from univariate regression
  if (compute_univariate_zscore)
    s$z = susieR:::calc_z(X,Y,center = intercept,scale = standardize)
  
  # FIXME: need a better approach
  s$variable_names  = colnames(X)
  s$condition_names = colnames(Y)
  return(s)
}
